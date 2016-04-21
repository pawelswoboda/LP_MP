
#ifndef __MINCOST_H__
#define __MINCOST_H__

#include <string.h>
#include <assert.h>
#include <stdlib.h>
#include <vector> // do zrobienia: used for bounds
#include <iostream>

// if GRAPH_ASSERT is defined then all calls to graph construction functions are assert'ed for correctness
// (e.g. that node_id's are valid id's and edge capacities are non-negative).
//#define GRAPH_ASSERT 

#define MINCOST_DEBUG


// do zrobienia: 
// - put MinCost inside appropriate namespace
// - we have many parallel edges next to each other each with same capacities. Handle this better (i.e. sort parallel edges), so that pushing flow along shortest path can be possibly done faster.
//   possibly hold such edges in a heap sorted by edge cost. Have two of them for each edge with one containing all edges with some residual capacity left and the other one with reversed residual capacity left
//   for this have in each arc a pointer to left and right element of heap.
// - use double ended dijkstra search to speed up shortest path computation. Reuse? Global source and sink node? Search from all nodes with positive excess simultaneously?
// - improve function names: currently "R" can stand for residual or reduced. Substitute R by correct meaning
// - static_assert for integral type for FlowType and for either integral or floating point type f or CostType in destructor
// - arcs and reverse arcs are stored next to each other in the arcs array. one can get sister arc either by explicitly storing its address or by computing its address, given that we know that adresses are aligned in some way. Is getting away with sister pointer advantageous for runtime? 
// - also it is possible to do some bittricks (as in Knuth) to make some fast modulo arithmetic to get sister arc. Possibly look into "Hacker's Delight".
// - cpp file seems not worth the issue. Templatize code and put everything in one hxx.
// - could other priority queues be better than the one used for shortest path computation?
// - replace NULL by nullptr everywhere.
// - make ArcIterator compatible with standard iterators. Still needed? What about nonSaturatedIterator?
// - support various heaps for Dijkstra via template
// - reorder arcs, such that outgoing arcs are consecutive. Better runtime?

template <typename FlowType, typename CostType> class MinCost
{
public:
	typedef int NodeId;
	typedef int EdgeId;

	struct Node;
	struct Arc;

	MinCost(int NodeNum, int edgeNumMax, void (*err_function)(const char *) = NULL);
   MinCost(const MinCost<FlowType, CostType>& other); // copy constructor

	// Destructor
	~MinCost();

	void AddNodeExcess(NodeId i, FlowType excess);

	// first call returns 0, second 1, and so on.
	// lower_bound < upper_bound
	// cost can be negative.
	EdgeId AddEdge(NodeId i, NodeId j, FlowType lower_bound, FlowType upper_bound, CostType cost);

	CostType Solve();

	///////////////////////////////////////////////////

	FlowType GetRCap(EdgeId e) const;
	void SetRCap(EdgeId e, FlowType new_rcap);
	FlowType GetReverseRCap(EdgeId e) const;
	void SetReverseRCap(EdgeId e, FlowType new_rcap);
	void PushFlow(EdgeId e, FlowType delta);
	void UpdateCost(EdgeId e, FlowType cap_orig, CostType delta);

   // functions added by Paul Swoboda //
   FlowType GetFlow(EdgeId e) const
   {
      assert(0 <= e && e < edgeNumMax);
      const FlowType lb = bounds[e].lower;
      const FlowType ub = bounds[e].upper;
      const FlowType r_cap = GetRCap(e);
      //std::cout << "ub = " << ub << ", residual capacity = " << r_cap << ", lb = " << lb << ", reverse residual capacity = " << reverse_r_cap << std::endl;
      FlowType delta;
      if(ub < 0) {
         delta += ub;
      } else if(lb > 0) {
         delta -= lb;
      } else {
         delta = 0;
      }
      return ub + delta - r_cap;
   }
   FlowType ExcessSum() const
   {
      FlowType sum = 0;
      for(int i=0; i<nodeNum; ++i) {
         sum += nodes[i].excess;
      }
      return sum;
   }

   void SetCost(const EdgeId e, const FlowType cap_orig, const CostType c) { UpdateCost(e, cap_orig, c - GetCost(e)); }
	CostType ShortestPath(const NodeId start_node, const NodeId end_node);
   FlowType GetNodeExcess(const NodeId i) const { return nodes[i].excess; }
   NodeId GetTailNodeId(const EdgeId e) const { assert(&arcs[2*e+1] == arcs[2*e].sister); return arcs[2*e+1].head - nodes; }
   NodeId GetHeadNodeId(const EdgeId e) const { return arcs[2*e].head - nodes; }
   CostType GetCost(const EdgeId e) const { return arcs[2*e].cost; }
   CostType GetReducedCost(const EdgeId e) const { assert(e<edgeNum); return arcs[2*e].GetRCost(); }
   void SetPotential(const NodeId i, const FlowType pi) { nodes[i].pi = pi; }
   CostType GetPotential(const NodeId i) const { return nodes[i].pi; }
   const Node& GetNode(const NodeId i) const { return nodes[i]; }
   EdgeId GetEdgeId(const Arc* const arc) const { assert(arc >= arcs && arc - arcs < 2*edgeNum); return (arc - arcs)/2; }

   FlowType GetUpperBound(const EdgeId e) const { return bounds[e].upper; }
   FlowType GetLowerBound(const EdgeId e) const { return bounds[e].lower; }

   int GetNodeNum() const { return nodeNum; }
   int GetEdgeNum() const { return edgeNum; }

   // delete all edges and reset all nodes
   void Reset() 
   {
	  edgeNum = 0;
	  counter = 0;
	  cost = 0;
     memset(nodes, 0, nodeNum*sizeof(Node));
     memset(arcs, 0, 2*edgeNumMax*sizeof(Arc));
     firstActive = &nodes[nodeNum];
#ifdef MINCOST_DEBUG
     for (int i=0; i<nodeNum; i++) nodes[i].id = i;
#endif
   }

/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
	
	// internal variables and functions

   struct Node
   {
      // I doubt that this structure is of much use, as writing any costs to the arcs is not possible while iterating, due to writing changing the linking structure. Possibly remove iterator again.
      struct ArcIterator // iterate over all arcs emanating from some node. Note that edges may be given in arbitrary order, and this order may change
      {
         ArcIterator(Arc* const arc) : arc_(arc) {}
         Arc* arc_;
         bool inNonsaturatedList_ = true;

         ArcIterator& operator++() 
         {
            assert(arc_ != NULL); 
            if(arc_->next == NULL && inNonsaturatedList_ == true) { arc_ = arc_->sister->head->firstSaturated; inNonsaturatedList_ = false; }
            else { assert(arc_->next != NULL); arc_ = arc_->next; }
            return *this; 
         }
         bool operator==(const ArcIterator& other) { return other.arc_ == arc_; }
         bool operator!=(const ArcIterator& other) { return !operator==(other); } // do zrobienia: direct implementation faster?
         const Arc* operator*() { return arc_; } // note that it is dangerous to write anything to arc, as this will change linking structure of arcs.
         const Arc* operator->() { return arc_; }

      };
      ArcIterator begin() const { return ArcIterator(this->firstNonsaturated); }
      ArcIterator end() const { return ArcIterator(NULL); }

		Arc			*firstNonsaturated;
		Arc			*firstSaturated;

		Arc			*parent;
		Node		*next; // list of nodes with positive excesses

		FlowType	excess;
		CostType	pi;
		int			flag;
		union
		{
			int		heap_ptr;
			Node*	next_permanent;
		};
#ifdef MINCOST_DEBUG
		int			id;
#endif
	};

	struct Arc
	{
		Node		*head;
		Arc			*prev;
		Arc			*next;
		Arc			*sister;	// reverse arc

		FlowType	r_cap;		// residual capacity
#ifdef MINCOST_DEBUG
		FlowType	cap_orig;
#endif
		CostType	cost;
		CostType GetRCost() const { return cost + head->pi - sister->head->pi; }
	};

private:

	int		nodeNum, edgeNum, edgeNumMax;
	Node	*nodes;
	Arc		*arcs;
	Node*	firstActive;
	int		counter;
	CostType cost;
   struct bound { FlowType lower, upper; };
   std::vector<bound> bounds; // used to hold original lower and upper bound, such that primal flow can be recomputed. To zrobienia: templatize this, such that this information is not held unless wanted. Default: false


	void	(*error_function)(const char *);	// this function is called if a error occurs,
										// with a corresponding error message
										// (or exit(1) is called if it's NULL)

	/////////////////////////////////////////////////////////////////////////

	struct PriorityQueue
	{
		PriorityQueue();
		~PriorityQueue();
		void Reset();
		CostType GetKey(Node* i);
		void Add(Node* i, CostType key);
		void DecreaseKey(Node* i, CostType key);
		Node* RemoveMin(CostType& key);

	private:
		struct Item
		{
			Node*		i;
			CostType	key;
		}* array;
		int N, arraySize;
		void Swap(int k1, int k2);
	};

	PriorityQueue queue;

	/////////////////////////////////////////////////////////////////////////

	void SetRCap(Arc* a, FlowType new_rcap);
	void PushFlow(Arc* a, FlowType delta);

	void Init();
	void DecreaseRCap(Arc* a, FlowType delta);
	void IncreaseRCap(Arc* a, FlowType delta);
	FlowType Augment(Node* start, Node* end);
	void Dijkstra(Node* start);

	void TestOptimality();
#ifdef MINCOST_DEBUG
	void TestCosts();
#endif
};











///////////////////////////////////////
// Implementation - inline functions //
///////////////////////////////////////



template <typename FlowType, typename CostType> 
	inline void MinCost<FlowType, CostType>::AddNodeExcess(NodeId _i, FlowType excess)
{
	assert(_i>=0 && _i<nodeNum);
	nodes[_i].excess += excess;
	if (nodes[_i].excess > 0 && !nodes[_i].next)
	{
		nodes[_i].next = firstActive;
		firstActive = &nodes[_i];
	}
}

// alternative definition in terms of lower and upper bound
template <typename FlowType, typename CostType> 
	inline typename MinCost<FlowType, CostType>::EdgeId MinCost<FlowType, CostType>::AddEdge(const NodeId _i, const NodeId _j, const FlowType lb, const FlowType ub, CostType cost)
{
	assert(_i>=0 && _i<nodeNum);
	assert(_j>=0 && _j<nodeNum);
	assert(_i!=_j && edgeNum<edgeNumMax);
	assert(lb < ub);

   bounds[edgeNum].lower = lb;
   bounds[edgeNum].upper = ub;

   FlowType cap = ub;
   FlowType rev_cap = -lb;

	Arc *a = &arcs[2*edgeNum];
	Arc *a_rev = a+1;
	edgeNum ++;

   // if either lb > 0 or ub < 0, shift both by the minimal amount so that the previous constraints are satisfied and modify excess on either end node
   FlowType delta = 0;
   if(ub < 0) {
      delta = -ub;
   } else if(lb > 0) {
      delta = -lb;
   }
   if(delta != 0) {
      cap += delta;
      rev_cap -= delta;
      AddNodeExcess(_j, +delta);
      AddNodeExcess(_i, -delta);
      cost += -delta*cost;
   }
   assert(cap >= 0 && rev_cap >= 0);

	Node* i = nodes + _i;
	Node* j = nodes + _j;

	a -> sister = a_rev;
	a_rev -> sister = a;

	if (cap > 0)
	{
		if (i->firstNonsaturated) i->firstNonsaturated->prev = a;
		a -> next = i -> firstNonsaturated;
		i -> firstNonsaturated = a;
	}
	else
	{
		if (i->firstSaturated) i->firstSaturated->prev = a;
		a -> next = i -> firstSaturated;
		i -> firstSaturated = a;
	}
	a->prev = NULL;
	if (rev_cap > 0)
	{
		if (j->firstNonsaturated) j->firstNonsaturated->prev = a_rev;
		a_rev -> next = j -> firstNonsaturated;
		j -> firstNonsaturated = a_rev;
	}
	else
	{
		if (j->firstSaturated) j->firstSaturated->prev = a_rev;
		a_rev -> next = j -> firstSaturated;
		j -> firstSaturated = a_rev;
	}
	a_rev->prev = NULL;

	a -> head = j;
	a_rev -> head = i;
	a -> r_cap = cap;
	a_rev -> r_cap = rev_cap;
	a -> cost = cost;
	a_rev -> cost = -cost;
#ifdef MINCOST_DEBUG
	a->cap_orig = cap;
	a_rev->cap_orig = rev_cap;
#endif

	if (a->r_cap > 0 && a->GetRCost() < 0) PushFlow(a, a->r_cap);
	if (a_rev->r_cap > 0 && a_rev->GetRCost() < 0) PushFlow(a_rev, a_rev->r_cap);

	return edgeNum-1;
}

///////////////////////////////////////
///////////////////////////////////////
///////////////////////////////////////

template <typename FlowType, typename CostType> 
	inline void MinCost<FlowType, CostType>::DecreaseRCap(Arc* a, FlowType delta)
{
	a->r_cap -= delta;
	if (a->r_cap == 0)
	{
		Node* i = a->sister->head;
		if (a->next) a->next->prev = a->prev;
		if (a->prev) a->prev->next = a->next;
		else         i->firstNonsaturated = a->next;
		a->next = i->firstSaturated;
		if (a->next) a->next->prev = a;
		a->prev = NULL;
		i->firstSaturated = a;
	}
}

template <typename FlowType, typename CostType> 
	inline void MinCost<FlowType, CostType>::IncreaseRCap(Arc* a, FlowType delta)
{
	if (a->r_cap == 0)
	{
		Node* i = a->sister->head;
		if (a->next) a->next->prev = a->prev;
		if (a->prev) a->prev->next = a->next;
		else         i->firstSaturated = a->next;
		a->next = i->firstNonsaturated;
		if (a->next) a->next->prev = a;
		a->prev = NULL;
		i->firstNonsaturated = a;
	}
	a->r_cap += delta;
}

template <typename FlowType, typename CostType> 
	inline FlowType MinCost<FlowType, CostType>::GetRCap(EdgeId e) const
{
   assert(0<=e && e<edgeNum);
	Arc* a = &arcs[2*e];
	return a->r_cap;
}

template <typename FlowType, typename CostType> 
	inline void MinCost<FlowType, CostType>::SetRCap(Arc* a, FlowType new_rcap)
{
	assert(new_rcap >= 0);
#ifdef MINCOST_DEBUG
	a->cap_orig += new_rcap - a->r_cap;
#endif
	if (a->r_cap == 0)
	{
		Node* i = a->sister->head;
		if (a->next) a->next->prev = a->prev;
		if (a->prev) a->prev->next = a->next;
		else         i->firstSaturated = a->next;
		a->next = i->firstNonsaturated;
		if (a->next) a->next->prev = a;
		a->prev = NULL;
		i->firstNonsaturated = a;
	}
	a->r_cap = new_rcap;
	if (a->r_cap == 0)
	{
		Node* i = a->sister->head;
		if (a->next) a->next->prev = a->prev;
		if (a->prev) a->prev->next = a->next;
		else         i->firstNonsaturated = a->next;
		a->next = i->firstSaturated;
		if (a->next) a->next->prev = a;
		a->prev = NULL;
		i->firstSaturated = a;
	}
}

template <typename FlowType, typename CostType> 
	inline void MinCost<FlowType, CostType>::SetRCap(EdgeId e, FlowType new_rcap)
{
	SetRCap(&arcs[2*e], new_rcap);
}

template <typename FlowType, typename CostType> 
	inline FlowType MinCost<FlowType, CostType>::GetReverseRCap(EdgeId e) const
{
	Arc* a = &arcs[2*e+1];
	return a->r_cap;
}

template <typename FlowType, typename CostType> 
	inline void MinCost<FlowType, CostType>::SetReverseRCap(EdgeId e, FlowType new_rcap)
{
	SetRCap(&arcs[2*e+1], new_rcap);
}

template <typename FlowType, typename CostType> 
	inline void MinCost<FlowType, CostType>::PushFlow(Arc* a, FlowType delta)
{
	if (delta < 0) { a = a->sister; delta = -delta; }
	DecreaseRCap(a, delta);
	IncreaseRCap(a->sister, delta);
	a->head->excess += delta;
	a->sister->head->excess -= delta;
	cost += delta*a->cost;
	if (a->head->excess > 0 && !a->head->next)
	{
		a->head->next = firstActive;
		firstActive = a->head;
	}
}

template <typename FlowType, typename CostType> 
	inline void MinCost<FlowType, CostType>::PushFlow(EdgeId e, FlowType delta)
{
	PushFlow(&arcs[2*e], delta);
}

template <typename FlowType, typename CostType> 
	inline void MinCost<FlowType, CostType>::UpdateCost(EdgeId e, FlowType cap_orig, CostType delta)
{
   assert(0<=e && e<edgeNum);
	Arc* a = &arcs[2*e];
	cost += delta*(cap_orig-a->r_cap);
	a->cost += delta;
	a->sister->cost = -a->cost;

	if (a->GetRCost() > 0) a = a->sister;
	if (a->r_cap > 0 && a->GetRCost() < 0) PushFlow(a, a->r_cap);
}

///////////////////////////////////////
///////////////////////////////////////
///////////////////////////////////////

template <typename FlowType, typename CostType> 
	inline MinCost<FlowType, CostType>::PriorityQueue::PriorityQueue()
{
	N = 0;
	arraySize = 16;
	array = (Item*) malloc(arraySize*sizeof(Item));
}

template <typename FlowType, typename CostType> 
	inline MinCost<FlowType, CostType>::PriorityQueue::~PriorityQueue()
{
	free(array);
}

template <typename FlowType, typename CostType> 
	inline void MinCost<FlowType, CostType>::PriorityQueue::Reset()
{
	N = 0;
}

template <typename FlowType, typename CostType> 
	inline CostType MinCost<FlowType, CostType>::PriorityQueue::GetKey(Node* i)
{
	return array[i->heap_ptr].key;
}

template <typename FlowType, typename CostType> 
	inline void MinCost<FlowType, CostType>::PriorityQueue::Swap(int k1, int k2)
{
	Item* a = array+k1;
	Item* b = array+k2;
	a->i->heap_ptr = k2;
	b->i->heap_ptr = k1;
	Node* i = a->i;   a->i   = b->i;   b->i   = i;
	CostType key = a->key; a->key = b->key; b->key = key;
}

template <typename FlowType, typename CostType> 
	inline void MinCost<FlowType, CostType>::PriorityQueue::Add(Node* i, CostType key)
{
	if (N == arraySize)
	{
		arraySize *= 2;
		array = (Item*) realloc(array, arraySize*sizeof(Item));
	}
	int k = i->heap_ptr = N ++;
	array[k].i = i;
	array[k].key = key;
	while (k > 0)
	{
		int k2 = (k-1)/2;
		if (array[k2].key <= array[k].key) break;
		Swap(k, k2);
		k = k2;
	}
}

template <typename FlowType, typename CostType> 
	inline void MinCost<FlowType, CostType>::PriorityQueue::DecreaseKey(Node* i, CostType key)
{
	int k = i->heap_ptr;
	array[k].key = key;
	while (k > 0)
	{
		int k2 = (k-1)/2;
		if (array[k2].key <= array[k].key) break;
		Swap(k, k2);
		k = k2;
	}
}

template <typename FlowType, typename CostType> 
	inline typename MinCost<FlowType, CostType>::Node* MinCost<FlowType, CostType>::PriorityQueue::RemoveMin(CostType& key)
{
	if (N == 0) return NULL;

	Swap(0, N-1);
	N --;

	int k = 0;
	while ( 1 )
	{
		int k1 = 2*k + 1, k2 = k1 + 1;
		if (k1 >= N) break;
		int k_min = (k2 >= N || array[k1].key <= array[k2].key) ? k1 : k2;
		if (array[k].key <= array[k_min].key) break;
		Swap(k, k_min);
		k = k_min;
	}

	key = array[N].key;
	return array[N].i;
}

#endif
