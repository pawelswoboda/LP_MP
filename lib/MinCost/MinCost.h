
#ifndef __MINCOST_H__
#define __MINCOST_H__

#include <string.h>
#include <assert.h>
#include <stdlib.h>
#include <vector> // do zrobienia: used for cap. Replace cap by normal array
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <algorithm>
#include <memory>

// if GRAPH_ASSERT is defined then all calls to graph construction functions are assert'ed for correctness
// (e.g. that node_id's are valid id's and edge capacities are non-negative).
//#define GRAPH_ASSERT 

//#define MINCOST_DEBUG

namespace MCF {

// do zrobienia: 
// - we have many parallel edges next to each other each with same capacities. Handle this better (i.e. sort parallel edges), so that pushing flow along shortest path can be possibly done faster.
//   possibly hold such edges in a heap sorted by edge cost. Have two of them for each edge with one containing all edges with some residual capacity left and the other one with reversed residual capacity left
//   for this have in each arc a pointer to left and right element of heap.
// - use double ended dijkstra search to speed up shortest path computation. Reuse? Global source and sink node? Search from all nodes with positive excess simultaneously?
// - improve function names: currently "R" can stand for residual or reduced. Substitute R by correct meaning
// - static_assert for integral type for FlowType and for either integral or floating point type f or CostType in destructor
// - could other priority queues be better than the one used for shortest path computation?
// - replace NULL by nullptr everywhere.
// - support various heaps for Dijkstra via template
// - reorder arcs, such that outgoing arcs are consecutive. Better runtime?

template <typename FlowType, typename CostType> class SSP
{
public:
	typedef int NodeId;
	typedef int EdgeId; // do zrobienia: remove, as there will be no such thing anymore
	typedef int ArcId;

	struct Node;
	struct Arc;

	SSP(int NodeNum, int edgeNumMax, void (*err_function)(const char *) = NULL);
   SSP(const SSP<FlowType, CostType>& other); // copy constructor

	// Destructor
	~SSP();

	void AddNodeExcess(NodeId i, FlowType excess);

	// first call returns 0, second 1, and so on.
	// lower_bound < upper_bound
	// cost can be negative.
	EdgeId AddEdge(NodeId i, NodeId j, FlowType lower_bound, FlowType upper_bound, CostType cost);

	CostType Solve();

	///////////////////////////////////////////////////

   void SetCap(ArcId e, FlowType new_cap);
	FlowType GetRCap(EdgeId e) const;
	void SetRCap(EdgeId e, FlowType new_rcap);
	FlowType GetReverseRCap(EdgeId e) const;
	void SetReverseRCap(EdgeId e, FlowType new_rcap);
	void PushFlow(EdgeId e, FlowType delta);
	void UpdateCost(EdgeId e, CostType delta);

   // functions added by Paul Swoboda //
   FlowType GetFlow(EdgeId e) const
   {
      assert(0 <= e && e < edgeNumMax);
      const FlowType ub = cap[e];
      const FlowType lb = cap[N_arc(arcs[e].sister)];
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

   ArcId StartingArc(NodeId i) const {
      assert(0 <= i && i < nodeNum);
      ArcId idx = std::numeric_limits<ArcId>::max();
      for(Arc* a=nodes[i].firstSaturated; a!=nullptr; a=a->next) {
         idx = std::min(N_arc(a),idx);
      }
      for(Arc* a=nodes[i].firstNonsaturated; a!=nullptr; a=a->next) {
         idx = std::min(N_arc(a),idx);
      }
      return idx;
      
   }
   FlowType NoArcs(NodeId i) const {
      if(i < nodeNum-1) {
         return StartingArc(i+1) - StartingArc(i);
      } else {
         return 2*edgeNum - StartingArc(i);
      }
   }

   FlowType GetCap(const ArcId a) const { return cap[a]; }
   void SetCost(const EdgeId e, const CostType c) { UpdateCost(e, c - GetCost(e)); }
	CostType ShortestPath(const NodeId start_node, const NodeId end_node);
   FlowType GetNodeExcess(const NodeId i) const { return nodes[i].excess; }
   NodeId GetTailNodeId(const ArcId e) const { return N_node(arcs[e].sister->head); }
   NodeId GetHeadNodeId(const ArcId e) const { return N_node(arcs[e].head); }
   CostType GetCost(const EdgeId e) const { return arcs[e].cost; }
   CostType GetReducedCost(const EdgeId e) const { assert(e<2*edgeNum); return arcs[e].GetRCost(); }
   void SetPotential(const NodeId i, const FlowType pi) { nodes[i].pi = pi; }
   CostType GetPotential(const NodeId i) const { return nodes[i].pi; }
   const Node& GetNode(const NodeId i) const { return nodes[i]; }

   FlowType GetUpperBound(const EdgeId e) const { return cap[e].upper; }
   FlowType GetLowerBound(const EdgeId e) const { return cap[e].lower; }

   NodeId N_node(Node* i) const { assert(i - nodes >= 0 && i - nodes < nodeNum); return i - nodes; }
   ArcId N_arc(Arc* a) const { assert(a - arcs >= 0 && a-arcs < 2*edgeNum); return a - arcs; }
   int GetNodeNum() const { return nodeNum; }
   int GetEdgeNum() const { return edgeNum; }
   int GetArcNum() const { return 2*edgeNum; }


   // do zrobienia: make private
   void SortArcs();
   void ExchangeArcs(Arc& a, Arc&b);

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
      // do zrobienia: possibly remove prev and next and always iterate through full list. arcs are stored contiguously. Then no separate saturated and non satured arc list is held
		Node* head;
		Arc*  prev;
		Arc*  next;
		Arc*  sister; // reverse arc

		FlowType	r_cap;		// residual capacity
#ifdef MINCOST_DEBUG
		FlowType	cap_orig; // not needed anymore: cap is an array
#endif
		CostType	cost;
		CostType GetRCost() const { return cost + head->pi - sister->head->pi; }
      Node* tail() const { return sister->head; }
	};

private:

	int		nodeNum, edgeNum, edgeNumMax;
   // make std::unique_ptr out of nodes and arcs
	Node	*nodes;
	Arc		*arcs;
	Node*	firstActive;
	int		counter;
	CostType cost;
   std::unique_ptr<FlowType[]> cap; // used to hold original lower and upper bound, such that primal flow can be recomputed. To zrobienia: templatize this, such that this information is not held unless wanted. 


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
	inline void SSP<FlowType, CostType>::AddNodeExcess(NodeId _i, FlowType excess)
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
	inline typename SSP<FlowType, CostType>::EdgeId SSP<FlowType, CostType>::AddEdge(const NodeId _i, const NodeId _j, const FlowType lb, const FlowType ub, CostType cost)
{
	assert(_i>=0 && _i<nodeNum);
	assert(_j>=0 && _j<nodeNum);
	assert(_i!=_j && edgeNum<edgeNumMax);
	assert(lb < ub);

   cap[2*edgeNum] = ub;
   cap[2*edgeNum+1] = -lb;

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

   a->next = nullptr;
   a->prev = nullptr;
   a_rev->next = nullptr;
   a_rev->prev = nullptr;

   // do zrobienia: delete
   /*
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
   */

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

	//if (a->r_cap > 0 && a->GetRCost() < 0) PushFlow(a, a->r_cap);
	//if (a_rev->r_cap > 0 && a_rev->GetRCost() < 0) PushFlow(a_rev, a_rev->r_cap);

	return edgeNum-1;
}

///////////////////////////////////////
///////////////////////////////////////
///////////////////////////////////////

template <typename FlowType, typename CostType> 
	inline void SSP<FlowType, CostType>::DecreaseRCap(Arc* a, FlowType delta)
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
	inline void SSP<FlowType, CostType>::IncreaseRCap(Arc* a, FlowType delta)
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
void SSP<FlowType, CostType>::SetCap(const ArcId a, const FlowType lb, const FlowType ub)
{
   assert(0<=a && a<2*edgeNum);
   assert(lb < ub);
   Arc* a = arcs[a];
   Arc* arc_rev = a->sister;
   const ArcId a_rev = N_arc(arc_rev);
   const FlowType ub_delta = cap[a] - ub;
   const FlowType lb_delta = cap[a_rev] + lb;
   // note: temporarily bounds may be infeasible. If this is so, the reverse arc will need to be changed as well. Wait until this happens and then set residual capacity etc. correctly
   if(new_cap < -cap[e]) { // bounds infeasible, wait until reverse bound will make it correct again

   } else { // set residual capacity correctly
      SetRCap(e, std::max(0, cap[e] - new_cap));
   }
   cap[a] = ub;
   cap[a_rev] = -lb;
}

template <typename FlowType, typename CostType> 
	inline FlowType SSP<FlowType, CostType>::GetRCap(EdgeId e) const
{
   assert(0<=e && e<2*edgeNum);
	Arc* a = &arcs[e];
	return a->r_cap;
}

template <typename FlowType, typename CostType> 
	inline void SSP<FlowType, CostType>::SetRCap(Arc* a, FlowType new_rcap)
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
	inline void SSP<FlowType, CostType>::SetRCap(EdgeId e, FlowType new_rcap)
{
	SetRCap(&arcs[e], new_rcap);
}

template <typename FlowType, typename CostType> 
	inline FlowType SSP<FlowType, CostType>::GetReverseRCap(EdgeId e) const
{
	Arc* a = arcs[e].sister;
	return a->r_cap;
}

template <typename FlowType, typename CostType> 
	inline void SSP<FlowType, CostType>::SetReverseRCap(EdgeId e, FlowType new_rcap)
{
	SetRCap(arcs[e].sister, new_rcap);
}

template <typename FlowType, typename CostType> 
	inline void SSP<FlowType, CostType>::PushFlow(Arc* a, FlowType delta)
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
	inline void SSP<FlowType, CostType>::PushFlow(EdgeId e, FlowType delta)
{
	PushFlow(&arcs[2*e], delta);
}

template <typename FlowType, typename CostType> 
	inline void SSP<FlowType, CostType>::UpdateCost(EdgeId e, CostType delta)
{
   assert(0<=e && e<edgeNum);
	Arc* a = &arcs[e];
	cost += delta*(cap[e]-a->r_cap);
	a->cost += delta;
	a->sister->cost = -a->cost;

	if (a->GetRCost() > 0) a = a->sister;
	if (a->r_cap > 0 && a->GetRCost() < 0) PushFlow(a, a->r_cap);
}

///////////////////////////////////////
///////////////////////////////////////
///////////////////////////////////////

template <typename FlowType, typename CostType> 
	inline SSP<FlowType, CostType>::PriorityQueue::PriorityQueue()
{
	N = 0;
	arraySize = 16;
	array = (Item*) malloc(arraySize*sizeof(Item));
}

template <typename FlowType, typename CostType> 
	inline SSP<FlowType, CostType>::PriorityQueue::~PriorityQueue()
{
	free(array);
}

template <typename FlowType, typename CostType> 
	inline void SSP<FlowType, CostType>::PriorityQueue::Reset()
{
	N = 0;
}

template <typename FlowType, typename CostType> 
	inline CostType SSP<FlowType, CostType>::PriorityQueue::GetKey(Node* i)
{
	return array[i->heap_ptr].key;
}

template <typename FlowType, typename CostType> 
	inline void SSP<FlowType, CostType>::PriorityQueue::Swap(int k1, int k2)
{
	Item* a = array+k1;
	Item* b = array+k2;
	a->i->heap_ptr = k2;
	b->i->heap_ptr = k1;
	Node* i = a->i;   a->i   = b->i;   b->i   = i;
	CostType key = a->key; a->key = b->key; b->key = key;
}

template <typename FlowType, typename CostType> 
	inline void SSP<FlowType, CostType>::PriorityQueue::Add(Node* i, CostType key)
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
	inline void SSP<FlowType, CostType>::PriorityQueue::DecreaseKey(Node* i, CostType key)
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
	inline typename SSP<FlowType, CostType>::Node* SSP<FlowType, CostType>::PriorityQueue::RemoveMin(CostType& key)
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


template <typename FlowType, typename CostType> 
	inline SSP<FlowType, CostType>::SSP(int _nodeNum, int _edgeNumMax, void (*err_function)(const char *))
	: nodeNum(_nodeNum),
	  edgeNum(0),
	  edgeNumMax(_edgeNumMax),
	  counter(0),
	  cost(0),
	  error_function(err_function)
{
	nodes = (Node*) malloc(nodeNum*sizeof(Node));
	arcs = (Arc*) malloc(2*edgeNumMax*sizeof(Arc));
   cap = std::unique_ptr<FlowType[]>( new FlowType[2*edgeNumMax] );
	if (!nodes || !arcs ||!cap.get()) { if (error_function) (*error_function)("Not enough memory!"); exit(1); }

	memset(nodes, 0, nodeNum*sizeof(Node));
	memset(arcs, 0, 2*edgeNumMax*sizeof(Arc));
	firstActive = &nodes[nodeNum];
#ifdef MINCOST_DEBUG
	for (int i=0; i<nodeNum; i++) nodes[i].id = i;
#endif
}

template <typename FlowType, typename CostType> 
SSP<FlowType, CostType>::SSP(const SSP<FlowType, CostType>& other)
	: nodeNum(other.nodeNum),
	  edgeNum(other.edgeNum),
	  edgeNumMax(other.edgeNumMax),
	  counter(0),
	  cost(0),
	  error_function(other.error_function),
     cap(other.cap)
{
   assert(false); // not working currently
	nodes = (Node*) malloc(nodeNum*sizeof(Node));
	arcs = (Arc*) malloc(2*edgeNumMax*sizeof(Arc));
	if (!nodes || !arcs) { if (error_function) (*error_function)("Not enough memory!"); exit(1); }

	memset(arcs+edgeNum, 0, 2*(edgeNumMax-edgeNum)*sizeof(Arc));

   cap = std::unique_ptr<FlowType[]>( new FlowType[2*edgeNumMax] );
   
	firstActive = &nodes[nodeNum];

	for (int i=0; i<nodeNum; i++) {
      AddNodeExcess(i, other.GetNodeExcess(i));
      SetPotential(i, other.GetPotential(i));
   }
   for(int e=0; e<edgeNum; e++) {
      // AddEdge has been modified.
      exit(1);
      AddEdge(other.GetTailNodeId(e), other.GetHeadNodeId(e), other.GetRCap(e), other.GetReverseRCap(e), other.GetCost(e));
   }
#ifdef MINCOST_DEBUG
	for (int i=0; i<nodeNum; i++) nodes[i].id = i;
#endif
}

template <typename FlowType, typename CostType> 
	SSP<FlowType, CostType>::~SSP()
{
	free(nodes);
	free(arcs);
}

template<typename FlowType, typename CostType>
void SSP<FlowType, CostType>::ExchangeArcs(Arc& a, Arc&b)
{
   std::swap(a, b);
   std::swap(cap[N_arc(&a)], cap[N_arc(&b)]);
   a.sister->sister = &a;
   b.sister->sister = &b;
}

// sort arcs lexicographically.
template<typename FlowType, typename CostType>
void SSP<FlowType, CostType>::SortArcs()
{
   auto perm = std::unique_ptr<ArcId[]>({ new FlowType[2*edgeNum] });
   for(int c=0; c<2*edgeNum; ++c) {
      perm[c] = c;
   }
   std::sort(perm.get(), perm.get()+2*edgeNum, [this](ArcId i, ArcId j) {
         auto tail_i = GetTailNodeId(i);
         auto tail_j = GetTailNodeId(j);
         if(tail_i != tail_j) {
            return tail_i < tail_j;
         } 
         return GetHeadNodeId(i) < GetHeadNodeId(j);
         });
   // follow cycles in permutation. negative permutation entries signify visited indices
   for(int c=0; c<2*edgeNum; ++c) {
      int next_idx = perm[c];
      if(next_idx == c || next_idx < 0) {
         continue;
      }
      int cur_idx = c;
      while(perm[next_idx] >= 0) {
         ExchangeArcs(arcs[cur_idx], arcs[next_idx]);
         perm[cur_idx] -= 2*edgeNum; // mark as visited
         cur_idx = next_idx;
         next_idx = perm[next_idx];
      }
   }
   // set firstSaturated and firstNonSaturated in nodes correctly.
   // set next and prev fields in Arc correctly
   for(Arc* a = arcs; a<arcs+2*edgeNum; ++a) {
      Node* tail = a->tail();
      if(a->r_cap > 0) { // put into non-saturated list
         if(tail->firstNonsaturated) {
            tail->firstNonsaturated->prev = a;
         }
         a->next = tail->firstNonsaturated;
         tail->firstNonsaturated = a;
      } else { // put into saturated list
         if(tail->firstSaturated) {
            tail->firstSaturated->prev = a;
         }
         a->next = tail->firstSaturated;
         a->prev = nullptr;
         tail->firstSaturated = a;
      }
   }
   // push flow
   for(Arc* a = arcs; a<arcs+2*edgeNum; ++a) {
      if (a->r_cap > 0 && a->GetRCost() < 0) PushFlow(a, a->r_cap);
   }
}

template <typename FlowType, typename CostType> 
	void SSP<FlowType, CostType>::Init()
{
	Node* i;
	Arc* a;

	for (a=arcs; a<arcs+2*edgeNum; a++)
	{
		if (a->r_cap > 0 && a->GetRCost() < 0) PushFlow(a, a->r_cap);
	}

	Node** lastActivePtr = &firstActive;
	for (i=nodes; i<nodes+nodeNum; i++)
	{
		if (i->excess > 0)
		{
			*lastActivePtr = i;
			lastActivePtr = &i->next;
		}
		else i->next = NULL;
	}
	*lastActivePtr = &nodes[nodeNum];
}


template <typename FlowType, typename CostType> 
	FlowType SSP<FlowType, CostType>::Augment(Node* start, Node* end)
{
	FlowType delta = (start->excess < -end->excess) ? start->excess : -end->excess;
	Arc* a;

	for (a=end->parent; a; a=a->sister->head->parent)
	{
		if (delta > a->r_cap) delta = a->r_cap;
	}
	assert(delta > 0);

	end->excess += delta;
	for (a=end->parent; a; a=a->head->parent)
	{
		DecreaseRCap(a, delta);
		a = a->sister;
		IncreaseRCap(a, delta);
	}
	start->excess -= delta;

	return delta;
}


// function which computes cost of a shortest path between specified nodes given optimal primal/dual values (i.e. after solve) for computing marginals
// do zrobienia: not tested
template <typename FlowType, typename CostType> 
	CostType SSP<FlowType, CostType>::ShortestPath(const NodeId start_node, const NodeId end_node)
{
   Node* start = &nodes[start_node];
   Node* end = &nodes[end_node];

   int FLAG0 = ++ counter; // permanently labelled
   int FLAG1 = ++ counter; // temporarily labelled

   start->parent = NULL;
   start->flag = FLAG1;
   queue.Reset();
   queue.Add(start, 0.0);

   Node* i;

	CostType d; // the current minimum distance
	while ( (i=queue.RemoveMin(d)) )
   {
      if(i == end) break; // do zrobienia: possibly directly return d - start->pi + end-pi (albo -+ odwrotnie)
      i->flag = FLAG0;

      for(Arc* a=i->firstNonsaturated; a; a=a->next)
      {
         assert(a->r_cap > 0);
         Node* j = a->head;
         if (j->flag == FLAG0) continue;
         CostType reduced_cost = a->GetRCost(); // must use reduced cost, otherwise cost might be negative and then Dijkstra's algorithm would not work
         assert(reduced_cost > -1e-10);
         if (j->flag == FLAG1)
         {
            if (reduced_cost + d >= queue.GetKey(j)) continue;
            queue.DecreaseKey(j, reduced_cost + d);
         }
         else
         {
            queue.Add(j, reduced_cost + d);
            j->flag = FLAG1;
         }
         j->parent = a;
      }
   }
   assert(i == end);

   // trace back to start node via parent pointers and record cost of shortest path
   CostType path_cost = 0.0;
   while(i != start) 
   {
      Arc* a = i->parent;
      // do zrobienia: possibly use reduced cost. but this has to be done everywhere
      //path_cost += a->GetRCost();
      path_cost += a->cost;
      i = a->sister->head; // the tail
   }
   /*
   assert(i->parent == NULL);
   assert(path_cost >= -1e-7);
   if(path_cost < -1e-7) { // kwaskwas
      printf("error: negative path cost = %f\n", path_cost);  
      exit(1);
   }
   assert((d - path_cost) < 1e-7 || (path_cost - d) < 1e-7);
   */

   return path_cost;
}

template <typename FlowType, typename CostType> 
	void SSP<FlowType, CostType>::Dijkstra(Node* start)
{
	assert(start->excess > 0);

	Node* i;
	Node* j;
	Arc* a;
	CostType d;
	Node* permanentNodes;

	int FLAG0 = ++ counter; // permanently labeled nodes
	int FLAG1 = ++ counter; // temporarily labeled nodes

	start->parent = NULL;
	start->flag = FLAG1;
	queue.Reset();
	queue.Add(start, 0);

	permanentNodes = NULL;

	while ( (i=queue.RemoveMin(d)) )
	{
		if (i->excess < 0)
		{
			FlowType delta = Augment(start, i);
			cost += delta*(d - i->pi + start->pi);
			for (i=permanentNodes; i; i=i->next_permanent) i->pi += d;
			break;
		}

		i->pi -= d;
		i->flag = FLAG0;
		i->next_permanent = permanentNodes;
		permanentNodes = i;

		for (a=i->firstNonsaturated; a; a=a->next)
		{
			j = a->head;
			if (j->flag == FLAG0) continue;
			d = a->GetRCost();
			if (j->flag == FLAG1)
			{
				if (d >= queue.GetKey(j)) continue;
				queue.DecreaseKey(j, d);
			}
			else
			{
				queue.Add(j, d);
				j->flag = FLAG1;
			}
			j->parent = a;
		}

	}
}


template <typename FlowType, typename CostType> 
	CostType SSP<FlowType, CostType>::Solve()
{
	Node* i;
	//Init();
	while ( 1 )
	{
		i = firstActive;
		if (i == &nodes[nodeNum]) break;
		firstActive = i->next;
		i->next = NULL;
		if (i->excess > 0)
		{
			Dijkstra(i);
			if (i->excess > 0 && !i->next) 
			{ 
            assert(i != firstActive);
				i->next = firstActive; 
				firstActive = i; 
			}
		}
	}
#ifdef MINCOST_DEBUG
	TestOptimality();
	TestCosts();
#endif

	return cost;
}


template <typename FlowType, typename CostType> 
	void SSP<FlowType, CostType>::TestOptimality()
{
	Node* i;
	Arc* a;

	for (i=nodes; i<nodes+nodeNum; i++)
	{
		if (i->excess != 0)
		{
			assert(0);
		}
		for (a=i->firstSaturated; a; a=a->next)
		{
			if (a->r_cap != 0)
			{
				assert(0);
			}
		}
		for (a=i->firstNonsaturated; a; a=a->next)
		{
			CostType c = a->GetRCost();
			if (a->r_cap <= 0 || a->GetRCost() < -1e-5)
			{
				assert(0);
			}
		}
	}
}

#ifdef MINCOST_DEBUG

template <typename FlowType, typename CostType> 
	void SSP<FlowType, CostType>::TestCosts()
{
	Arc* a;

	CostType _cost = 0;

   // do zrobienia: this will not work anymore due to different arc ordering
	for (a=arcs; a<arcs+2*edgeNum; a+=2)
	{
		assert(a->r_cap + a->sister->r_cap == a->cap_orig + a->sister->cap_orig);
		_cost += a->cost*(a->cap_orig - a->r_cap);
	}

	CostType delta = cost - _cost;
	if (delta < 0) delta = -delta;
	if (delta >= 1e-5)
	{
		assert(0);
	}
}

#endif

} // end namespace MCF

#endif
