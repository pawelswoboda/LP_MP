
#ifndef __MINCOST_H__
#define __MINCOST_H__

#include <string.h>
#include <assert.h>

#include <algorithm>
#include <memory>

// if GRAPH_ASSERT is defined then all calls to graph construction functions are assert'ed for correctness
// (e.g. that node_id's are valid id's and edge capacities are non-negative).
//#define GRAPH_ASSERT 

//#define MINCOST_DEBUG

template <typename FlowType, typename CostType> class MinCost
{
public:
	typedef int NodeId;
	typedef int EdgeId;

	MinCost(int NodeNum, int edgeNumMax, void (*err_function)(const char *) = NULL);

	// Destructor
	~MinCost();

	void AddNodeExcess(NodeId i, FlowType excess);

	// first call returns 0, second 1, and so on.
	// cap, rev_cap must be non-negative. 
	// cost can be negative.
	EdgeId AddEdge(NodeId i, NodeId j, FlowType cap, FlowType rev_cap, CostType cost);

	CostType Solve();

	///////////////////////////////////////////////////

	FlowType GetRCap(EdgeId e);
	void SetRCap(EdgeId e, FlowType new_rcap);
	FlowType GetReverseRCap(EdgeId e);
	void SetReverseRCap(EdgeId e, FlowType new_rcap);
	void PushFlow(EdgeId e, FlowType delta);
	void UpdateCost(EdgeId e, CostType delta);


/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
	
private:
	// internal variables and functions

	struct Node;
	struct Arc;

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
		Node		*head;
		Arc			*prev;
		Arc			*next;
		Arc			*sister;	// reverse arc

		FlowType	r_cap;		// residual capacity
#ifdef MINCOST_DEBUG
		FlowType	cap_orig;
#endif
		CostType	cost;
		CostType GetRCost() { return cost + head->pi - sister->head->pi; }
	};

	int		nodeNum, edgeNum, edgeNumMax;
	Node	*nodes;
	Arc		*arcs;
  std::unique_ptr<FlowType[]> cap; // used to hold original lower and upper bound, such that primal flow can be recomputed. To zrobienia: templatize this, such that this information is not held unless wanted. 
	Node*	firstActive;
	int		counter;
	CostType cost;


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

public:
   // functions added by Paul Swoboda //
   FlowType GetFlow(EdgeId e) const
   {
      assert(0 <= e && e < edgeNumMax);
      return cap[e] - arcs[2*e].r_cap;
   }
   FlowType GetCap(EdgeId e) const
   {
      assert(0 <= e && e < edgeNumMax);
      return cap[2*e];
   }
   FlowType GetRCap(EdgeId e) const
   {
      assert(0 <= e && e < edgeNumMax);
      return cap[2*e+1];
   }

   void SetCost(const EdgeId e, const CostType c);
   NodeId GetTailNodeId(const EdgeId e) const { return N_node(arcs[2*e].sister->head); }
   NodeId GetHeadNodeId(const EdgeId e) const { return N_node(arcs[2*e].head); }
   EdgeId GetReverseEdgeId(const EdgeId e) const { return N_arc(arcs[2*e].sister); }
   CostType GetCost(const EdgeId e) const { return arcs[2*e].cost; }
   CostType GetReducedCost(const EdgeId e) const { assert(e<edgeNum); return arcs[2*e].GetRCost(); }
   NodeId N_node(Node* i) const { assert(i - nodes >= 0 && i - nodes < nodeNum); return i - nodes; }
   EdgeId N_arc(Arc* a) const { assert(a - arcs >= 0 && a-arcs < 2*edgeNum); return a - arcs; }
   int GetNodeNum() const { return nodeNum; }
   int GetEdgeNum() const { return edgeNum; }

   CostType Objective() const
   {
      CostType c = 0.0;
      NodeId a;
      Arc* arc;
      for(arc=arcs, a=0; a<edgeNum; ++arc, ++a) {
         c += GetFlow(2*a)*(arc->cost);
      }
      return 0.5*c; // we have counted objective for arcs and reverse arcs
   }
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

template <typename FlowType, typename CostType> 
	inline typename MinCost<FlowType, CostType>::EdgeId MinCost<FlowType, CostType>::AddEdge(NodeId _i, NodeId _j, FlowType cap_orig, FlowType rev_cap, CostType cost)
{
	assert(_i>=0 && _i<nodeNum);
	assert(_j>=0 && _j<nodeNum);
	assert(_i!=_j && edgeNum<edgeNumMax);
	assert(cap >= 0);
	assert(rev_cap >= 0);

	Arc *a = &arcs[2*edgeNum];
	Arc *a_rev = a+1;
  cap[2*edgeNum] = cap_orig;
  cap[2*edgeNum+1] = rev_cap;
	edgeNum ++;

	Node* i = nodes + _i;
	Node* j = nodes + _j;

	a -> sister = a_rev;
	a_rev -> sister = a;
	if (cap_orig > 0)
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
	a -> r_cap = cap_orig;
	a_rev -> r_cap = rev_cap;
	a -> cost = cost;
	a_rev -> cost = -cost;
#ifdef MINCOST_DEBUG
	a->cap_orig = cap_orig;
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
	inline FlowType MinCost<FlowType, CostType>::GetRCap(EdgeId e)
{
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
	inline FlowType MinCost<FlowType, CostType>::GetReverseRCap(EdgeId e)
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
inline void MinCost<FlowType, CostType>::SetCost(EdgeId e, CostType c)
{
   assert(0<=e && e<edgeNum);
   Arc* arc = &arcs[2*e];
   arc->cost = c;
   arc->sister->cost = -c;

	if (arc->GetRCost() > 0) arc = arc->sister;
	if (arc->r_cap > 0 && arc->GetRCost() < 0) PushFlow(arc, arc->r_cap);
}

template <typename FlowType, typename CostType> 
	inline void MinCost<FlowType, CostType>::UpdateCost(EdgeId e, CostType delta)
{
	Arc* a = &arcs[2*e];
	cost += delta*(cap[2*e]-a->r_cap);
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


template <typename FlowType, typename CostType> 
	inline MinCost<FlowType, CostType>::MinCost(int _nodeNum, int _edgeNumMax, void (*err_function)(const char *))
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
	if (!nodes || !arcs || !cap) { if (error_function) (*error_function)("Not enough memory!"); exit(1); }

	memset(nodes, 0, nodeNum*sizeof(Node));
	memset(arcs, 0, 2*edgeNumMax*sizeof(Arc));
	firstActive = &nodes[nodeNum];
#ifdef MINCOST_DEBUG
	for (int i=0; i<nodeNum; i++) nodes[i].id = i;
#endif
}

template <typename FlowType, typename CostType> 
	MinCost<FlowType, CostType>::~MinCost()
{
	free(nodes);
	free(arcs);
}

template <typename FlowType, typename CostType> 
	void MinCost<FlowType, CostType>::Init()
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
	FlowType MinCost<FlowType, CostType>::Augment(Node* start, Node* end)
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

template <typename FlowType, typename CostType> 
	void MinCost<FlowType, CostType>::Dijkstra(Node* start)
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
	CostType MinCost<FlowType, CostType>::Solve()
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
	void MinCost<FlowType, CostType>::TestOptimality()
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
	void MinCost<FlowType, CostType>::TestCosts()
{
	Arc* a;

	CostType _cost = 0;

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

#endif
