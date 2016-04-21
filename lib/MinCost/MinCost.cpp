#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "MinCost.h"


template <typename FlowType, typename CostType> 
	inline MinCost<FlowType, CostType>::MinCost(int _nodeNum, int _edgeNumMax, void (*err_function)(const char *))
	: nodeNum(_nodeNum),
	  edgeNum(0),
	  edgeNumMax(_edgeNumMax),
	  counter(0),
	  cost(0),
	  error_function(err_function),
     bounds(_edgeNumMax)
{
	nodes = (Node*) malloc(nodeNum*sizeof(Node));
	arcs = (Arc*) malloc(2*edgeNumMax*sizeof(Arc));
	if (!nodes || !arcs) { if (error_function) (*error_function)("Not enough memory!"); exit(1); }

	memset(nodes, 0, nodeNum*sizeof(Node));
	memset(arcs, 0, 2*edgeNumMax*sizeof(Arc));
	firstActive = &nodes[nodeNum];
#ifdef MINCOST_DEBUG
	for (int i=0; i<nodeNum; i++) nodes[i].id = i;
#endif
}

template <typename FlowType, typename CostType> 
MinCost<FlowType, CostType>::MinCost(const MinCost<FlowType, CostType>& other)
	: nodeNum(other.nodeNum),
	  edgeNum(other.edgeNum),
	  edgeNumMax(other.edgeNumMax),
	  counter(0),
	  cost(0),
	  error_function(other.error_function),
     bounds(other.bounds)
{
	nodes = (Node*) malloc(nodeNum*sizeof(Node));
	arcs = (Arc*) malloc(2*edgeNumMax*sizeof(Arc));
	if (!nodes || !arcs) { if (error_function) (*error_function)("Not enough memory!"); exit(1); }

	memset(arcs+edgeNum, 0, 2*(edgeNumMax-edgeNum)*sizeof(Arc));
   
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


// function which computes cost of a shortest path between specified nodes given optimal primal/dual values (i.e. after solve) for computing marginals
// do zrobienia: not tested
template <typename FlowType, typename CostType> 
	CostType MinCost<FlowType, CostType>::ShortestPath(const NodeId start_node, const NodeId end_node)
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

#include "config.hxx"
