#ifndef LP_MP_MULTICUT_CONSTRUCTOR_HXX
#define LP_MP_MULTICUT_CONSTRUCTOR_HXX

#include "LP_MP.h"
#include "multicut_unary_factor.hxx"
#include "multicut_triplet_factor.hxx"
#include "multicut_odd_wheel.hxx"
#include "lifted_multicut_factors_messages.hxx"

#include "union_find.hxx"
#include "max_flow.hxx"

#include <unordered_map>
#include <unordered_set>
#include <queue>
#include <list>

#include "andres/graph/graph.hxx"
#include "andres/graph/grid-graph.hxx"
#include "andres/graph/multicut/kernighan-lin.hxx"
#include "andres/graph/multicut/greedy-additive.hxx"
#include "andres/graph/multicut-lifted/kernighan-lin.hxx"
#include "andres/graph/multicut-lifted/greedy-additive.hxx"


namespace LP_MP {


// hash function for maps used in constructors. Do zrobienia: define hash functions used somewhere globally in config.hxx

template<class FACTOR_MESSAGE_CONNECTION, INDEX UNARY_FACTOR_NO, INDEX TRIPLET_FACTOR_NO, INDEX UNARY_TRIPLET_MESSAGE_NO, INDEX CONSTANT_FACTOR_NO>
class MulticutConstructor {
public:
   using FMC = FACTOR_MESSAGE_CONNECTION;

   using UnaryFactorContainer = meta::at_c<typename FMC::FactorList, UNARY_FACTOR_NO>;
   using TripletFactorContainer = meta::at_c<typename FMC::FactorList, TRIPLET_FACTOR_NO>;
   using ConstantFactorContainer = meta::at_c<typename FMC::FactorList, CONSTANT_FACTOR_NO>;
   //using GlobalFactorContainer = meta::at_c<typename FMC::FactorList, GLOBAL_FACTOR_NO>;
   using UnaryTripletMessageContainer = typename meta::at_c<typename FMC::MessageList, UNARY_TRIPLET_MESSAGE_NO>::MessageContainerType;
   using UnaryTripletMessageType = typename UnaryTripletMessageContainer::MessageType;
   //using UnaryGlobalMessageContainer = typename meta::at_c<typename FMC::MessageList, UNARY_GLOBAL_MESSAGE_NO>::MessageContainerType;
   //using UnaryGlobalMessageType = typename UnaryGlobalMessageContainer::MessageType;

   // efficient graph structure for shortest path computation.
   // arcs are held contiguously w.r.t. tail node and descending w.r.t. cost.
   // Hence, iterating over all outgoing arcs from a given node is fast. 
   // Also early stopping, when we go over a thresholded graph is possible, as arcs are ordered w.r.t. cost when tail node is fixed.
   class Graph2 {
      struct arc;
      struct node {
         arc* first = nullptr;
         arc* last = nullptr;

         arc* begin() const { return first; }
         arc* end() const { return last; }

         void push_back(arc& a) {
            assert(&a == last);
            last += 1;
         }
            
      };
      struct arc {
         REAL cost;
         node* head;
      };

      public:
      template<typename VEC>
      Graph2(const INDEX no_nodes, const INDEX no_arcs, const VEC& number_outgoing_arcs)
      : nodes_(no_nodes,{nullptr,nullptr}),
         arcs_(no_arcs) 
      {
         assert(no_nodes > 0);
         nodes_[0].first = &arcs_[0];
         nodes_[0].last = nodes_[0].first;
         for(INDEX i=1; i<no_nodes; ++i) {
            nodes_[i].first = nodes_[i-1].first + number_outgoing_arcs[i-1];
            nodes_[i].last = nodes_[i].first;
         }
      }
      Graph2(Graph2& o) = delete;

      INDEX size() const { return nodes_.size(); }
      void add_arc(INDEX i, INDEX j, REAL cost) {
         assert(i<nodes_.size());
         assert(j<nodes_.size());
         *nodes_[i].last = arc({cost, &nodes_[j]});
         if(i+1<nodes_.size()) {
            assert(nodes_[i].last < nodes_[i+1].first);
         }
         nodes_[i].last += 1;
      }

      INDEX operator[](node* n) const { return n - &nodes_[0]; }
      const node& operator[](INDEX i) const { return nodes_[i]; }

      void sort()
      {
         for(auto node : nodes_) {
            std::sort(node.begin(), node.end(), [](auto& a, auto& b) { return a.cost > b.cost; });
         }
         for(INDEX i=0; i<nodes_.size()-1; ++i) {
            assert(nodes_[i].last == nodes_[i+1].first);
         }
      }

      private:
         std::vector<node> nodes_;
         std::vector<arc> arcs_;
   };


   // graph for edges with positive cost
   struct Arc;
   struct Node {
      Arc* first; // first outgoing arc from node
   };
   struct Arc {
		Node* head;
		Arc* next;
      REAL cost; 
		Arc* prev; // do zrobienia: not needed currently
      Arc* sister; // do zrobienia: not really needed
   };
   struct Graph {
      Graph(const INDEX noNodes, const INDEX noEdges) : nodes_(noNodes,{nullptr}), arcs_(2*noEdges), noArcs_(0) {}
      INDEX size() const { return nodes_.size(); }
      // note: this is possibly not very intuitive: we have two ways to index the graph: either by Node* structors or by node numbers, returning the other structure
      // also mixture of references and pointers is ugly
      Node& operator[](const INDEX i) { assert(i<nodes_.size()); return nodes_[i]; }
      INDEX operator[](const Node* i) { assert(i - &*nodes_.begin() < nodes_.size()); return i - &*nodes_.begin(); } // do zrobienia: this is ugly pointer arithmetic
      const Node& operator[](const INDEX i) const { return nodes_[i]; }
      const INDEX operator[](const Node* i) const { assert(i - &*nodes_.begin() < nodes_.size()); return i - &*nodes_.begin(); } // do zrobienia: this is ugly pointer arithmetic
      void AddEdge(INDEX i, INDEX j, REAL cost) {
         //logger->info() << "Add edge (" << i << "," << j << ") with cost = " << cost;
         assert(noArcs_ + 1 < arcs_.size());
         Arc* a = &arcs_[noArcs_];
         ++noArcs_;
         Arc* a_rev = &arcs_[noArcs_];
         ++noArcs_;

         a->head = &(*this)[j]; a_rev->head = &(*this)[i];
         a->sister = a_rev; a_rev->sister = a;
         a->cost = cost; a_rev->cost = cost; // note: this can be made more efficient by only storing one copy of costs in a separate vector

         a->prev = nullptr;
         a->next = nodes_[i].first;
         // put away with doubly linked list
         if(nodes_[i].first != nullptr) {
            nodes_[i].first->prev = a;
         }
         nodes_[i].first = a;

         a_rev->prev = nullptr;
         a_rev->next = nodes_[j].first;
         if(nodes_[j].first != nullptr) {
            nodes_[j].first->prev = a_rev;
         }
         nodes_[j].first = a_rev;
      }
      std::vector<Node> nodes_;
      std::vector<Arc> arcs_;
      INDEX noArcs_;
   };

   // for parallelization, we must store all information for the path computations separately.
   // do zrobienia: make shortest path computation simultaneously from both ends
   struct MostViolatedPathData {
      struct Item { INDEX parent; INDEX heapPtr; INDEX flag; };
      MostViolatedPathData(const Graph& g) 
      {
         d.resize(g.size());
         for(INDEX i=0; i<d.size(); ++i) {
            d[i].flag = 0;
         }
         flagTmp = 0;
         flagPerm = 1;
      }
      void Reset() { flagTmp += 2; flagPerm += 2; pq.resize(0); }
      Item& operator[](const INDEX i) { return d[i]; }
      void LabelTemporarily(const INDEX i) { d[i].flag = flagTmp; }
      void LabelPermanently(const INDEX i) { d[i].flag = flagPerm; }
      bool LabelledTemporarily(const INDEX i) const { return d[i].flag == flagTmp; }
      bool LabelledPermanently(const INDEX i) const { return d[i].flag == flagPerm; }
      INDEX& Parent(const INDEX i) { return d[i].parent; }

      std::tuple<REAL,std::vector<INDEX>> FindPath(const INDEX startNode, const INDEX endNode, const Graph& g)
      {
         Reset();
         //logger->info() << "in finding most violated path from " << startNode << " to " << endNode;
         //PriorityQueue p;
         //p.Add(&g[startNode], std::numeric_limits<REAL>::max());
         AddPQ(startNode, std::numeric_limits<REAL>::max());
         LabelTemporarily(startNode);
         REAL curCost;
         while(!EmptyPQ()) {
            const INDEX i = RemoveMaxPQ(curCost);
            LabelPermanently(i);
            //logger->info() << "Investigating node " << g[i] << " with cost " << curCost;

            // found path
            if(i == endNode) { // trace back path to startNode
               //logger->info() << "Found shortest cycle with cost = " << curCost;
               std::vector<INDEX> path({i});
               INDEX j = i;
               while(Parent(j) != startNode) {
                  j = Parent(j);
                  path.push_back(j);
               }
               path.push_back(startNode);
               //logger->info() << "Found cycle ";
               for(INDEX i=0;i<path.size();++i) {
                  //logger->info() << path[i] << ",";
               }
               //logger->info() << " with violation " << curCost;
               return std::make_tuple(curCost,path);;
            } 

            for(Arc* a=g[i].first; a!=nullptr; a = a->next) { // do zrobienia: make iteration out of this?
               const REAL edgeCost = a->cost;
               const REAL pathCost = std::min(curCost,edgeCost); // cost of path until head
               Node* head = a->head;
               //logger->info() << " Investigating edge to " << g[head] << " with cost = " << edgeCost;
               if(LabelledPermanently(g[head])) { // permanently labelled
                  //logger->info() << "         permanently labelled";
                  continue; 
               } else if(LabelledTemporarily(g[head])) { // temporarily labelled
                  //logger->info() << "         temporarily labelled, cost = " << p.GetKey(head);
                  if(GetKeyPQ(g[head]) < pathCost) {
                     //logger->info() << "    Increasing value of node " << g[head] << " to value " << pathCost;
                     IncreaseKeyPQ(g[head], pathCost);
                     Parent(g[head]) = i;
                  }
               } else { // seen for first time
                  //logger->info() << "    Inserting node " << g[head] << " with cost " << pathCost;
                  Parent(g[head]) = i;
                  AddPQ(g[head], pathCost);
                  LabelTemporarily(g[head]);
               }
            }
         }
         return std::make_tuple(0.0,std::vector<INDEX>(0)); // no path found
      }


      private:
      std::vector<Item> d;
      INDEX flagTmp, flagPerm;

      // the priority queue used in FindPath
      struct pqItem { INDEX i; REAL key; };
      std::vector<pqItem> pq; // the priority queue
      void AddPQ(const INDEX i, const REAL key)
      {
         INDEX k = pq.size();
         pq.push_back(pqItem({i,key}));
         d[i].heapPtr = k;
         // do zrobienia: make operation heapify out of this
         while (k > 0)
         {
            INDEX k2 = (k-1)/2;
            if (pq[k2].key >= pq[k].key) break; 
            SwapPQ(k, k2);
            k = k2;
         }
      }
      REAL GetKeyPQ(const INDEX i) const
      { 
         return pq[d[i].heapPtr].key;
      }
      void IncreaseKeyPQ(const INDEX i, const REAL key) 
      {
         INDEX k = d[i].heapPtr;
         pq[k].key = key;
         while (k > 0)
         {
            INDEX k2 = (k-1)/2;
            if (pq[k2].key >= pq[k].key) break;
            SwapPQ(k, k2);
            k = k2;
         }
      }
      bool EmptyPQ() const 
      {
         return pq.size() == 0;
      }
      const INDEX RemoveMaxPQ(REAL& key) 
      {
         assert(pq.size() > 0);
         SwapPQ(0,pq.size()-1);
         INDEX k=0;

         while ( 1 )
         {
            INDEX k1 = 2*k + 1, k2 = k1 + 1;
            if (k1 >= pq.size()-1) break;
            INDEX k_max = (k2 >= pq.size()-1 || pq[k1].key >= pq[k2].key) ? k1 : k2; // do zrobienia: max-heap or min-heap
            if (pq[k].key >= pq[k_max].key) break;
            SwapPQ(k, k_max);
            k = k_max;
         }

         key = pq.back().key;
         const INDEX i = pq.back().i;
         pq.resize(pq.size()-1);
         return i;
      }
      void SwapPQ(const INDEX k1, const INDEX k2) 
      {
         pqItem& a = pq[k1];//.begin()+k1;
         pqItem& b = pq[k2];//.begin()+k2;
         d[k1].heapPtr = k2;
         d[k2].heapPtr = k1;
         std::swap(a,b);
      }
   };

   // do zrobienia: possibly one can reuse the parent structure in subsequent FindPath computations
   // However one then needs to store in union find data structures whether end node is in component. Similar can be done in most violated path search, more complicated, though
   struct BfsData {
      struct Item { INDEX parent; INDEX flag; };
      BfsData(const Graph& g) 
      {
         d.resize(g.size());
         for(INDEX i=0; i<d.size(); ++i) {
            d[i].flag = 0;
         }
         flag1 = 0;
         flag2 = 1;
      }
      void Reset() 
      {
	      flag1 += 2;
	      flag2 += 2; 
      }
      Item& operator[](const INDEX i) { return d[i]; }
      void Label1(const INDEX i) { d[i].flag = flag1; }
      void Label2(const INDEX i) { d[i].flag = flag2; }
      bool Labelled(const INDEX i) const { return Labelled1(i) || Labelled2(i); }
      bool Labelled1(const INDEX i) const { return d[i].flag == flag1; }
      bool Labelled2(const INDEX i) const { return d[i].flag == flag2; }
      INDEX& Parent(const INDEX i) { return d[i].parent; }
      INDEX Parent(const INDEX i) const { return d[i].parent; }
      std::vector<INDEX> TracePath(const INDEX i) const 
      {
	      std::vector<INDEX> path({i});
	      INDEX j=i;
	      while(Parent(j) != j) {
		      //minPathCost = std::min(minPathCost);
		      //logger->info() << j << ",";
		      j = Parent(j);
		      path.push_back(j);
	      }
         std::reverse(path.begin(),path.end()); // note: we reverse paths once too often. Possibly call this TracePathReverse. Do zrobienia.
	      return path;
      }
      
      // do bfs with thresholded costs and iteratively lower threshold until enough cycles are found
      // only consider edges that have cost equal or larger than th
      std::tuple<REAL,std::vector<INDEX>> FindPath(const INDEX startNode, const INDEX endNode, const Graph& g, const REAL th = 0, const INDEX max_length = std::numeric_limits<INDEX>::max()) 
      {
         Reset();
         std::queue<std::array<INDEX,2>> visit; // node number, distance from start or end // do zrobienia: do not allocate each time, make visit a member
         visit.push({startNode, 0});
         Label1(startNode);
         Parent(startNode) = startNode;
         visit.push({endNode, 0});
         Label2(endNode);
         Parent(endNode) = endNode;

         while(!visit.empty()) {
            const INDEX i = visit.front()[0];
            const INDEX distance = visit.front()[1];
            visit.pop();

            if(distance <= max_length) {

            if(Labelled1(i)) {
               for(Arc* a=g[i].first; a!=nullptr; a = a->next) { // do zrobienia: make iteration out of this?
                  if(a->cost >= th) {
                     Node* head = a->head;
                     const INDEX j = g[head];
                     if(!Labelled(j)) {
                        visit.push({j, distance+1});
                        Parent(j) = i;
                        Label1(j);
                     } else if(Labelled2(j)) { // shortest path found
                        // trace back path from j to endNode and from i to startNode
                        std::vector<INDEX> startPath = TracePath(i);
                        std::vector<INDEX> endPath = TracePath(j);
                        std::reverse(endPath.begin(), endPath.end());
                        startPath.insert( startPath.end(), endPath.begin(), endPath.end());
                        return std::make_tuple(th, startPath);
                     }
                  }
               }
            } else {
               assert(Labelled2(i));
               for(Arc* a=g[i].first; a!=nullptr; a = a->next) { // do zrobienia: make iteration out of this?
                  if(a->cost >= th) {
                     Node* head = a->head;
                     const INDEX j = g[head];
                     if(!Labelled(j)) {
                        visit.push({j, distance+1});
                        Parent(j) = i;
                        Label2(j);
                     } else if(Labelled1(j)) { // shortest path found
                        // trace back path from j to endNode and from i to startNode
                        std::vector<INDEX> startPath = TracePath(j);
                        std::vector<INDEX> endPath = TracePath(i);
                        std::reverse(endPath.begin(), endPath.end());
                        startPath.insert( startPath.end(), endPath.begin(), endPath.end());
                        return std::make_tuple(th, startPath);
                     }
                  }
               }
            }

            }
         }
         return std::make_tuple(th,std::vector<INDEX>(0));
      }

private:
   std::vector<Item> d;
   INDEX flag1, flag2;
};

struct BfsData2 {
      struct Item { INDEX parent; REAL cost; INDEX flag; };
      BfsData2(const Graph2& g) 
      {
         d.resize(g.size());
         for(INDEX i=0; i<d.size(); ++i) {
            d[i].flag = 0;
         }
         flag1 = 0;
         flag2 = 1;
      }
      void Reset() 
      {
         visit.clear();
	      flag1 += 2;
	      flag2 += 2; 
      }
      Item& operator[](const INDEX i) { return d[i]; }

      void Label1(const INDEX i) { d[i].flag = flag1; }
      void Label2(const INDEX i) { d[i].flag = flag2; }
      bool Labelled(const INDEX i) const { return Labelled1(i) || Labelled2(i); }
      bool Labelled1(const INDEX i) const { return d[i].flag == flag1; }
      bool Labelled2(const INDEX i) const { return d[i].flag == flag2; }

      INDEX& Parent(const INDEX i) { return d[i].parent; }
      INDEX Parent(const INDEX i) const { return d[i].parent; }
      REAL& Cost(const INDEX i) { return d[i].cost; }
      REAL Cost(const INDEX i) const { return d[i].cost; }

      std::tuple<REAL,std::vector<INDEX>> TracePath(const INDEX i1, const INDEX i2, const REAL cost_i1i2) const 
      {
         REAL c = cost_i1i2;
	      std::vector<INDEX> path({i1});
	      INDEX j=i1;
	      while(Parent(j) != j) {
            c = std::min(c, Cost(j) );
		      j = Parent(j);
		      path.push_back(j);
	      }
         std::reverse(path.begin(),path.end());
         path.push_back(i2);
         j=i2;
         while(Parent(j) != j) {
            c = std::min(c, Cost(j) );
		      j = Parent(j);
		      path.push_back(j);
         }

	      return std::move(std::make_tuple(c,std::move(path)));
      }

      // do bfs with thresholded costs and iteratively lower threshold until enough cycles are found
      // only consider edges that have cost equal or larger than th
      std::tuple<REAL,std::vector<INDEX>> FindPath(const INDEX startNode, const INDEX endNode, const Graph2& g, const REAL th = 0, const INDEX max_length = std::numeric_limits<INDEX>::max()) 
      {
         Reset();
         visit.push_back({startNode, 0});
         Label1(startNode);
         Parent(startNode) = startNode;
         Cost(startNode) = std::numeric_limits<REAL>::infinity();
         visit.push_back({endNode, 0});
         Label2(endNode);
         Parent(endNode) = endNode;
         Cost(endNode) = std::numeric_limits<REAL>::infinity();

         while(!visit.empty()) {
            const INDEX i = visit.front()[0];
            const INDEX distance = visit.front()[1];
            visit.pop_front();

            if(distance <= max_length) {

            if(Labelled1(i)) {
               for(auto* a=g[i].begin(); a->cost>=th && a!=g[i].end(); ++a) { 
                  auto* head = a->head;
                  const INDEX j = g[head];
                  if(!Labelled(j)) {
                     visit.push_back({j, distance+1});
                     Parent(j) = i;
                     Cost(j) = a->cost;
                     Label1(j);
                  } else if(Labelled2(j)) { // shortest path found
                     // trace back path from j to endNode and from i to startNode
                     return std::move(TracePath(i,j, a->cost));
                  }
               }
            } else {
               assert(Labelled2(i));
               for(auto* a=g[i].begin(); a->cost>=th && a!=g[i].end(); ++a) { 
                  auto* head = a->head;
                  const INDEX j = g[head];
                  if(!Labelled(j)) {
                     visit.push_back({j, distance+1});
                     Parent(j) = i;
                     Cost(j) = a->cost;
                     Label2(j);
                  } else if(Labelled1(j)) { // shortest path found
                     // trace back path from j to endNode and from i to startNode
                     return std::move(TracePath(i,j, a->cost));
                  }
               }
            }

            }
         }
         return std::make_tuple(th,std::vector<INDEX>(0));
      }

private:
   std::vector<Item> d;
   std::deque<std::array<INDEX,2>> visit; // node number, distance from start or end 
   INDEX flag1, flag2;
};


public:
template<typename SOLVER>
MulticutConstructor(SOLVER& pd) 
   : lp_(&pd.GetLP()),
      //pd_(pd),
   //unaryFactors_(100,hash::array2),
   tripletFactors_(100,hash::array3)
   {
      //unaryFactors_.max_load_factor(0.7);
      tripletFactors_.max_load_factor(0.7);
      constant_factor_ = new ConstantFactorContainer(0.0);
      pd.GetLP().AddFactor(constant_factor_);
   }
   ~MulticutConstructor()
   {
      static_assert(std::is_same<typename UnaryFactorContainer::FactorType, MulticutUnaryFactor>::value,"");
      static_assert(std::is_same<typename TripletFactorContainer::FactorType, MulticutTripletFactor>::value,"");
      //static_assert(std::is_same<typename MessageContainer::MessageType, MulticutUnaryTripletMessage<MessageSending::SRMP>>::value,"");
   }
   void AddToConstant(const REAL delta) { constant_factor_->GetFactor()->AddToOffset(delta); }

   virtual UnaryFactorContainer* AddUnaryFactor(const INDEX i1, const INDEX i2, const REAL cost) // declared virtual so that derived class notices when unary factor is added
   {
      //if(globalFactor_ == nullptr) {
         //globalFactor_ = new GlobalFactorContainer(MulticutGlobalFactor(), 1); // we have one element currently
      //   lp_->AddFactor(globalFactor_);
      //} else {
      //   globalFactor_->ResizeRepam(unaryFactors_.size()+1);
      //}
      assert(i1 < i2);
      assert(!HasUnaryFactor(i1,i2));
      
      auto* u = new UnaryFactorContainer(cost);
      lp_->AddFactor(u);
      auto it = unaryFactors_.insert(std::make_pair(std::array<INDEX,2>{i1,i2}, u)).first;
      unaryFactorsVector_.push_back(std::make_pair(std::array<INDEX,2>{i1,i2}, u));

      //std::cout << "current edge: (" << i1 << "," << i2 << ")";
      if(it != unaryFactors_.begin()) {
         auto prevIt = it;
         --prevIt;
         //std::cout << ", prev edge: (" << prevIt->first.operator[](0) << "," << prevIt->first.operator[](1) << ")";
         assert(prevIt->second != u);
         lp_->AddFactorRelation(prevIt->second, u);
      }
      auto nextIt = it;
      ++nextIt;
      if(nextIt != unaryFactors_.end()) {
         assert(nextIt->second != u);
         //std::cout << ", next edge: (" << nextIt->first.operator[](0) << "," << nextIt->first.operator[](1) << ")";
         lp_->AddFactorRelation(u, nextIt->second);
      }
      //std::cout << "\n";

      noNodes_ = std::max(noNodes_,std::max(i1,i2)+1);

      //LinkUnaryGlobal(u,globalFactor_,i1,i2);
      
      //logger->info() << "Add unary factor (" << i1 << "," << i2 << ") with cost = " << cost;
      return u;
   }
   UnaryFactorContainer* GetUnaryFactor(const INDEX i1, const INDEX i2) const {
      assert(HasUnaryFactor(i1,i2));
      return unaryFactors_.find(std::array<INDEX,2>{i1,i2})->second;
   }
   UnaryTripletMessageContainer* LinkUnaryTriplet(UnaryFactorContainer* u, TripletFactorContainer* t, const INDEX i) // argument i denotes which edge the unary factor connects to
   {
      assert(i < 3);
      auto* m = new UnaryTripletMessageContainer(UnaryTripletMessageType(i), u, t);
      lp_->AddMessage(m);
      return m;
   }
   //UnaryGlobalMessageContainer* LinkUnaryGlobal(UnaryFactorContainer* u, GlobalFactorContainer* g, const INDEX i1, const INDEX i2)
   //{
   //   auto* m = new UnaryGlobalMessageContainer(UnaryGlobalMessageType( g->GetFactor()->AddEdge(i1,i2) ), u, g, 0);
   //   lp_->AddMessage(m);
   //   return m;
   //}
   virtual TripletFactorContainer* AddTripletFactor(const INDEX i1, const INDEX i2, const INDEX i3) // declared virtual so that derived constructor notices when triplet factor is added
   {
      assert(i1 < i2 && i2 < i3);
      assert(!HasTripletFactor(i1,i2,i3));
      assert(HasUnaryFactor(i1,i2) && HasUnaryFactor(i1,i3) && HasUnaryFactor(i2,i3));
      auto* t = new TripletFactorContainer();
      lp_->AddFactor(t);
      tripletFactors_.insert(std::make_pair( std::array<INDEX,3>{i1,i2,i3}, t ));
      // use following ordering of unary and triplet factors: triplet comes after edge factor (i1,i2) and before (i2,i3)
      auto* before = GetUnaryFactor(i1,i2);
      lp_->AddFactorRelation(before,t);
      auto* middle = GetUnaryFactor(i1,i3);
      lp_->AddFactorRelation(middle,t);
      auto* after = GetUnaryFactor(i2,i3);
      lp_->AddFactorRelation(t,after);
      // get immediate predeccessor and successor and place new triplet in between
      //auto succ = tripletFactors_.upper_bound(std::make_tuple(i1,i2,i3));
      //if(succ != tripletFactors_.end()) {
      //   assert(t != succ->second);
      //   lp_->AddFactorRelation(t,succ->second);
      //}
      auto tripletEdges = MulticutTripletFactor::SortEdges(i1,i2,i3);
      // link with all three unary factors
      LinkUnaryTriplet(GetUnaryFactor(tripletEdges[0][0],tripletEdges[0][1]), t, 0);
      LinkUnaryTriplet(GetUnaryFactor(tripletEdges[1][0],tripletEdges[1][1]), t, 1);
      LinkUnaryTriplet(GetUnaryFactor(tripletEdges[2][0],tripletEdges[2][1]), t, 2);
      //LinkUnaryTriplet(GetUnaryFactor(i1,i2), t, 0);
      //LinkUnaryTriplet(GetUnaryFactor(i1,i3), t, 1);
      //LinkUnaryTriplet(GetUnaryFactor(i2,i3), t, 2);
      return t;
   }
   bool HasUnaryFactor(const std::tuple<INDEX,INDEX> e) const 
   {
      assert(std::get<0>(e) < std::get<1>(e));
      return unaryFactors_.find(std::array<INDEX,2>{std::get<0>(e),std::get<1>(e)}) != unaryFactors_.end();
   }
   bool HasUnaryFactor(const INDEX i1, const INDEX i2) const 
   {
      //logger->info() << "Has Unary factor with: " << i1 << ":" << i2;
      assert(i1 < i2);
      return (unaryFactors_.find(std::array<INDEX,2>{i1,i2}) != unaryFactors_.end());
   }
   bool HasTripletFactor(const INDEX i1, const INDEX i2, const INDEX i3) const 
   {
      assert(i1 < i2 && i2 < i3);
      return (tripletFactors_.find(std::array<INDEX,3>{i1,i2,i3}) != tripletFactors_.end());
   }

   TripletFactorContainer* GetTripletFactor(const INDEX i1, const INDEX i2, const INDEX i3) const 
   {
      assert(HasTripletFactor(i1,i2,i3));
      return tripletFactors_.find(std::array<INDEX,3>{i1,i2,i3})->second;
   }


   std::tuple<INDEX,INDEX> GetEdge(const INDEX i1, const INDEX i2) const
   {
      return std::make_tuple(std::min(i1,i2), std::max(i1,i2));
   }

   REAL get_edge_cost(const INDEX i1, const INDEX i2) const
   {
      assert(HasUnaryFactor(i1,i2));
      return *(unaryFactors_.find(std::array<INDEX,2>{i1,i2})->second->GetFactor());
   }

   INDEX AddCycle(std::vector<INDEX> cycle)
   {
      //return true, if cycle was not already present
      // we triangulate the cycle by looking for the smallest index -> do zrobienia!
      
      assert(cycle.size() >= 3);
      // shift cycle so that minNode becomes first element
      std::rotate(cycle.begin(), std::min_element(cycle.begin(), cycle.end()), cycle.end());
      const INDEX minNode = cycle[0];
      assert(minNode == *std::min_element(cycle.begin(), cycle.end()));
      // first we assert that the edges in the cycle are present
      for(INDEX i=0; i<cycle.size(); ++i) {
         assert(HasUnaryFactor(GetEdge(cycle[i], cycle[(i+1)%cycle.size()])));
      }
      // now we add all triplets with triangulation edge. Possibly, a better triangulation scheme would be possible
      INDEX noTripletsAdded = 0;
      for(INDEX i=2; i<cycle.size(); ++i) {
         if(!HasUnaryFactor(minNode, cycle[i])) {
            AddUnaryFactor(minNode,cycle[i],0.0);
         }
         const INDEX secondNode = std::min(cycle[i], cycle[i-1]);
         const INDEX thirdNode = std::max(cycle[i], cycle[i-1]);
         if(!HasTripletFactor(minNode, secondNode, thirdNode)) {
            AddTripletFactor(minNode, secondNode, thirdNode);
            ++noTripletsAdded;
         }
      }
      return noTripletsAdded;
   }

   // search for cycles to add such that coordinate ascent will be possible
   INDEX Tighten(const INDEX maxCuttingPlanesToAdd)
   {
      std::cout << "Search for violated triplet constraints\n";
      INDEX tripletsAdded = FindViolatedTriplets(maxCuttingPlanesToAdd);
      std::cout << "Added " << tripletsAdded << " triplet(s) out of " <<  maxCuttingPlanesToAdd << " by searching for triplets\n"; 
      //tripletsAdded = FindViolated4Cycles(maxCuttingPlanesToAdd - tripletsAdded);
      if(tripletsAdded < 0.6*maxCuttingPlanesToAdd) {
         std::cout << "Additionally search via shortest paths for violated constraints\n";
         tripletsAdded += FindViolatedCycles2(maxCuttingPlanesToAdd - tripletsAdded);
         //tripletsAdded += FindViolatedCycles(maxCuttingPlanesToAdd - tripletsAdded);
         std::cout << "Added " << tripletsAdded << " triplet(s) out of " <<  maxCuttingPlanesToAdd << " in total\n";
      }
      return tripletsAdded;
      //if(tripletsAdded == 0) {
      //   const INDEX tripletsAdded_old = FindNegativeCycles(1e-8, maxCuttingPlanesToAdd);
      //   std::cout << "Added " << tripletsAdded_old << " triplet(s) out of " <<  maxCuttingPlanesToAdd << " by old method,error\n";
      //   return tripletsAdded_old + tripletsAdded;
      //}
      //const INDEX tripletsAdded = FindNegativeCycles(minDualIncrease,maxCuttingPlanesToAdd);
      //const REAL th = FindNegativeCycleThreshold(maxCuttingPlanesToAdd);
      //if(th >= 0.0) { // otherwise no constraint can be added
      //   const INDEX tripletsAdded = FindNegativeCycles(th,maxCuttingPlanesToAdd);
      //   std::cout << "Added " << tripletsAdded << " triplet(s) out of " <<  maxCuttingPlanesToAdd << "\n";
      //   return tripletsAdded;
      //} else {
      //   std::cout << "could not find any violated cycle\n";
      //   return 0;
      //}
   }

template<
class InputIt1, class InputIt2,
      class OutputIt, class Compare, class Merge
      >
   OutputIt set_intersection_merge
(
 InputIt1 first1, InputIt1 last1,
 InputIt2 first2, InputIt2 last2,
 OutputIt d_first, Compare comp, Merge merge
 )
{
   while (first1 != last1 && first2 != last2)
   {
      if (comp(*first1, *first2))
         ++first1;
      else
      {
         if (!comp(*first2, *first1))
            *d_first++ = merge(*first1++, *first2);
         ++first2;
      }
   }
   return d_first;
}

   // search for violated triplets, e.g. triplets with one negative edge and two positive ones.
   INDEX FindViolatedTriplets(const INDEX max_triplets_to_add)
   {
      std::vector<INDEX> adjacency_list_count(noNodes_,0);
      // first determine size for adjacency_list
      for(auto& it : unaryFactorsVector_) {
         const INDEX i = std::get<0>(it.first);
         const INDEX j = std::get<1>(it.first);
         adjacency_list_count[i]++;
         adjacency_list_count[j]++; 
      }
      two_dim_variable_array<std::tuple<INDEX,REAL>> adjacency_list(adjacency_list_count);
      std::fill(adjacency_list_count.begin(), adjacency_list_count.end(), 0);
      for(auto& it : unaryFactorsVector_) {
         const INDEX i = std::get<0>(it.first);
         const INDEX j = std::get<1>(it.first);
         const REAL cost_ij = *(it.second->GetFactor());
         assert(i<j);
         adjacency_list[i][adjacency_list_count[i]] = std::make_tuple(j,cost_ij);
         adjacency_list_count[i]++;
         adjacency_list[j][adjacency_list_count[j]] = std::make_tuple(i,cost_ij);
         adjacency_list_count[j]++;
      }

      // Sort the adjacency list, for fast intersections later
      auto adj_sort = [](const auto a, const auto b) { return std::get<0>(a) < std::get<0>(b); };
//#pragma omp parallel for schedule(static)
      for(int i=0; i < adjacency_list.size(); i++) {
         std::sort(adjacency_list[i].begin(), adjacency_list[i].end(), adj_sort);
      }

      // Iterate over all of the edge intersection sets
      // do zrobienia: parallelize
      // we will intersect two adjacency list by head node, but we want to preserve the costs of either edge pointing to head node
      using intersection_type = std::tuple<INDEX,REAL,REAL>;
      auto merge = [](const auto a, const auto b) -> intersection_type { 
         assert(std::get<0>(a) == std::get<0>(b));
         return std::make_tuple(std::get<0>(a), std::get<1>(a), std::get<1>(b)); 
      };
      std::vector<intersection_type> commonNodes(noNodes_);
      std::vector<std::tuple<INDEX,INDEX,INDEX,REAL>> triplet_candidates;
      //for(auto& it : unaryFactorsVector_) 
//#pragma omp parallel for schedule(static)
      for(INDEX c=0; c<unaryFactorsVector_.size(); ++c) {
      //for(auto it=unaryFactorsVector_.begin(); it!=unaryFactorsVector_.end(); ++it) 
         const REAL cost_ij = *(unaryFactorsVector_[c].second->GetFactor());
         const INDEX i = std::get<0>(unaryFactorsVector_[c].first);
         const INDEX j = std::get<1>(unaryFactorsVector_[c].first);

         // Now find all neighbors of both i and j to see where the triangles are
         // TEMP TEMP -- fails at i=0, j=1, on i==3.
         auto intersects_iter_end = set_intersection_merge(adjacency_list[i].begin(), adjacency_list[i].end(), adjacency_list[j].begin(), adjacency_list[j].end(), commonNodes.begin(), adj_sort, merge);

         for(auto n=commonNodes.begin(); n != intersects_iter_end; ++n) {
            const INDEX k = std::get<0>(*n);

            // Since a triplet shows up three times as an edge plus
            // a node, we only consider it for the case when i<j<k 
            if(!(j<k))
               continue;
            const REAL cost_ik = std::get<1>(*n);
            const REAL cost_jk = std::get<2>(*n);

            const REAL lb = std::min(0.0, cost_ij) + std::min(0.0, cost_ik) + std::min(0.0, cost_jk);
            const REAL best_labeling = std::min({0.0, cost_ij+cost_ik, cost_ij+cost_jk, cost_ij+cost_jk, cost_ij+cost_ik+cost_jk});
            assert(lb <= best_labeling+eps);
            const REAL guaranteed_dual_increase = best_labeling - lb;
            /*
            std::array<REAL,3> c({cost_ij, cost_ik, cost_jk});
            const REAL gdi1 = std::min({ -c[0], c[1], c[2] });
            const REAL gdi2 = std::min({ c[0], -c[1], c[2] });
            const REAL gdi3 = std::min({ c[0], c[1], -c[2] });
            const REAL guaranteed_dual_increase = std::max({ gdi1, gdi2, gdi3 });
            */
            if(guaranteed_dual_increase > 0.0) {
               triplet_candidates.push_back(std::make_tuple(i,j,k,guaranteed_dual_increase));
            } 
         }
      }
      std::sort(triplet_candidates.begin(), triplet_candidates.end(), [](auto& a, auto& b) { return std::get<3>(a) > std::get<3>(b); });

      INDEX triplets_added = 0;
      for(const auto& triplet_candidate : triplet_candidates) {
         const INDEX i = std::get<0>(triplet_candidate);
         const INDEX j = std::get<1>(triplet_candidate);
         const INDEX k = std::get<2>(triplet_candidate);
         if(!HasTripletFactor(i,j,k)) {
            AddTripletFactor(i,j,k);
            triplets_added++;
            if(triplets_added > max_triplets_to_add) {
               break;
            } 
         }
      }

      return triplets_added;
   }

   INDEX FindViolated4Cycles(const INDEX max_triplets_to_add)
   {
      // other approach: explicitly construct graph, where edges are pairs of edges in the original graph. Then iterate over all parallel edges

      // first generate all simple paths of length 2. Record them in a structure indexed by end points.
      //std::unordered_multimap<std::array<INDEX,2>, INDEX, decltype(hash::array2)> two_paths(unaryFactors_.size(), hash::array2);
      //std::multimap<std::array<INDEX,2>, INDEX> two_paths;
      std::map<std::array<INDEX,2>, std::vector<INDEX>> two_paths; // possibly use two_paths_positive and two_paths_mixed to get two_paths with all positive edge and two_paths with positive and negative edge
      // build adjacency list
      std::vector<std::vector<int> > adjacency_list(noNodes_);

      // Construct adjacency list for the graph
      for(auto& it : unaryFactors_) {
         const INDEX i = std::get<0>(it.first);
         const INDEX j = std::get<1>(it.first);
         assert(i<j);
         adjacency_list[i].push_back(j);
         adjacency_list[j].push_back(i);
      }

      for(INDEX i=0; i<adjacency_list.size(); ++i) {
         for(INDEX j_idx=0; j_idx<adjacency_list[i].size(); ++j_idx) {
            const INDEX j=adjacency_list[i][j_idx];
            for(INDEX k_idx=j_idx+1; k_idx<adjacency_list[i].size(); ++k_idx) {
               const INDEX k = adjacency_list[i][k_idx];
               const REAL cost_ij = get_edge_cost(std::min(i,j), std::max(i,j));
               const REAL cost_ik = get_edge_cost(std::min(i,k), std::max(i,k));
               if(!(cost_ij < 0.0 && cost_ik < 0.0)) {
                  auto it = two_paths.find({std::min(j,k), std::max(j,k)});
                  if(it != two_paths.end()) {
                     (it->second).push_back(i);
                  } else {
                     two_paths.insert(std::make_pair(std::array<INDEX,2>{std::min(j,k), std::max(j,k)}, std::vector<INDEX>(i)));
                  } 
                  //two_paths.insert(std::make_pair(std::array<INDEX,2>{std::min(j,k), std::max(j,k)}, i)); 
               }
            }
         } 
      }

      for(const auto& end_points : two_paths) {
         const INDEX i = end_points.first[0];
         const INDEX j = end_points.first[1];
         const auto &mid_points = end_points.second;
         for(INDEX k1_idx=0; k1_idx<mid_points.size(); ++k1_idx) {
            const INDEX k1 = mid_points[k1_idx];
            for(INDEX k2_idx=k1_idx+1; k2_idx<mid_points.size(); ++k2_idx) {
               const INDEX k2 = mid_points[k2_idx];
               assert(k1 != k2);
               // check whether there is exactly one negative entry in 4-cycle
               const REAL cost_ik1 = get_edge_cost(std::min(i,k1), std::max(i,k1));
               const REAL cost_ik2 = get_edge_cost(std::min(i,k2), std::max(i,k2));
               const REAL cost_jk1 = get_edge_cost(std::min(j,k1), std::max(j,k1));
               const REAL cost_jk2 = get_edge_cost(std::min(j,k2), std::max(j,k2));


            }
         }
      }
      // Then for each two fixed end points, enumearte all pairs of distinct mid points. This gives all four cycles.
      //for (auto it=two_paths.begin(); it!=two_paths.end(); ) {
      //   const INDEX i = *(it->first)[0];
      //   const INDEX j = *(it->first)[1];
      //   auto equal_it = two_paths.equal_range(it->first);
      //   for(; equal_it.first!=equal_it.second; ++equal_it.first) {
      //      //const INDEX i = *(equal_it.first);
      //      for(auto sec_it=std::next(equal_it.first); sec_it!=equal_it.second; ++sec_it) {
      //         //const INDEX j = *(sec_it);
      //
      //          } 
      //     }

      //   /* Now, go skip to the first entry with a new key. */
      //   auto cur = it;
      //   while(it!=two_paths.end() && it->first == cur->first) ++it;
      //}
      std::cout << two_paths.size() << "\n";
      return 0;
   }

   INDEX FindViolatedCycles2(const INDEX maxTripletsToAdd)
   {
      std::vector<std::tuple<INDEX,INDEX,REAL> > negativeEdges;
      // we can speed up compution by skipping path searches for node pairs which lie in different connected components. Connectedness is stored in a union find structure
      REAL pos_th = 0.0;
      std::vector<INDEX> number_outgoing_arcs(noNodes_,0); // number of arcs outgoing arcs of each node
      INDEX number_outgoing_arcs_total = 0;
      for(auto& it : unaryFactorsVector_) {
         const REAL v = *(it.second->GetFactor());
         const INDEX i = std::get<0>(it.first);
         const INDEX j = std::get<1>(it.first);
         if(v >= 0.0) {
            pos_th = std::max(pos_th, v);
            number_outgoing_arcs[i]++;
            number_outgoing_arcs[j]++;
            number_outgoing_arcs_total += 2;
         } else {
            negativeEdges.push_back(std::make_tuple(i,j,v));
         }
      }
      if(negativeEdges.size() == 0 || negativeEdges.size() == unaryFactorsVector_.size()) { return 0; }

      Graph2 posEdgesGraph2(noNodes_, number_outgoing_arcs_total, number_outgoing_arcs); // graph consisting of positive edges
      for(auto& it : unaryFactorsVector_) {
         const REAL v = *(it.second->GetFactor());
         const INDEX i = std::get<0>(it.first);
         const INDEX j = std::get<1>(it.first);
         assert(i<j);
         if(v >= 0.0) {
            posEdgesGraph2.add_arc(i,j,v);
            posEdgesGraph2.add_arc(j,i,v);
         }
      }
      posEdgesGraph2.sort();

      // do zrobienia: possibly add reparametrization of triplet factors additionally

      std::sort(negativeEdges.begin(), negativeEdges.end(), [](const std::tuple<INDEX,INDEX,REAL>& e1, const std::tuple<INDEX,INDEX,REAL>& e2)->bool {
            return std::get<2>(e1) < std::get<2>(e2);
            });


      // now search for every negative edge for most negative path from end point to starting point. Do zrobienia: do this in parallel
      // here, longest path is sought after only the edges with positive taken into account
      // the cost of the path is the minimum of the costs of its edges. The guaranteed increase in the dual objective is v_min > -v ? -v : v_min

      //MostViolatedPathData mp(posEdgesGraph);
      //BfsData mp(posEdgesGraph);
      BfsData2 mp2(posEdgesGraph2);

      UnionFind uf(noNodes_);
      INDEX tripletsAdded = 0;
      const REAL initial_th = 0.6*std::min(-std::get<2>(negativeEdges[0]), pos_th);
      bool zero_th_iteration = true;
      //for(REAL th=initial_th; th>=eps || zero_th_iteration; th*=0.1) 
      for(REAL th=initial_th; th>=eps; th*=0.1) {
         if(th < eps) {
            if(tripletsAdded <= 0.01*maxTripletsToAdd) {
               std::cout << "additional separation with no guaranteed dual increase, i.e. th = 0\n";
               th = 0.0;
               zero_th_iteration = false;
            } else {
               break;
            }
         }
         // first update union find datastructure
         for(auto& it : unaryFactorsVector_) {
            const REAL v = *(it.second->GetFactor());
            const INDEX i = std::get<0>(it.first);
            const INDEX j = std::get<1>(it.first);
            if(v >= th) {
               uf.merge(i,j);   
            }
         }
         using CycleType = std::tuple<REAL, std::vector<INDEX>>;
         std::vector<CycleType > cycles;
         for(auto& it : negativeEdges) {
            const INDEX i = std::get<0>(it);
            const INDEX j = std::get<1>(it);
            const REAL v = std::get<2>(it);
            if(-v <= th) break;
            if(uf.connected(i,j)) {
               //auto cycle = mp.FindPath(i,j,posEdgesGraph);
               auto cycle = mp2.FindPath(i,j,posEdgesGraph2, th);
               const REAL dualIncrease = std::min(-v, std::get<0>(cycle));
               assert(std::get<1>(cycle).size() > 0);
               if(std::get<1>(cycle).size() > 0) {
                  cycles.push_back( std::make_tuple(dualIncrease, std::move(std::get<1>(cycle))) );
                  //tripletsAdded += AddCycle(std::get<1>(cycle));
                  //if(tripletsAdded > maxTripletsToAdd) {
                  //   return tripletsAdded;
                  //}
               } else {
                  throw std::runtime_error("No path found although there should be one"); 
               }
            }
         }
         // sort by guaranteed increase in decreasing order
         std::sort(cycles.begin(), cycles.end(), [](const CycleType& i, const CycleType& j) { return std::get<0>(i) > std::get<0>(j); });
         for(auto& cycle : cycles) {
            if(std::get<1>(cycle).size() > 2) {
               tripletsAdded += AddCycle(std::get<1>(cycle));
               if(tripletsAdded > maxTripletsToAdd) {
                  return tripletsAdded;
               }
            }
         }
      }

      return tripletsAdded;
   }

   // search for cycles with one negative edge and all else positive edges in descending order of guaranteed dual increase.
   INDEX FindViolatedCycles(const INDEX maxTripletsToAdd)
   {
      assert(maxTripletsToAdd > 0);
      // keep here negative edge and associated maximal decrease that we can achieve using it in cycle.
      std::multimap<REAL,std::array<INDEX,2>,std::greater<REAL>> negEdgeCandidates;

      // initialize data structures for violated cycle search
      UnionFind uf(noNodes_);
      std::vector<std::tuple<INDEX,INDEX,REAL> > posEdges;
      // do zrobienia: we contract graphs a lot. Look into Bjoern Andres graph package for efficient implementation of this operation.
      std::vector<std::vector<std::tuple<INDEX,REAL,INDEX,INDEX>>> negEdges(this->noNodes_); // forward_list would also be possible, but is slower. node entries here are the labels of connected components. Original endpoints are recorded in last two indexes.
      //std::vector<std::list<std::tuple<INDEX,REAL,INDEX,INDEX>>> negEdges(this->noNodes_); // forward_list would also be possible, but is slower. node entries here are the labels of connected components. Original endpoints are recorded in last two indexes.
      //Graph posEdgesGraph(noNodes_,unaryFactors_.size()); // graph consisting of positive edges
      std::vector<INDEX> number_outgoing_arcs(noNodes_,0); // number of arcs outgoing arcs of each node
      // sort edges 
      for(auto& it : unaryFactors_) {
         const REAL v = *(it.second->GetFactor());
         const INDEX i = std::get<0>(it.first);
         const INDEX j = std::get<1>(it.first);
         assert(i<j);
         if(v < 0) {
            negEdges[i].push_back(std::make_tuple(j,v,i,j));
            negEdges[j].push_back(std::make_tuple(i,v,i,j));
         } else if(v > 0) {
            number_outgoing_arcs[i]++;
            number_outgoing_arcs[j]++;
            posEdges.push_back(std::make_tuple(i,j,v));
            //posEdgesGraph.AddEdge(i,j,v);
         }
      }
      Graph2 posEdgesGraph2(noNodes_, 2*unaryFactors_.size(), number_outgoing_arcs);
      for(auto& it : unaryFactors_) {
         const REAL v = *(it.second->GetFactor());
         const INDEX i = std::get<0>(it.first);
         const INDEX j = std::get<1>(it.first);
         assert(i<j);
         if(v > 0) {
            posEdgesGraph2.add_arc(i,j,v);
            posEdgesGraph2.add_arc(j,i,v);
         }
      }
      posEdgesGraph2.sort();

      //BfsData mp(posEdgesGraph);
      BfsData2 mp2(posEdgesGraph2);

      std::sort(posEdges.begin(), posEdges.end(), [](const std::tuple<INDEX,INDEX,REAL>& e1, const std::tuple<INDEX,INDEX,REAL>& e2)->bool {
            return std::get<2>(e1) > std::get<2>(e2); // descending order
            });
      auto neg_edge_sort = [] (const auto& a, const auto& b)->bool { 
         if(std::get<0>(a) != std::get<0>(b)) {
            return std::get<0>(a) < std::get<0>(b); 
         } else {
            return std::get<1>(a) < std::get<1>(b); 
         }
      };
      for(INDEX i=0; i<negEdges.size(); ++i) {
         std::sort(negEdges[i].begin(), negEdges[i].end(), neg_edge_sort);
         //negEdges[i].sort(neg_edge_sort);
      }
      //std::cout << "sorted positive and negative edges for violated cycles\n";

      // now contract negative edges in descending order and check for self loops of negative edges.
      // Whenever there is one, add it to negEdgeCandidates, or replace negative edge that will not be needed anymore.
      // if maximal guaranteed increase of newly added edges is lower than some guaranteed increase of an already present edge, find cycle for the latter.
      INDEX tripletsAdded = 0;
      for(INDEX posE=0; posE<posEdges.size(); posE++) {
         const INDEX i = std::get<0>(posEdges[posE]);
         const INDEX j = std::get<1>(posEdges[posE]);
         const REAL posTh = std::get<2>(posEdges[posE]);
         //std::cout << "positive edge " << posE << ", value = " << posTh << "\n" << std::flush;
         /*
         const REAL bestCandidateTh = negEdgeCandidates.begin()->first;
         if(posTh < bestCandidateTh) {
            // search for cycles already. Note that posTh is an upper bound on the minimum dual increase of all remaining edges
            for(const auto& negEdgeIt : negEdgeCandidates) {
               const REAL th = negEdgeIt.first;
               if(th > posTh) {
                  const INDEX i = negEdgeIt.second[0];
                  const INDEX j = negEdgeIt.second[1];
                  //std::cout << "find path ..." << std::flush;
                  tripletsAdded += FindPositivePath(posEdgesGraph, mp,th,i,j, 400000000); // do not search for path longer than 200
                  if(tripletsAdded >= maxTripletsToAdd) {
                     return tripletsAdded;
                  }
                  //std::cout << "done\n" << std::flush;
               } else {
                  break;
               }
            }
            // now erase all added negative edges
            negEdgeCandidates.erase(negEdgeCandidates.begin(), negEdgeCandidates.lower_bound(posTh));
         }
         std::cout << "added " << tripletsAdded << " triplets; " << std::endl;
         */
         if(!uf.connected(i,j)) {
            //std::cout << "merge nodes..." << std::flush;
            uf.merge(i,j);
            const INDEX c = uf.find(i); // the new cc node number
            // merge the edges with bases at i and j and detect edge from i to j
            if(i != c) {
               std::vector<std::tuple<INDEX,REAL,INDEX,INDEX>> merged_vec;
               merged_vec.reserve(negEdges[c].size() + negEdges[i].size());
               //auto it = std::set_union(negEdges[c].begin(),negEdges[c].end(), negEdges[i].begin(), negEdges[i].end(), std::back_inserter(merged_vec), neg_edge_sort);
               std::merge(negEdges[c].begin(),negEdges[c].end(), negEdges[i].begin(), negEdges[i].end(), std::back_inserter(merged_vec), neg_edge_sort);
               std::swap(negEdges[c], merged_vec);
               //std::merge(negEdges[c].begin(),negEdges[c].end(), negEdges[i].begin(), begEdges[i].end(), std::back_inserter(...), neg_edge_sort);
               //negEdges[c].merge(negEdges[i],neg_edge_sort);
            } 
            if(j != c) {
               std::vector<std::tuple<INDEX,REAL,INDEX,INDEX>> merged_vec;
               merged_vec.reserve(negEdges[c].size() + negEdges[j].size());
               //auto it = std::set_union(negEdges[c].begin(),negEdges[c].end(), negEdges[j].begin(), negEdges[j].end(), std::back_inserter(merged_vec), neg_edge_sort);
               std::merge(negEdges[c].begin(),negEdges[c].end(), negEdges[j].begin(), negEdges[j].end(), std::back_inserter(merged_vec), neg_edge_sort);
               std::swap(negEdges[c], merged_vec);
               //negEdges[c].merge(negEdges[j],neg_edge_sort);
            }
            std::transform(negEdges[c].begin(), negEdges[c].end(), negEdges[c].begin(), 
                  [&uf] (auto a) {
                  std::get<0>(a) = uf.find(std::get<0>(a));
                  return a; 
                  });
            //std::cout << " ... " << negEdges[c].size() << std::flush;

            // do zrobienia: mainly unnecessary: edges have already been sorted in merge operations, only parallel edges need to be sorted here.
            //negEdges[c].sort([] (const auto& a, const auto& b)->bool { 
            //      if(std::get<0>(a) != std::get<0>(b)) {
            //      return std::get<0>(a) < std::get<0>(b); 
            //      }
            //      return std::get<1>(a) > std::get<1>(b); // this ensures that remove deleted parallel copies with smaller weight. Thus, the largest one only remains.
            //      });
            std::sort(negEdges[c].begin(), negEdges[c].end(), neg_edge_sort);
            // now go through edge list and search for self loops. Add such to negEdgeCandidates
            for(auto it=negEdges[c].begin(); it!=negEdges[c].end(); ++it) {
               const INDEX cc = std::get<0>(*it);
               std::get<0>(*it) = cc;
               assert(uf.find(cc) == cc);
               if(cc == c) { // we have found a positive edge for which a negative path exists. record the dual increase possible.
                  const REAL th = std::min(-std::get<1>(*it), std::get<2>(posEdges[posE]));
                  assert(th > 0.0);
                  const INDEX negI = std::get<2>(*it);
                  const INDEX negJ = std::get<3>(*it);
                  // insert negative edge if either there are not yet enough candidates or the candidate has a better bound than already inserted ones.
                  assert(tripletsAdded < maxTripletsToAdd);
                  if(negEdgeCandidates.size() <= (maxTripletsToAdd - tripletsAdded)) {
                     negEdgeCandidates.insert(std::make_pair(th, std::array<INDEX,2>({negI,negJ})));
                  } else { 
                     const REAL leastBestCandidateTh = negEdgeCandidates.rbegin()->first;
                     // may take too much time
                     //assert(leastBestCandidateTh == std::min_element(negEdgeCandidates.begin(), negEdgeCandidates.end(), [](auto& a, auto&b) { return a.first < b.first; })->first);
                     if(th > leastBestCandidateTh) {
                        auto last = negEdgeCandidates.end();
                        --last;
                        negEdgeCandidates.erase(last);
                        negEdgeCandidates.insert(std::make_pair(th,std::array<INDEX,2>({negI,negJ})));
                     }
                  }
               }
            }
            //std::cout << " ..." << std::flush;
            // do zrobienia: remove edges that have too small threshold already.
            //negEdges[c].remove_if([&uf,c,maxTh] (const std::tuple<INDEX,REAL>& a)->bool { return (std::get<1>(a) > -maxTh || uf.find(std::get<0>(a)) == c); });
            // remove all self loops
            auto it = std::remove_if(negEdges[c].begin(), negEdges[c].end(), [&uf, c] (const auto& a)->bool { return std::get<0>(a) == c; });
            negEdges[c].resize(std::distance(negEdges[c].begin(), it));
            //std::cout << " done\n" << std::flush;
            //negEdges[c].remove_if([&uf,c] (const auto& a)->bool { return std::get<0>(a) == c; });
            if(negEdgeCandidates.size() > maxTripletsToAdd) { 
               std::cout << "enough negative edge candidates found" << std::endl;
               break;
            }
         }
      }
      std::cout << "found negative edge candidates in violated cycles\n";
      // now go through all negative edge candidates and find path
      for(const auto& negEdge : negEdgeCandidates) {
         const REAL th = negEdge.first;
         const INDEX i = negEdge.second[0];
         const INDEX j = negEdge.second[1];
         //tripletsAdded += FindPositivePath(posEdgesGraph, mp, th, i, j);
         tripletsAdded += FindPositivePath(posEdgesGraph2, mp2, th, i, j);
         if(tripletsAdded > maxTripletsToAdd) {
            return tripletsAdded;
         }
      }
      return tripletsAdded;
   }

// find and add violated cycle with given th
template<typename GRAPH, typename BFS_STRUCT>
INDEX FindPositivePath(const GRAPH& g, BFS_STRUCT& mp, const REAL th, const INDEX i, const INDEX j, const INDEX max_length = std::numeric_limits<INDEX>::max())
{
   // this is not nice: mp must be reinitalized every time
   //MostViolatedPathData mp(g);
   //auto cycle = mp.FindPath(i,j,g);
   auto cycle = mp.FindPath(i,j,g,th, max_length);
   if(std::get<1>(cycle).size() > 1) {
      assert(std::get<1>(cycle)[0] == i);
      assert(std::get<1>(cycle).back() == j);
      //std::cout << " add path of length = " << std::get<1>(cycle).size() << "; " << std::flush;
      return AddCycle(std::get<1>(cycle));
   } else {
      return 0;
   }
}


// find threshold such that at least maxTripletsToAdd factors will be included. Assumes that graph is connected.
REAL FindNegativeCycleThreshold(const INDEX maxTripletsToAdd)
{
   // do zrobienia: reuse data structures in constraint search. make one function.
   UnionFind uf(noNodes_);
   std::vector<std::tuple<INDEX,INDEX,REAL> > posEdges;
   std::vector<std::list<std::tuple<INDEX,REAL>>> negEdges(this->noNodes_); // forward_list would also be possible, but is slower.
   // sort edges 
   for(auto& it : unaryFactors_) {
      const REAL v = it.second->operator[](0);
      const INDEX i = std::get<0>(it.first);
      const INDEX j = std::get<1>(it.first);
      if(v < 0) {
         negEdges[i].push_front(std::make_tuple(j,v));
         negEdges[j].push_front(std::make_tuple(i,v));
      } else if(v > 0) {
            posEdges.push_back(std::make_tuple(i,j,v));
         }
      }
      
      std::sort(posEdges.begin(), posEdges.end(), [](const std::tuple<INDEX,INDEX,REAL>& e1, const std::tuple<INDEX,INDEX,REAL>& e2)->bool {
            return std::get<2>(e1) > std::get<2>(e2); // descending order
            });
      auto edge_sort = [] (const std::tuple<INDEX,REAL>& a, const std::tuple<INDEX,REAL>& b)->bool { return std::get<0>(a) < std::get<0>(b); };
      for(INDEX i=0; i<negEdges.size(); ++i) {
         negEdges[i].sort(edge_sort);
      }

      // maximum violated inequality consists of cycle with one negative edge and rest positive edges. Find negative edge of smallest weight and path with minimum positive weight maximum
      // connect positive edge in descending order and see, whether there is some negative edge that has endpoint both in the same component. If yes, record the value of the cycle that can be obtained thisw ay.
      REAL maxTh = -std::numeric_limits<REAL>::infinity();
      INDEX maxThIndex = 0;
      for(INDEX posE=0; posE<posEdges.size(); posE++) {
         const INDEX i = std::get<0>(posEdges[posE]);
         const INDEX j = std::get<1>(posEdges[posE]);
         const REAL posTh = std::get<2>(posEdges[posE]);
         if(posTh < maxTh) { // all the remaining edges will not lead to maximum violated constraint
            return std::get<2>(posEdges[std::min(posE + maxTripletsToAdd,INDEX(posEdges.size()-1))]);
         }
         if(!uf.connected(i,j)) {
            uf.merge(i,j);
            const INDEX c = uf.find(i); // the new cc node number
            // merge the edges with bases at i and j and detect edge from i to j
            if(i != c) {
               negEdges[c].merge(negEdges[i],edge_sort);
            } 
            if(j != c) {
               negEdges[c].merge(negEdges[j],edge_sort);
            }
            std::transform(negEdges[c].begin(), negEdges[c].end(), negEdges[c].begin(), 
                  [&uf] (std::tuple<INDEX,REAL> a)->std::tuple<INDEX,REAL> {
                  std::get<0>(a) = uf.find(std::get<0>(a));
                  return a; 
                  });
            // do zrobienia: mainly unnecessary: edges have already been sorted in merge operations, only parallel edges need to be sorted here.
            negEdges[c].sort([] (const std::tuple<INDEX,REAL>& a, const std::tuple<INDEX,REAL>& b)->bool { 
                  if(std::get<0>(a) != std::get<0>(b)) {
                  return std::get<0>(a) < std::get<0>(b); 
                  }
                  return std::get<1>(a) > std::get<1>(b); // this ensures that remove removes copies with smaller weight. Thus, the largest one only remains.
                  });
            // now go through edge list and remove duplicates
            for(auto it=negEdges[c].begin(); it!=negEdges[c].end(); ++it) {
               const INDEX cc = std::get<0>(*it);
               std::get<0>(*it) = cc;
               if(uf.find(cc) == c) { // we have found a positive edge for which a negative path exists. remove this element and record that the dual increase possible.
                  const REAL th = std::min(-std::get<1>(*it), std::get<2>(posEdges[posE]));
                  if(th > maxTh) {
                     maxTh = std::max(th, maxTh);
                     maxThIndex = posE;
                  }
                  continue;
               }
            }
            negEdges[c].remove_if([&uf,c,maxTh] (const std::tuple<INDEX,REAL>& a)->bool { return (std::get<1>(a) > -maxTh || uf.find(std::get<0>(a)) == c); });
            negEdges[c].unique([] (const std::tuple<INDEX,REAL>& a, const std::tuple<INDEX,REAL>& b)->bool { return std::get<0>(a) == std::get<0>(b); });
         }
      }
      if(maxTh == -std::numeric_limits<REAL>::infinity()) { // this means that no 
         return maxTh;
      } else {
         // return threshold that will guarantee maximum element. Small factor allows more constraints to be included.
         assert(false);
         return 0.001*maxTh;
      }
   }

   // returns number of triplets added
   INDEX FindNegativeCycles(const REAL minDualIncrease, const INDEX maxTripletsToAdd)
   {
      //assert(minDualIncrease > 0.0);
      std::vector<std::tuple<INDEX,INDEX,REAL> > negativeEdges;
      // we can speed up compution by skipping path searches for node pairs which lie in different connected components. Connectedness is stored in a union find structure
      UnionFind uf(noNodes_);
      std::vector<INDEX> number_outgoing_arcs(noNodes_,0); // number of arcs outgoing arcs of each node
      INDEX number_outgoing_arcs_total = 0;
      //Graph posEdgesGraph(noNodes_,unaryFactors_.size()); // graph consisting of positive edges
      for(auto& it : unaryFactors_) {
         const REAL v = *(it.second->GetFactor());
         const INDEX i = std::get<0>(it.first);
         const INDEX j = std::get<1>(it.first);
         if(v >= minDualIncrease) {
            //posEdgesGraph.AddEdge(i,j,v);
            uf.merge(i,j);
            number_outgoing_arcs[i]++;
            number_outgoing_arcs[j]++;
            number_outgoing_arcs_total += 2;
         }
      }
      Graph2 posEdgesGraph2(noNodes_, number_outgoing_arcs_total, number_outgoing_arcs); // graph consisting of positive edges
      for(auto& it : unaryFactors_) {
         const REAL v = *(it.second->GetFactor());
         const INDEX i = std::get<0>(it.first);
         const INDEX j = std::get<1>(it.first);
         assert(i<j);
         if(v >= minDualIncrease) {
            posEdgesGraph2.add_arc(i,j,v);
            posEdgesGraph2.add_arc(j,i,v);
         }
      }
      posEdgesGraph2.sort();


      for(auto& it : unaryFactors_) {
         const REAL v = *(it.second->GetFactor());
         const INDEX i = std::get<0>(it.first);
         const INDEX j = std::get<1>(it.first);
         if(v < -minDualIncrease && uf.connected(i,j)) {
            negativeEdges.push_back(std::make_tuple(i,j,v));
         }
      }
      // do zrobienia: possibly add reparametrization of triplet factors additionally

      std::sort(negativeEdges.begin(), negativeEdges.end(), [](const std::tuple<INDEX,INDEX,REAL>& e1, const std::tuple<INDEX,INDEX,REAL>& e2)->bool {
            return std::get<2>(e1) < std::get<2>(e2);
            });

      // now search for every negative edge for most negative path from end point to starting point. Do zrobienia: do this in parallel
      // here, longest path is sought after only the edges with positive taken into account
      // the cost of the path is the minimum of the costs of its edges. The guaranteed increase in the dual objective is v_min > -v ? -v : v_min

      //MostViolatedPathData mp(posEdgesGraph);
      //BfsData mp(posEdgesGraph);
      BfsData2 mp2(posEdgesGraph2);

      // better: collect all cycles, then sort them according to violation, and only then add them

      INDEX tripletsAdded = 0;
      using CycleType = std::tuple<REAL, std::vector<INDEX>>;
      std::vector<CycleType > cycles;
      for(auto& it : negativeEdges) {
         const INDEX i = std::get<0>(it);
         const INDEX j = std::get<1>(it);
         const REAL v = std::get<2>(it);
         //if(-v < minDualIncrease) break;
         //auto cycle = mp.FindPath(i,j,posEdgesGraph);
         auto cycle = mp2.FindPath(i,j,posEdgesGraph2);
         const REAL dualIncrease = std::min(-v, std::get<0>(cycle));
         assert(std::get<1>(cycle).size() > 0);
         if(std::get<1>(cycle).size() > 0) {
            cycles.push_back( std::make_tuple(dualIncrease, std::move(std::get<1>(cycle))) );
            //tripletsAdded += AddCycle(std::get<1>(cycle));
            if(tripletsAdded > maxTripletsToAdd) {
               return tripletsAdded;
            }
         } else {
            throw std::runtime_error("No path found although there should be one"); 
         }
      }
      // sort by guaranteed increase in decreasing order
      std::sort(cycles.begin(), cycles.end(), [](const CycleType& i, const CycleType& j) { return std::get<0>(i) > std::get<0>(j); });
      for(auto& cycle : cycles) {
         const REAL cycleDualIncrease = std::get<0>(cycle);
         if(std::get<1>(cycle).size() > 2) {
            tripletsAdded += AddCycle(std::get<1>(cycle));
            if(tripletsAdded > maxTripletsToAdd) {
               return tripletsAdded;
            }
         }
      }

      return tripletsAdded;
   }

   bool CheckPrimalConsistency(PrimalSolutionStorage::Element primal) const
   {
      //std::cout << "checking primal feasibility for multicut\n";
      UnionFind uf(noNodes_);
      for(const auto& e : unaryFactorsVector_) {
         UnaryFactorContainer* f = e.second; 
         assert(primal[f->GetPrimalOffset()] != unknownState);
         if(primal[f->GetPrimalOffset()] == false) {
            // connect components 
            const INDEX i = std::get<0>(e.first);
            const INDEX j = std::get<1>(e.first);
            uf.merge(i,j);
         }
      }
      for(const auto& e : unaryFactorsVector_) {
         UnaryFactorContainer* f = e.second; 
         if(primal[f->GetPrimalOffset()] == true) {
            const INDEX i = std::get<0>(e.first);
            const INDEX j = std::get<1>(e.first);
            // there must not be a path from i1 to i2 consisting of edges with primal value false only
            if(uf.connected(i,j)) {
               std::cout << "solution infeasible: (" << i << "," << j << ") = true, yet there exists a path with false values only\n";
               return false;
            }
         }
      }

      std::cout << "solution feasible\n";
      return true;
   }

/*
   void ComputePrimal() const
   {
      // do zrobienia: templatize for correct graph type
      // do zrobienia: put original graph into multicut constructor and let it be constant, i.e. not reallocate it every time for primal computation. Problem: When adding additional edges, we may not add them to the lifted multicut solver, as extra edges must not participate in cut inequalities

      // use GAEC and Kernighan&Lin algorithm of andres graph package to compute primal solution
      const INDEX noNodes = noNodes_;
      andres::graph::Graph<> graph(noNodes);
      std::vector<REAL> edgeValues;
      edgeValues.reserve(unaryFactors_.size());

      for(const auto& e : unaryFactors_) {
         graph.insertEdge(e.first[0], e.first[1]);
         edgeValues.push_back(*(e.second->GetFactor()));
      }

      std::vector<char> labeling(unaryFactors_.size(),0);
      andres::graph::multicut::greedyAdditiveEdgeContraction(graph,edgeValues,labeling);
      andres::graph::multicut::kernighanLin(graph,edgeValues,labeling,labeling);

      // now write back primal solution and evaluate cost. 
      INDEX i=0; // the index in labeling
      for(const auto& e : unaryFactors_) {
         const auto* f = e.second;
         assert(false);
         //f->SetAndPropagatePrimal(primal, labeling.begin()+i);
         ++i;
      }
   }
   */



protected:
   //GlobalFactorContainer* globalFactor_;
   // do zrobienia: replace this by unordered_map, provide hash function.
   // possibly dont do this, but use sorting to provide ordering for LP
   std::map<std::array<INDEX,2>, UnaryFactorContainer*> unaryFactors_; // actually unary factors in multicut are defined on edges. assume first index < second one
   std::vector<std::pair<std::array<INDEX,2>, UnaryFactorContainer*>> unaryFactorsVector_; // we store a second copy of unary factors for faster iterating
   //std::unordered_map<std::array<INDEX,2>, UnaryFactorContainer*, decltype(hash::array2)> unaryFactors_; // actually unary factors in multicut are defined on edges. assume first index < second one
   // sort triplet factors as follows: Let indices be i=(i1,i2,i3) and j=(j1,j2,j3). Then i<j iff i1+i2+i3 < j1+j2+j3 or for ties sort lexicographically
   struct tripletComp {
      bool operator()(const std::tuple<INDEX,INDEX,INDEX> i, const std::tuple<INDEX,INDEX,INDEX> j) const
      {
         const INDEX iSum = std::get<0>(i) + std::get<1>(i) + std::get<2>(i);
         const INDEX jSum = std::get<0>(j) + std::get<1>(j) + std::get<2>(j);
         if(iSum < jSum) return true;
         else if(iSum > jSum) return false;
         else return i<j; // lexicographic comparison
      }
   };
   std::unordered_map<std::array<INDEX,3>, TripletFactorContainer*, decltype(hash::array3)> tripletFactors_; // triplet factors are defined on cycles of length three
   INDEX noNodes_ = 0;
   ConstantFactorContainer* constant_factor_;

   LP* lp_;
};



template<class MULTICUT_CONSTRUCTOR, INDEX TRIPLET_PLUS_SPOKE_FACTOR_NO, INDEX TRIPLET_PLUS_SPOKE_MESSAGE_NO, INDEX TRIPLET_PLUS_SPOKE_COVER_MESSAGE_NO>
class MulticutOddWheelConstructor : public MULTICUT_CONSTRUCTOR {
   using FMC = typename MULTICUT_CONSTRUCTOR::FMC;
   using BaseConstructor = MULTICUT_CONSTRUCTOR;

   using TripletPlusSpokeFactorContainer = meta::at_c<typename FMC::FactorList, TRIPLET_PLUS_SPOKE_FACTOR_NO>;
   using TripletPlusSpokeMessageContainer = typename meta::at_c<typename FMC::MessageList, TRIPLET_PLUS_SPOKE_MESSAGE_NO>::MessageContainerType;
   using TripletPlusSpokeCoverMessageContainer = typename meta::at_c<typename FMC::MessageList, TRIPLET_PLUS_SPOKE_COVER_MESSAGE_NO>::MessageContainerType;
public:
   template<typename SOLVER>
   MulticutOddWheelConstructor(SOLVER& pd) : BaseConstructor(pd), tripletPlusSpokeFactors_(100,hash::array4) { 
      tripletPlusSpokeFactors_.max_load_factor(0.7); 
   }

   // add triplet indices additionally to tripletIndices_
   virtual typename BaseConstructor::UnaryFactorContainer* AddUnaryFactor(const INDEX i1, const INDEX i2, const REAL cost)
   {
      //if(BaseConstructor::noNodes_ > unaryIndices_.size()) {
      //   unaryIndices_.resize(BaseConstructor::noNodes_);
      //}
      if(BaseConstructor::noNodes_ > tripletByIndices_.size()) {
         tripletByIndices_.resize(BaseConstructor::noNodes_);
      }
      typename BaseConstructor::UnaryFactorContainer* u = BaseConstructor::AddUnaryFactor(i1,i2,cost);   
      //unaryIndices_[i1].push_back(std::make_tuple(i2,u));
      //unaryIndices_[i2].push_back(std::make_tuple(i1,u));
      return u;
   }

   // add triplet indices additionally to tripletIndices_
   virtual typename BaseConstructor::TripletFactorContainer* AddTripletFactor(const INDEX i1, const INDEX i2, const INDEX i3)
   {
      assert(i1 < i2 && i2 < i3);
      assert(i3 < BaseConstructor::noNodes_);
      if(BaseConstructor::noNodes_ > tripletByIndices_.size()) {
         tripletByIndices_.resize(BaseConstructor::noNodes_);
      }
      auto* t = BaseConstructor::AddTripletFactor(i1,i2,i3);
      tripletByIndices_[i1].push_back(std::make_tuple(i2,i3,t));
      tripletByIndices_[i2].push_back(std::make_tuple(i1,i3,t));
      tripletByIndices_[i3].push_back(std::make_tuple(i1,i2,t));
      return t;
   }

   void CycleNormalForm(std::vector<INDEX>& cycle) const
   {
      assert(cycle.size() >= 3);
      assert(cycle.size()%2 == 1);
      // first search for smallest entry and make it first
      std::rotate(cycle.begin(), std::min_element(cycle.begin(), cycle.end()), cycle.end());
      // now two choices left: we can traverse cycle in forward or backward direction. Choose direction such that second entry is smaller than in reverse directoin.
      if(cycle[1] > cycle.back()) {
         std::reverse(cycle.begin()+1, cycle.end());
      }
   }

   TripletPlusSpokeFactorContainer* AddTripletPlusSpokeFactor(const INDEX n1, const INDEX n2, const INDEX centerNode, const INDEX spokeNode)
   {
      assert(!HasTripletPlusSpokeFactor(n1,n2,centerNode,spokeNode));
      auto* tps = new TripletPlusSpokeFactorContainer();
      BaseConstructor::lp_->AddFactor(tps);
      assert(n1<n2);
      tripletPlusSpokeFactors_.insert(std::make_pair(std::array<INDEX,4>({n1,n2,centerNode,spokeNode}),tps));
      std::array<INDEX,3> tripletIndices{n1,n2,centerNode};
      std::sort(tripletIndices.begin(), tripletIndices.end());
      auto* t = BaseConstructor::GetTripletFactor(tripletIndices[0], tripletIndices[1], tripletIndices[2]);
      auto* m = new TripletPlusSpokeCoverMessageContainer(MulticutTripletPlusSpokeCoverMessage(n1,n2,centerNode,spokeNode), t, tps);
      BaseConstructor::lp_->AddMessage(m);
      //BaseConstructor::lp_->AddFactorRelation(t,tps);
      return tps;
   }
   bool HasTripletPlusSpokeFactor(const INDEX n1, const INDEX n2, const INDEX centerNode, const INDEX spokeNode) const
   {
      assert(n1<n2);
      return tripletPlusSpokeFactors_.find(std::array<INDEX,4>({n1,n2,centerNode,spokeNode})) != tripletPlusSpokeFactors_.end();
   }
   TripletPlusSpokeFactorContainer* GetTripletPlusSpokeFactor(const INDEX n1, const INDEX n2, const INDEX centerNode, const INDEX spokeNode) const
   {
      assert(n1<n2);
      assert(HasTripletPlusSpokeFactor(n1,n2,centerNode,spokeNode));
      return tripletPlusSpokeFactors_.find(std::array<INDEX,4>({n1,n2,centerNode,spokeNode}))->second;
   }
   // tripletNode is either n1 or n2, depending on which triplet we want to take
   TripletPlusSpokeMessageContainer* LinkTripletPlusSpokeFactor(const INDEX n1, const INDEX n2, const INDEX centerNode, const INDEX spokeNode, const INDEX tripletNode)
   {
      assert(n1<n2);
      assert(tripletNode == n1 || tripletNode == n2);
      std::array<INDEX,3> tripletIndices{tripletNode,centerNode,spokeNode};
      std::sort(tripletIndices.begin(), tripletIndices.end());
      auto* t = BaseConstructor::GetTripletFactor(tripletIndices[0], tripletIndices[1], tripletIndices[2]);
      auto* tps = GetTripletPlusSpokeFactor(n1,n2,centerNode,spokeNode);
      
      // the edges in TripletPlusSpokeFactor are ordered as (n1,n2),(n1,centerNode),(n2,centerNode),(centerNode,spokeNode). We seek a permutation mapping the triplet edges onto the first three nodes of the TripletPlusSpoke edges
      auto te = MulticutTripletFactor::SortEdges(centerNode,spokeNode,tripletNode);
      // do zrobienia: make this a static function in TripletPlusSpokeFactor
      auto tspE = MulticutTripletPlusSpokeFactor::SortEdges(n1,n2,centerNode,spokeNode);
      // now find last spokeEdge in triplet:
      const INDEX spokeEdgeTriplet = std::find(te.begin(), te.end(), tspE[3]) - te.begin();
      assert(spokeEdgeTriplet < 3);
      INDEX sharedTripletEdgeTriplet = 3;
      INDEX sharedTripletEdgeTripletPlusSpoke;
      for(INDEX i=0; i<3; ++i) {
         const INDEX tripletIdx = std::find(te.begin(), te.end(), tspE[i]) - te.begin();
         if(tripletIdx < 3) {
            sharedTripletEdgeTriplet = tripletIdx;
            sharedTripletEdgeTripletPlusSpoke = i;
            break;
         }
      }
      assert(sharedTripletEdgeTriplet < 3);

      auto* m = new TripletPlusSpokeMessageContainer(MulticutTripletPlusSpokeMessage(n1,n2,centerNode,spokeNode,tripletIndices[0],tripletIndices[1],tripletIndices[2]), t, tps);
      //auto* m = new TripletPlusSpokeMessageContainer(MulticutTripletPlusSpokeMessage(sharedTripletEdgeTriplet,spokeEdgeTriplet,sharedTripletEdgeTripletPlusSpoke),t,tps,MulticutTripletPlusSpokeMessage::size());
      BaseConstructor::lp_->AddMessage(m);
      if(tripletNode == n1) {
         //BaseConstructor::lp_->AddFactorRelation(t,tps);
      } else {
         assert(tripletNode == n2);
         //BaseConstructor::lp_->AddFactorRelation(tps,t);
      }
      return m;
   }

   INDEX EnforceOddWheel(const INDEX centerNode, std::vector<INDEX> cycle)
   {
      assert(cycle.size()%2 == 1);
      CycleNormalForm(cycle);
      //logger->info() << "Enforce odd wheel with center node " << centerNode << " and cycle nodes ";
      //for(auto i : cycle) logger->info() << i << ",";
      for(auto i : cycle) { assert(i != centerNode); }
      // add emppty triplets emanating from center node
      for(INDEX i=0; i<cycle.size(); ++i) {
         std::array<INDEX,3> triplet{centerNode, cycle[i], cycle[(i+1)%cycle.size()]};
         std::sort(triplet.begin(), triplet.end()); // do zrobienia: use faster sorting
         if(!BaseConstructor::HasTripletFactor(triplet[0],triplet[1],triplet[2])) {
            AddTripletFactor(triplet[0],triplet[1],triplet[2]);
         }
      }
   
      INDEX tripletPlusSpokesAdded = 0;
      // we enforce the odd wheel constraint by quadrangulating the odd wheel. This results in triplets and triplets with an additional spoke.
      // The spoke is chosen to be the edge (centerNode,cycle[0]);
      std::array<INDEX,2> spoke{centerNode,cycle[0]};
      std::sort(spoke.begin(), spoke.end());
      // for quadrangulating, we need to artificially create triangles with two edges on the cycle plus the spoke
      for(INDEX i=1; i<cycle.size()-1; ++i) {
         auto e = BaseConstructor::GetEdge(cycle[0],cycle[i]);
         if(!BaseConstructor::HasUnaryFactor(e)) {
            AddUnaryFactor(std::get<0>(e), std::get<1>(e),0.0);
         }
         std::array<INDEX,3> tripletIndices{centerNode,cycle[i],cycle[0]};
         std::sort(tripletIndices.begin(), tripletIndices.end());
         if(!BaseConstructor::HasTripletFactor(tripletIndices[0],tripletIndices[1],tripletIndices[2])) {
            AddTripletFactor(tripletIndices[0],tripletIndices[1],tripletIndices[2]);
         }
      }
      // now add all needed TripletPlusSpoke factors and join them to the two triplets that are relevant
      for(INDEX i=1; i<cycle.size()-1; ++i) {
         // the edge {centerNode,cycle[0]} is the spoke
         INDEX n1 = cycle[i];
         INDEX n2 = cycle[i+1];
         if(n1>n2) {
            std::swap(n1,n2); 
         }
         if(!HasTripletPlusSpokeFactor(n1,n2,centerNode, cycle[0])) {
            AddTripletPlusSpokeFactor(n1,n2,centerNode,cycle[0]);
            LinkTripletPlusSpokeFactor(n1,n2,centerNode,cycle[0],n1);
            LinkTripletPlusSpokeFactor(n1,n2,centerNode,cycle[0],n2);
            ++tripletPlusSpokesAdded;
         }
      }

      return tripletPlusSpokesAdded;
   }

   // node i is center node, (j,k) is cycle edge
   REAL ComputeTriangleTh(const INDEX i, const INDEX j, const INDEX k)
   {
      std::array<INDEX,3> triplet{i,j,k};
      std::sort(triplet.begin(),triplet.end()); // do zrobienia: use faster sorting
      assert(this->HasUnaryFactor(triplet[0],triplet[1]) && this->HasUnaryFactor(triplet[0],triplet[2]) && this->HasUnaryFactor(triplet[1],triplet[2]));
      std::array<REAL,4> cost;
      std::fill(cost.begin(),cost.end(),0.0);
      if(this->HasTripletFactor(triplet[0],triplet[1],triplet[2])) {
         auto* t = this->GetTripletFactor(triplet[0],triplet[1],triplet[2])->GetFactor();
         assert(t->size() == 4);
         cost[0] += (*t)[0];
         cost[1] += (*t)[1];
         cost[2] += (*t)[2];
         cost[3] += (*t)[3];
      } 
      // get cost directly from edge factors
      const REAL c01 = *(this->GetUnaryFactor(triplet[0], triplet[1])->GetFactor());
      const REAL c02 = *(this->GetUnaryFactor(triplet[0], triplet[2])->GetFactor());
      const REAL c12 = *(this->GetUnaryFactor(triplet[1], triplet[2])->GetFactor());
      cost[0] += c02 + c12;
      cost[1] += c01 + c12;
      cost[2] += c01 + c02;
      cost[3] += c01 + c02 + c12;

      assert(j<k); // if not, below computation is not valid
      // compute difference between cost such that exactly one edge incident to center node is 1 againt cost when when zero or two incident to it are 1
      if(std::min(j,k) == triplet[0] && std::max(j,k) == triplet[1]) { // jk is first edge
         return std::min(0.0,std::min(cost[3],cost[0])) - std::min(cost[1],cost[2]);
      } else if(std::min(j,k) == triplet[0] && std::max(j,k) == triplet[2]) { // jk is second edge
         return std::min(0.0,std::min(cost[3],cost[1])) - std::min(cost[0],cost[2]);
      } else { // jk is third edge
         assert(std::min(j,k) == triplet[1] && std::max(j,k) == triplet[2]);
         return std::min(0.0,std::min(cost[3],cost[2])) - std::min(cost[0],cost[1]);
      }
      assert(false);
   }

   template<typename ADJ_LIST>
   void ComputeTriangles(const INDEX i, const ADJ_LIST& adjacencyList, const REAL minTh, 
         std::unordered_map<INDEX,INDEX>& origToCompressedNode, 
         std::vector<INDEX>& compressedToOrigNode, 
         std::vector<std::tuple<INDEX,INDEX,REAL>>& compressedEdges)
   {
      //std::unordered_map<INDEX,INDEX> origToCompressedNode {tripletByIndices_[i].size()}; // compresses node indices
      origToCompressedNode.max_load_factor(0.7);
      //std::vector<INDEX> compressedToOrigNode; // compressed nodes to original
      //std::vector<std::tuple<INDEX,INDEX,REAL>> compressedEdges;
      origToCompressedNode.clear();
      compressedToOrigNode.clear();
      compressedEdges.clear();

      // find all triangles ijk
      std::vector<INDEX> commonNodes(this->noNodes_); // for detecting triangles
      for(INDEX j : adjacencyList[i]) {
         auto commonNodesEnd = std::set_intersection(adjacencyList[i].begin(), adjacencyList[i].end(), adjacencyList[j].begin(), adjacencyList[j].end(), commonNodes.begin());
         for(auto it=commonNodes.begin(); it!=commonNodesEnd; ++it) {
            const INDEX k = *it; 
            if(j<k) { // edge is encountered twice
               const REAL th = ComputeTriangleTh(i,j,k);
               if(th >= minTh) {
                  // add j and k to compressed nodes
                  if(origToCompressedNode.find(j) == origToCompressedNode.end()) {
                     origToCompressedNode.insert(std::make_pair(j, origToCompressedNode.size()));
                     compressedToOrigNode.push_back(j);
                  }
                  if(origToCompressedNode.find(k) == origToCompressedNode.end()) {
                     origToCompressedNode.insert(std::make_pair(k, origToCompressedNode.size()));
                     compressedToOrigNode.push_back(k);
                  }
                  const INDEX jc = origToCompressedNode[j];
                  const INDEX kc = origToCompressedNode[k];
                  assert(jc != kc);
                  compressedEdges.push_back(std::make_tuple(jc,kc,th));
               }
            }
         }
      }
   }

   template<typename ADJ_LIST>
   REAL ComputeThreshold(const INDEX i, const ADJ_LIST& adjacencyList)
   {
      std::unordered_map<INDEX,INDEX> origToCompressedNode {tripletByIndices_[i].size()}; // compresses node indices
      std::vector<INDEX> compressedToOrigNode; // compressed nodes to original
      std::vector<std::tuple<INDEX,INDEX,REAL>> compressedEdges;
      ComputeTriangles(i, adjacencyList, 0.0, origToCompressedNode, compressedToOrigNode, compressedEdges);

      std::sort(compressedEdges.begin(), compressedEdges.end(), [](auto a, auto b) { return std::get<2>(a) > std::get<2>(b); });

      const INDEX noCompressedNodes = origToCompressedNode.size();
      const INDEX noBipartiteCompressedNodes = 2*noCompressedNodes;
      UnionFind uf(noBipartiteCompressedNodes);
      // construct bipartite graph based on triangles
      for(auto& e : compressedEdges) {
         const INDEX jc = std::get<0>(e);
         const INDEX kc = std::get<1>(e);
         uf.merge(jc,noCompressedNodes + kc);
         uf.merge(noCompressedNodes + jc,kc);
         if(uf.connected(jc, noCompressedNodes + jc)) {
            assert(uf.connected(kc, noCompressedNodes + kc));
            return std::get<2>(e);
         }
      }
      return -std::numeric_limits<REAL>::infinity(); // no constraint found
   }

   // returns nodes of odd wheel without center node
   template<typename ADJ_LIST>
   std::vector<INDEX> ComputeViolatedOddWheel(const INDEX i, const REAL minTh, const ADJ_LIST& adjacencyList)
   {
      std::unordered_map<INDEX,INDEX> origToCompressedNode {tripletByIndices_[i].size()}; // compresses node indices
      std::vector<INDEX> compressedToOrigNode; // compressed nodes to original
      std::vector<std::tuple<INDEX,INDEX,REAL>> compressedEdges;
      ComputeTriangles(i, adjacencyList, minTh, origToCompressedNode, compressedToOrigNode, compressedEdges);

      const INDEX noCompressedNodes = origToCompressedNode.size();
      const INDEX noBipartiteCompressedNodes = 2*noCompressedNodes;
      std::vector<INDEX> no_outgoing_arcs(2*noCompressedNodes, 0);
      for(const auto e : compressedEdges) {
         const INDEX i = std::get<0>(e);
         const INDEX j = std::get<1>(e);
         no_outgoing_arcs[i]++;
         no_outgoing_arcs[noCompressedNodes + i]++;
         no_outgoing_arcs[j]++;
         no_outgoing_arcs[noCompressedNodes + j]++;
      }
      typename BaseConstructor::Graph2 g(noBipartiteCompressedNodes,4*compressedEdges.size(), no_outgoing_arcs);
      typename BaseConstructor::BfsData2 mp(g);
      UnionFind uf(noBipartiteCompressedNodes);
      // construct bipartite graph based on triangles
      for(auto& e : compressedEdges) {
         assert(std::get<2>(e) >= minTh);
         const INDEX jc = std::get<0>(e);
         const INDEX  kc = std::get<1>(e);
         g.add_arc(jc, noCompressedNodes + kc,0.0);
         g.add_arc(noCompressedNodes + kc, jc,0.0);
         g.add_arc(noCompressedNodes + jc, kc,0.0);
         g.add_arc(kc, noCompressedNodes + jc,0.0);
         uf.merge(jc,noCompressedNodes + kc);
         uf.merge(noCompressedNodes + jc,kc);
      }
      g.sort();
      // now check whether path exists between any given edges on graph
      for(INDEX j=0; j<noCompressedNodes; ++j) { // not nice: this has to be original number of nodes and bipartiteNumberOfNodes
         // find path from node j to node noNodes+j in g
         if(uf.connected(j,noCompressedNodes+j)) {
            auto path = mp.FindPath(j,noCompressedNodes+j,g);
            auto& pathNormalized = std::get<1>(path);
            pathNormalized.resize(pathNormalized.size()-1); // first and last node coincide
            for(INDEX k=0; k<pathNormalized.size(); ++k) { // note: last node is copy of first one
               //assert(compressedToOrigNode.find(pathNormalized[k]%noCompressedNodes) != compressedToOrigNode.end());
               pathNormalized[k] = compressedToOrigNode[pathNormalized[k]%noCompressedNodes];
            }

            if(HasUniqueValues(pathNormalized)) { // possibly already add the subpath that is unique and do not search for it later. Indicate this with a std::vector<bool>
               //assert(HasUniqueValues(pathNormalized)); // if not, a shorter subpath has been found. This subpath will be detected or has been deteced and has been added
               CycleNormalForm(pathNormalized);
               //CycleNormalForm called unnecesarily in EnforceOddWheel
               return std::move(pathNormalized);
            } 
         }
      }
      return std::vector<INDEX>(0);
   }

   INDEX FindOddWheels(const INDEX maxCuttingPlanesToAdd)
   {
      // first prepare datastructures for threshold finding and violated constraint search

      // search for all triangles present in the graph, also in places where no triplet factor has been added
      // Construct adjacency list


      std::vector<INDEX> adjacency_list_count(this->noNodes_,0);
      // first determine size for adjacency_list
      for(auto& it : this->unaryFactorsVector_) {
         const INDEX i = std::get<0>(it.first);
         const INDEX j = std::get<1>(it.first);
         adjacency_list_count[i]++;
         adjacency_list_count[j]++; 
      }
      two_dim_variable_array<INDEX> adjacency_list(adjacency_list_count);
      std::fill(adjacency_list_count.begin(), adjacency_list_count.end(), 0);
      for(auto& it : this->unaryFactorsVector_) {
         const INDEX i = std::get<0>(it.first);
         const INDEX j = std::get<1>(it.first);
         assert(i<j);
         adjacency_list[i][adjacency_list_count[i]] = j;
         adjacency_list_count[i]++;
         adjacency_list[j][adjacency_list_count[j]] = i;
         adjacency_list_count[j]++;
      }

      // Sort the adjacency list, for fast intersections later
      // do zrobienia: parallelize
      for(int i=0; i < adjacency_list.size(); i++) {
         std::sort(adjacency_list[i].begin(), adjacency_list[i].end());
      }

      // find maximum threshold where still some cycle can be added for each node
      std::vector<std::tuple<INDEX,REAL>> threshold(this->noNodes_);
      // do zrobienia: parallelize
      for(INDEX i=0; i<this->noNodes_; ++i) {
         // compute cost difference of all triplets with i as one of its nodes. 
         // sort in descending order.
         // populate union find datastructure successively with them, checking whether opposite nodes in bipartite graph are connected. If so, threshold is found
         threshold[i] = std::make_tuple(i, ComputeThreshold(i,adjacency_list));
      }

      //std::cout << "maximal odd wheel threshold = " << std::get<1>(*std::max_element(threshold.begin(), threshold.end(), [](auto a, auto b) { return std::get<1>(a) > std::get<1>(b); })) << "\n\n";
      // find the maxCuttingPlanesToAdd nodes with largest guaranteed dual increase.
      assert( threshold.size() > 0 && maxCuttingPlanesToAdd > 0);
      const INDEX n = std::min(INDEX(threshold.size()), maxCuttingPlanesToAdd) - 1;
      auto sort_func = [](auto a, auto b) { return std::get<1>(a) > std::get<1>(b); };
      std::nth_element(threshold.begin(), threshold.begin() + n, threshold.end(), sort_func);
      std::sort(threshold.begin(), threshold.begin()+n, sort_func);

      // go over all nodes with large threshold, compute optimum odd wheel and try to add it to lp
      INDEX factorsAdded = 0;
      // do zrobienia: parallelize
      for(auto it=threshold.begin(); it!=threshold.begin()+n; ++it) {
         const INDEX i = std::get<0>(*it);
         const REAL th = std::get<1>(*it);
         if(th >= 0.0) {
            auto oddWheel = ComputeViolatedOddWheel(i,th, adjacency_list);
            assert(oddWheel.size() > 0);
            factorsAdded += EnforceOddWheel(i,oddWheel);
         }
      }
      return factorsAdded;

      /*

      assert(false);
      REAL minDualIncrease = 0.0; // for now, remove later

      INDEX oddWheelsAdded = 0;
      std::vector<INDEX> commonNodes(BaseConstructor::noNodes_); // for detecting triangles
      std::vector<std::tuple<REAL,std::vector<INDEX>>> oddWheels; // record violated odd wheels here for later sorting and inserting
      for(INDEX i=0; i<this->noNodes_; ++i) {
         // compressed nodes for later usage in bfs
         // do zrobienia: origToCompressedNode and reverse is reallocated all the time. Allocate once and reset.
         std::unordered_map<INDEX,INDEX> origToCompressedNode {tripletByIndices_[i].size()}; // compresses node indices
         origToCompressedNode.max_load_factor(0.7);
         std::vector<INDEX> compressedToOrigNode; // compressed nodes to original
         std::vector<std::array<INDEX,2>> compressedEdges;

         // find all triangles ijk
         for(INDEX j : adjacencyList[i]) {
            auto commonNodesEnd = std::set_intersection(adjacencyList[i].begin(), adjacencyList[i].end(), adjacencyList[j].begin(), adjacencyList[j].end(), commonNodes.begin());
            for(auto it=commonNodes.begin(); it!=commonNodesEnd; ++it) {
               const INDEX k = *it; 
               if(j<k) { // edge is encountered twice
                  // compute cost of all possible triangle assignments on ijk
                  bool addTriplet = false;
                  std::array<INDEX,3> triplet{i,j,k};
                  std::sort(triplet.begin(),triplet.end()); // do zrobienia: use faster sorting
                  std::array<REAL,4> tripletCost {0.0,0.0,0.0,0.0};
                  if(this->HasTripletFactor(triplet[0],triplet[1],triplet[2])) {
                     auto* t = this->GetTripletFactor(triplet[0],triplet[1],triplet[2])->GetFactor();
                     INDEX l1, l2;
                     INDEX l3, l4; // the other labelings
                     assert(t->size() == 4);
                     if(i < j && i < k) { // the cycle edge is the last one
                        l1 = 0; l2 = 1;
                        l3 = 2; l4 = 3;
                     } else if (i > j && i > k) { // the cycle edge is the first one
                        l1 = 1; l2 = 2;
                        l3 = 0; l4 = 3;
                     } else { // j < i < k, the cycle edge is the second one
                        l1 = 0; l2 = 2;
                        l3 = 1; l4 = 3;
                     }
                     if( std::min((*t)[l1], (*t)[l2]) <= std::min(0.0,std::min((*t)[l3],(*t)[l4])) - minDualIncrease) { // do zrobienia: check again
                        addTriplet = true;
                     }
                  } else { // get cost directly from edge factors
                     const REAL ij = *(this->GetUnaryFactor(std::min(i,j), std::max(i,j))->GetFactor());
                     const REAL ik = *(this->GetUnaryFactor(std::min(i,k), std::max(i,k))->GetFactor());
                     const REAL jk = *(this->GetUnaryFactor(j,k)->GetFactor());
                     if(std::min(ij+jk, ik+jk) <= std::min(std::min(ij+ik, ij+ik+jk),0.0) - minDualIncrease) { // do zrobienia: check again
                        addTriplet = true;
                     }
                  }
                  if(addTriplet == true) {
                     // add j and k to compressed nodes
                     if(origToCompressedNode.find(j) == origToCompressedNode.end()) {
                        origToCompressedNode.insert(std::make_pair(j, origToCompressedNode.size()));
                        compressedToOrigNode.push_back(j);
                     }
                     if(origToCompressedNode.find(k) == origToCompressedNode.end()) {
                        origToCompressedNode.insert(std::make_pair(k, origToCompressedNode.size()));
                        compressedToOrigNode.push_back(k);
                     }
                     const INDEX jc = origToCompressedNode[j];
                     const INDEX kc = origToCompressedNode[k];
                     assert(jc != kc);
                     compressedEdges.push_back(std::array<INDEX,2>{jc,kc});

                  }
               }
            }
         }
         const INDEX noCompressedNodes = origToCompressedNode.size();
         const INDEX noBipartiteCompressedNodes = 2*noCompressedNodes;
         typename BaseConstructor::Graph g(noBipartiteCompressedNodes,2*compressedEdges.size());
         typename BaseConstructor::BfsData mp(g);
         UnionFind uf(noBipartiteCompressedNodes);
         // construct bipartite graph based on triangles
         for(auto& e : compressedEdges) {
            const INDEX jc = e[0];
            const INDEX  kc = e[1];
            g.AddEdge(jc,noCompressedNodes + kc,0.0);
            g.AddEdge(noCompressedNodes + jc,kc,0.0);
            uf.merge(jc,noCompressedNodes + kc);
            uf.merge(noCompressedNodes + jc,kc);
         }
         // now check whether path exists between any given edges on graph
         //logger->info() << "built auxiliary graph for node " << i << " with " << origToCompressedNode.size() << " nodes";
         for(INDEX j=0; j<noCompressedNodes; ++j) { // not nice: this has to be original number of nodes and bipartiteNumberOfNodes
            // find path from node j to node noNodes+j in g
            if(uf.connected(j,noCompressedNodes+j)) {
               //logger->info() << "find path from " << j << " to " << j+noNodes << " and add corresponding wheel";
               auto path = mp.FindPath(j,noCompressedNodes+j,g);
               auto pathNormalized = std::get<1>(path);
               //logger->info() << "found compressed path ";
               //for(INDEX k=0; k<pathNormalized.size(); ++k) {
               //   logger->info() << pathNormalized[k] << ", ";
               //}
               pathNormalized.resize(pathNormalized.size()-1); // first and last node coincide
               for(INDEX k=0; k<pathNormalized.size(); ++k) { // note: last node is copy of first one
                  //assert(compressedToOrigNode.find(pathNormalized[k]%noCompressedNodes) != compressedToOrigNode.end());
                  pathNormalized[k] = compressedToOrigNode[pathNormalized[k]%noCompressedNodes];
               }
               //for(auto& k : pathNormalized) {
               //   assert(compressedToOrigNode.find(k) != compressedToOrigNode.end());
               //   k = compressedToOrigNode[k%noNodes];
               //}
               if(HasUniqueValues(pathNormalized)) { // possibly already add the subpath that is unique and do not search for it later. Indicate this with a std::vector<bool>
                  //assert(HasUniqueValues(pathNormalized)); // if not, a shorter subpath has been found. This subpath will be detected or has been deteced and has been added
                  CycleNormalForm(pathNormalized);
                  //CycleNormalForm called unnecesarily in EnforceOddWheel
                  oddWheelsAdded += EnforceOddWheel(i,pathNormalized);
               } else {
                  //spdlog::get("logger")->info() << "kwaskwas: add subcycles";
                  //assert(false); //
               }
            }
         }
      }

      return oddWheelsAdded;

      //logger->info() << "find odd wheel, " << BaseConstructor::tripletFactors_.size();
      // do zrobienia: use reparametrization of edge potentials or of triplet potentials or of both simultaneously?
      // currently we assume we have triplet edges, which may not hold true. Better use original edges
      for(INDEX i=0; i<BaseConstructor::noNodes_; ++i) {
         // get all triplet factors attached to node i and build bipartite subgraph with doubled edges as described in Nowozin's thesis
         std::unordered_map<INDEX,INDEX> origToCompressedNode {tripletByIndices_[i].size()}; // compresses node indices
         origToCompressedNode.max_load_factor(0.7);
         std::vector<INDEX> compressedToOrigNode; // compressed nodes to original
         for(INDEX j=0; j<tripletByIndices_[i].size(); ++j) {
            const INDEX i1 = std::get<0>(tripletByIndices_[i][j]);
            const INDEX i2 = std::get<1>(tripletByIndices_[i][j]);
            assert(i1 != i);
            assert(i2 != i);
            assert(i1 != i2);
            if(origToCompressedNode.find(i1) == origToCompressedNode.end()) {
               origToCompressedNode.insert(std::make_pair(i1, origToCompressedNode.size()));
               compressedToOrigNode.push_back(i1);
            }
            if(origToCompressedNode.find(i2) == origToCompressedNode.end()) {
               origToCompressedNode.insert(std::make_pair(i2, origToCompressedNode.size()));
               compressedToOrigNode.push_back(i2);
            }
         }
         const INDEX noCompressedNodes = origToCompressedNode.size();
         const INDEX noBipartiteCompressedNodes = 2*noCompressedNodes;
         typename BaseConstructor::Graph g(noBipartiteCompressedNodes,2*tripletByIndices_[i].size());
         typename BaseConstructor::BfsData mp(g);
         UnionFind uf(noBipartiteCompressedNodes);
         for(INDEX j=0; j<tripletByIndices_[i].size(); ++j) {
            auto*  t = std::get<2>(tripletByIndices_[i][j])->GetFactor();
            const INDEX i1 = std::get<0>(tripletByIndices_[i][j]);
            const INDEX i2 = std::get<1>(tripletByIndices_[i][j]);
            // check if 110 or 101 are among the minimal labelings, where the first edge is the outer cycle edge
            // the corresponding labeling numbers are
            INDEX l1, l2;
            INDEX l3, l4; // the other labelings
            assert(t->size() == 4);
            if(i < i1 && i < i2) { // the cycle edge is the last one
               l1 = 0; l2 = 1;
               l3 = 2; l4 = 3;
            } else if (i > i1 && i > i2) { // the cycle edge is the first one
               l1 = 1; l2 = 2;
               l3 = 0; l4 = 3;
            } else { // i1 < i < i2, the cycle edge is the second one
               l1 = 0; l2 = 2;
               l3 = 1; l4 = 3;
            }
            if( std::min((*t)[l1], (*t)[l2]) <= std::min(0.0,std::min((*t)[l3],(*t)[l4])) - minDualIncrease) {
               const INDEX i1c = origToCompressedNode[i1];
               const INDEX i2c = origToCompressedNode[i2];
               assert(i1c != i2c);
               g.AddEdge(i1c,noCompressedNodes + i2c,0.0);
               g.AddEdge(noCompressedNodes + i1c,i2c,0.0);
               uf.merge(i1c,noCompressedNodes + i2c);
               uf.merge(noCompressedNodes + i1c,i2c);
            }
         }
         // now check whether path exists between any given edges on graph
         //logger->info() << "built auxiliary graph for node " << i << " with " << origToCompressedNode.size() << " nodes";
         for(INDEX j=0; j<noCompressedNodes; ++j) { // not nice: this has to be original number of nodes and bipartiteNumberOfNodes
            // find path from node j to node noNodes+j in g
            if(uf.connected(j,noCompressedNodes+j)) {
               //logger->info() << "find path from " << j << " to " << j+noNodes << " and add corresponding wheel";
               auto path = mp.FindPath(j,noCompressedNodes+j,g);
               auto pathNormalized = std::get<1>(path);
               //logger->info() << "found compressed path ";
               //for(INDEX k=0; k<pathNormalized.size(); ++k) {
               //   logger->info() << pathNormalized[k] << ", ";
               //}
               pathNormalized.resize(pathNormalized.size()-1); // first and last node coincide
               for(INDEX k=0; k<pathNormalized.size(); ++k) { // note: last node is copy of first one
                  //assert(compressedToOrigNode.find(pathNormalized[k]%noCompressedNodes) != compressedToOrigNode.end());
                  pathNormalized[k] = compressedToOrigNode[pathNormalized[k]%noCompressedNodes];
               }
               //for(auto& k : pathNormalized) {
               //   assert(compressedToOrigNode.find(k) != compressedToOrigNode.end());
               //   k = compressedToOrigNode[k%noNodes];
               //}
               if(HasUniqueValues(pathNormalized)) { // possibly already add the subpath that is unique and do not search for it later. Indicate this with a std::vector<bool>
                  //assert(HasUniqueValues(pathNormalized)); // if not, a shorter subpath has been found. This subpath will be detected or has been deteced and has been added
                  CycleNormalForm(pathNormalized);
                  //CycleNormalForm called unnecesarily in EnforceOddWheel
                  oddWheelsAdded += EnforceOddWheel(i,pathNormalized);
               } else {
                  //spdlog::get("logger")->info() << "kwaskwas: add subcycles";
                  //assert(false); //
               }
            }
         }

      }
      
      return oddWheelsAdded;
      */
   }

   INDEX Tighten(const INDEX maxCuttingPlanesToAdd)
   {
      const INDEX tripletsAdded = BaseConstructor::Tighten(maxCuttingPlanesToAdd);
      if(tripletsAdded > 0) {
         return tripletsAdded;
      } else {
         // possibly require the odd wheels to have larger impact relative to violated cycles. This ensures that odd wheels are only added late in the optimziation, when no good violated cycles are present any more
         const INDEX oddWheelsAdded = FindOddWheels(maxCuttingPlanesToAdd);
         std::cout << "Added " << oddWheelsAdded << " factors for odd wheel constraints\n";
         return oddWheelsAdded;
      }
      assert(false);
   }

   

private:
   std::vector<std::vector<std::tuple<INDEX,INDEX,typename BaseConstructor::TripletFactorContainer*>>> tripletByIndices_; // of triplet factor with indices (i1,i2,i3) exists, then (i1,i2,i3) will be in the vector of index i1, i2 and i3
   // the format for TripletPlusSpoke is (node1,node2, centerNode, spokeNode) and we assume n1<n2
   // hash for std::array<INDEX,4>
   std::unordered_map<std::array<INDEX,4>,TripletPlusSpokeFactorContainer*,decltype(hash::array4)> tripletPlusSpokeFactors_;
};

template<class MULTICUT_CONSTRUCTOR, INDEX LIFTED_MULTICUT_CUT_FACTOR_NO, INDEX CUT_EDGE_LIFTED_MULTICUT_FACTOR_NO, INDEX LIFTED_EDGE_LIFTED_MULTICUT_FACTOR_NO>
class LiftedMulticutConstructor : public MULTICUT_CONSTRUCTOR {
public:
   using FMC = typename MULTICUT_CONSTRUCTOR::FMC;
   using LiftedMulticutCutFactorContainer = meta::at_c<typename FMC::FactorList, LIFTED_MULTICUT_CUT_FACTOR_NO>;
   using CutEdgeLiftedMulticutFactorMessageContainer = typename meta::at_c<typename FMC::MessageList, CUT_EDGE_LIFTED_MULTICUT_FACTOR_NO>::MessageContainerType;
   using LiftedEdgeLiftedMulticutFactorMessageContainer = typename meta::at_c<typename FMC::MessageList, LIFTED_EDGE_LIFTED_MULTICUT_FACTOR_NO>::MessageContainerType;

   // do zrobienia: use this everywhere instead of std::array<INDEX,2>
   struct Edge : public std::array<INDEX,2> {
      Edge(const INDEX i, const INDEX j) : std::array<INDEX,2>({std::min(i,j), std::max(i,j)}) {}
   };
   using CutId = std::vector<Edge>;

   template<typename SOLVER>
   LiftedMulticutConstructor(SOLVER& pd) : MULTICUT_CONSTRUCTOR(pd) {}

   virtual typename MULTICUT_CONSTRUCTOR::UnaryFactorContainer* AddUnaryFactor(const INDEX i1, const INDEX i2, const REAL cost)
   {
      assert(i1<i2);
      auto* f = MULTICUT_CONSTRUCTOR::AddUnaryFactor(i1,i2,cost);
      if(!addingTighteningEdges) {
         baseEdges_.push_back({i1,i2,f});
      }
      return f;
   }
   typename MULTICUT_CONSTRUCTOR::UnaryFactorContainer* AddLiftedUnaryFactor(const INDEX i1, const INDEX i2, const REAL cost)
   {
      auto* f = MULTICUT_CONSTRUCTOR::AddUnaryFactor(i1,i2,cost);
      liftedEdges_.push_back({i1,i2,f}); 
      return f;
   }

   bool HasCutFactor(const CutId& cut) 
   {
      assert(std::is_sorted(cut.begin(), cut.end()));
      return liftedMulticutFactors_.find(cut) != liftedMulticutFactors_.end();
   }

   bool HasLiftedEdgeInCutFactor(const CutId& cut, const INDEX i1, const INDEX i2)
   {
      assert(i1<i2);
      assert(HasCutFactor(cut));
      const auto& edgeList = liftedMulticutFactors_[cut].second;
      // speed this up by holding edge list sorted
      return std::find(edgeList.begin(),edgeList.end(),Edge({i1,i2})) != edgeList.end();
   }

   // do zrobienia: provide AddCutFactor(const CutId& cut, const INDEX i1, const INDEX i2) as well
   LiftedMulticutCutFactorContainer* AddCutFactor(const CutId& cut)
   {
      assert(!HasCutFactor(cut));
      //std::cout << "Add cut with edges ";
      //for(auto i : cut) { std::cout << "(" << std::get<0>(i) << "," << std::get<1>(i) << ");"; } std::cout << "\n";
      auto* f = new LiftedMulticutCutFactorContainer(cut.size());
      MULTICUT_CONSTRUCTOR::lp_->AddFactor(f);
      // connect the cut edges
      for(INDEX e=0; e<cut.size(); ++e) {
         auto* unaryFactor = MULTICUT_CONSTRUCTOR::GetUnaryFactor(cut[e][0],cut[e][1]);
         auto* m = new CutEdgeLiftedMulticutFactorMessageContainer(CutEdgeLiftedMulticutFactorMessage(e),unaryFactor,f);
         MULTICUT_CONSTRUCTOR::lp_->AddMessage(m);
      }
      liftedMulticutFactors_.insert(std::make_pair(cut,std::make_pair(f,std::vector<Edge>())));
      return f;
   }

   LiftedMulticutCutFactorContainer* GetCutFactor(const CutId& cut)
   {
      assert(HasCutFactor(cut));
      return liftedMulticutFactors_[cut].first;
   }

   void AddLiftedEdge(const CutId& cut, const INDEX i1, const INDEX i2)
   {
      assert(!HasLiftedEdgeInCutFactor(cut,i1,i2));
      assert(HasCutFactor(cut));
      assert(i1<i2);
      //std::cout << "Add lifted edge (" << i1 << "," << i2 << ") to cut\n";
      auto& c = liftedMulticutFactors_[cut];
      auto* f = c.first;
      auto& edgeList = c.second;
      auto* unaryFactor = MULTICUT_CONSTRUCTOR::GetUnaryFactor(i1,i2);
      f->GetFactor()->IncreaseLifted();
      auto* m = new LiftedEdgeLiftedMulticutFactorMessageContainer(LiftedEdgeLiftedMulticutFactorMessage(edgeList.size() + cut.size()), unaryFactor, f);
      MULTICUT_CONSTRUCTOR::lp_->AddMessage(m);
      c.second.push_back(Edge({i1,i2}));
   }



   INDEX Tighten(const INDEX maxCuttingPlanesToAdd)
   {
      const bool prevMode = addingTighteningEdges;
      addingTighteningEdges = true;
      assert(maxCuttingPlanesToAdd > 5); //otherwise the below arrangement makes little sense.
      const INDEX noBaseConstraints = MULTICUT_CONSTRUCTOR::Tighten(0.8*maxCuttingPlanesToAdd);
      //return noBaseConstraints;
      INDEX noLiftingConstraints = 0;
      std::cout << "number of cut constraints = " << liftedMulticutFactors_.size() << "\n";
      if(noBaseConstraints < maxCuttingPlanesToAdd) {
         REAL th = FindViolatedCutsThreshold(maxCuttingPlanesToAdd - noBaseConstraints);
         if(th >= 0.0) {
            noLiftingConstraints = FindViolatedCuts(th, maxCuttingPlanesToAdd - noBaseConstraints);
            std::cout << "added " << noLiftingConstraints << " lifted cut factors.\n";
         }
      }
      addingTighteningEdges = prevMode;
      return noBaseConstraints + noLiftingConstraints;
   }

   REAL FindViolatedCutsThreshold(const INDEX maxTripletsToAdd)
   {
      // make one function to reuse allocated datastructures.
      UnionFind uf(this->noNodes_);
      std::vector<std::tuple<INDEX,INDEX,REAL>> edges;
      edges.reserve(baseEdges_.size());
      for(const auto& e : baseEdges_) {
         const INDEX i = e.i;
         const INDEX j = e.j;
         const REAL weight = e.weight();
         if(weight > 0) {
            uf.merge(i,j);
         } else {
            edges.push_back(std::make_tuple(i,j,weight));
         }
      }
      std::sort(edges.begin(),edges.end(), [] (auto& a, auto& b)->bool { return std::get<2>(a) > std::get<2>(b); });
      std::vector<std::list<std::tuple<INDEX,REAL>>> liftedEdges(this->noNodes_);
      for(auto& e : liftedEdges_) {
         const INDEX i = e.i;
         const INDEX j = e.j;
         const REAL weight = e.weight();
         if(weight > 0.0) {
            liftedEdges[i].push_front(std::make_tuple(j,weight));
            liftedEdges[j].push_front(std::make_tuple(i,weight));
         }
      }
      auto edge_sort = [] (const std::tuple<INDEX,REAL>& a, const std::tuple<INDEX,REAL>& b)->bool { return std::get<0>(a) < std::get<0>(b); };
      for(INDEX i=0; i<liftedEdges.size(); ++i) {
         liftedEdges[i].sort(edge_sort);
      }
      REAL prevWeight = 0.0;
      REAL maxTh = -std::numeric_limits<REAL>::infinity();
      for(INDEX e=0; e<edges.size(); ++e) {
         const INDEX i = std::get<0>(edges[e]);
         const INDEX j = std::get<1>(edges[e]);
         const REAL weight = std::get<2>(edges[e]);
         if(uf.find(i) != uf.find(j)) {
            uf.merge(i,j);
            const INDEX c = uf.find(i);
            // (i) merge edges emanating from same connected component
            if(i != c) {
               liftedEdges[c].merge(liftedEdges[i],edge_sort);
            } 
            if(j != c) {
               liftedEdges[c].merge(liftedEdges[j],edge_sort);
            }
            std::transform(liftedEdges[c].begin(), liftedEdges[c].end(), liftedEdges[c].begin(), 
                  [&uf,&maxTh,weight] (std::tuple<INDEX,REAL> a)->std::tuple<INDEX,REAL> {
                  const INDEX cc = uf.find(std::get<0>(a));
                  if(cc != std::get<0>(a)) {
                  const REAL th = std::min(-weight, std::get<1>(a));
                  maxTh = std::max(th,maxTh); // do zrobienia: correct?
                  }
                  std::get<0>(a) = cc;
                  return a; 
                  });
            // do zrobienia: unnecessary sort here: has been sorted more or less after merge, only parallel edges need to be sorted
            liftedEdges[c].sort([] (const std::tuple<INDEX,REAL>& a, const std::tuple<INDEX,REAL>& b)->bool { 
                  if(std::get<0>(a) != std::get<0>(b)) {
                  return std::get<0>(a) < std::get<0>(b); 
                  }
                  return std::get<1>(a) > std::get<1>(b); // this ensures that remove removes copies with smaller weight. Thus, the largest one only remains.
                  });
            // (ii) remove lifted edges that have been collapsed. do zrobienia: record weight of those
            liftedEdges[c].remove_if([&uf,c,maxTh] (const auto& a)->bool { 
                  if(std::get<1>(a) < maxTh) return true; // we need not take this edge into account anymore, as it would lead to smaller minimum dual improvement.
                  // first check whether lifted edge belonged to different connected components before the last merge operation. If yes, update threshold
                  const INDEX cc = uf.find(std::get<0>(a));
                  return uf.find(std::get<0>(a)) == c; 
                  });
            // (iii) take maximal representative of edges that are parallel now
            liftedEdges[c].unique([] (const auto& a, const auto& b)->bool { return std::get<0>(a) == std::get<0>(b); }); 
         }

      }

      //if(maxTh != -std::numeric_limits<REAL>::infinity()) {
      //   std::cout << "\nnon-trivial cut factor with weight = " << maxTh << "\n\n";
      //}
      return 0.1*maxTh;
   }

   INDEX FindViolatedCuts(const INDEX minDualIncrease, const INDEX noConstraints)
   {
      // update weight of base edges
      //for(auto& e : baseEdges_) {
      //   e.w = MULTICUT_CONSTRUCTOR::GetUnaryFactor(e.i,e.j)->operator[](0);
      //}
      //std::sort(baseEdges_.begin(), baseEdges_.end(), [](const MulticutEdge& e1, const MulticutEdge& e2) { return e1.weight() < e2.weight(); });
      UnionFind uf(MULTICUT_CONSTRUCTOR::noNodes_);
      for(const auto& e : baseEdges_) {
         if(e.weight() >= -minDualIncrease) {
            uf.merge(e.i,e.j);
         }
      }

      // build reduced graph with connected components as nodes and edge weights as number of edges with weight < -minDualIncrease

      // union find indices are not contiguous. Make them so, to use them as identifiers for connected components
      std::map<INDEX,INDEX> ufIndexToContiguous;
      for(INDEX i=0; i<MULTICUT_CONSTRUCTOR::noNodes_; ++i) {
         const INDEX ufIndex = uf.find(i);
         if(ufIndexToContiguous.find(ufIndex) == ufIndexToContiguous.end()) {
            ufIndexToContiguous.insert(std::make_pair(ufIndex,ufIndexToContiguous.size()));
         }
      }
      const INDEX ccNodes = ufIndexToContiguous.size();

      std::map<INDEX,INDEX> origToCompressedNode; // union find index to compressed node indices // do zrobienia: use hash map
      for(INDEX i=0; i<MULTICUT_CONSTRUCTOR::noNodes_; ++i) {
         const INDEX ufIndex = uf.find(i);
         const INDEX collapsedIndex = ufIndexToContiguous[ufIndex];
         origToCompressedNode.insert(std::make_pair(i,collapsedIndex));
      }
      
      INDEX ccEdges = 0;
      std::map<Edge,std::vector<Edge>> ccToBaseEdges;
      for(const auto& e : baseEdges_) {
         if(!uf.connected(e.i,e.j)) {
            const INDEX i = ufIndexToContiguous[uf.find(e.i)];
            const INDEX j = ufIndexToContiguous[uf.find(e.j)];
            //if(ccEdgeCap.find({i,j}) == cc.EdgeCap.end()) {
            //   ccEdgeCap.insert(std::make_pair(std::array<INDEX,2>(i,j),1));
            //}
            //ccEdgeCap[std::array<INDEX,2>(i,j)]++;
            if(ccToBaseEdges.find(Edge(i,j)) == ccToBaseEdges.end()) {
               ++ccEdges;
               ccToBaseEdges.insert(std::make_pair(Edge(i,j),std::vector<Edge>()));
            }
            ccToBaseEdges[Edge(i,j)].push_back(Edge(e.i,e.j));
         }
      }
      BKMaxFlow::Graph<int,int,int> maxFlow(ccNodes, ccEdges); 
      maxFlow.add_node(ccNodes);
      for(auto& e : ccToBaseEdges) {
         const INDEX i = e.first.operator[](0);
         const INDEX j = e.first.operator[](1);
         const INDEX cap = e.second.size();
         maxFlow.add_edge(i,j,cap,cap);
      }
      
      // note: this can possibly be made faster by caching the weight
      std::sort(liftedEdges_.begin(), liftedEdges_.end(), [](const MulticutEdge& e1, const MulticutEdge& e2) { return e1.weight() > e2.weight(); });
      INDEX factorsAdded = 0;

      // note that currently possibly multiple max flow computations are performed, when lifted edges come from the same connected components. This is superfluous and searches could be remembered and reused.
      const int capacityMax = baseEdges_.size()+1;
      for(const auto& liftedEdge : liftedEdges_) {
         //spdlog::get("logger")->info() << "considering lifted edge " << liftedEdge.i << "," << liftedEdge.j;
         if(factorsAdded >= noConstraints) { 
            //spdlog::get("logger")->info() << "maximal number of constraints to add reached";
            break; 
         }
         if(liftedEdge.weight() > minDualIncrease) {
            const INDEX i = origToCompressedNode[liftedEdge.i];
            const INDEX j = origToCompressedNode[liftedEdge.j];
            if(!uf.connected(i,j)) {
               // find minimum cut in unweighted graph containing only base edges with weight < eps
               maxFlow.add_tweights(i,capacityMax,0);
               maxFlow.add_tweights(j,0,capacityMax);
               const INDEX noCutEdges = maxFlow.maxflow();
               assert(noCutEdges > 0 && noCutEdges < baseEdges_.size()); // otherwise there is no path from i to j or all paths were collapsed
               std::vector<Edge> minCut;
               minCut.reserve(noCutEdges);
               // now do a dfs from i on those vertices which are in the same segment as i. The edges from those to vertices in segment of j form a minimum cut.
               std::stack<INDEX> q;
               std::vector<bool> visited(ccNodes,false);
               q.push(i);
               //spdlog::get("logger")->info() << "finding max flow between node " << i << " and " << j;
               while(!q.empty()) {
                  const INDEX v = q.top();
                  q.pop();
                  if(visited[v]) {
                     continue;
                  }
                  visited[v] = true;
                  auto* a = maxFlow.get_first_arc(v);
                  while(a != nullptr) {
                     int v_test, w; // do zrobienia: use proper type
                     maxFlow.get_arc_ends(a, v_test, w);
                     assert(v != w);
                     assert(v_test == v);
                     if(maxFlow.what_segment(INDEX(v)) != maxFlow.what_segment(INDEX(w))) {
                        // expand all edges that were collapsed into (v,w) in the original graph
                           //spdlog::get("logger")->info() << "edge in mincut: " << v << "," << w;
                        for(const auto& e : ccToBaseEdges[Edge(INDEX(v),INDEX(w))]) {
                           //spdlog::get("logger")->info() << " expanded edge : " << e[0] << "," << e[1];
                           minCut.push_back(Edge(e[0],e[1]));
                        }
                     } else if(!visited[INDEX(w)]) {
                        q.push(INDEX(w));
                     }
                     a = a->next;
                  }
               }
               assert(minCut.size() == noCutEdges);
               std::sort(minCut.begin(),minCut.end()); // unique form of cut

               // check if minimum cut is already present
               if(!HasCutFactor(minCut)) {
                  auto* f = AddCutFactor(minCut);
                  AddLiftedEdge(minCut,liftedEdge.i,liftedEdge.j);
                  ++factorsAdded;
               } else {
                  auto* MinCutFactor = GetCutFactor(minCut);
                  if(!HasLiftedEdgeInCutFactor(minCut,liftedEdge.i,liftedEdge.j)) {
                     AddLiftedEdge(minCut,liftedEdge.i,liftedEdge.j);
                     ++factorsAdded;
                  }
               }

               // restore original terminal weights
               maxFlow.add_tweights(i,-capacityMax,0);
               maxFlow.add_tweights(j,0,-capacityMax);
            }
         }
      }

      return factorsAdded;
   }

   // check if all lifted edges are primally consistent by asserting that a path of zero values exists in the ground graph whenever lifted edge is zero
   bool CheckPrimalConsistency(PrimalSolutionStorage::Element primal) const
   {
      const bool multicutConsistent = MULTICUT_CONSTRUCTOR::CheckPrimalConsistency(primal);
      if(!multicutConsistent) {
         return false;
      }
      
      //collect connectivity information with union find w.r.t. base edges
      UnionFind uf(MULTICUT_CONSTRUCTOR::noNodes_);
      for(const auto& e : baseEdges_) {
         if(primal[e.f->GetPrimalOffset()] == false) {
            uf.merge(e.i,e.j);
         }
      }
      for(const auto& e : liftedEdges_) {
        if(primal[e.f->GetPrimalOffset()] == false) {
           if(!uf.connected(e.i,e.j)) {
              return false;
           }
        }
      }
      return true;
   }

   /* primals are different now!
   void ComputePrimal(PrimalSolutionStorage::Element primal) const
   {
      // do zrobienia: templatize for correct graph type
      // do zrobienia: put original graph into multicut constructor and let it be constant, i.e. not reallocate it every time for primal computation. Problem: When adding additional edges, we may not add them to the lifted multicut solver, as extra edges must not participate in cut inequalities

      // use GAEC and Kernighan&Lin algorithm of andres graph package to compute primal solution
      const INDEX noNodes = MULTICUT_CONSTRUCTOR::noNodes_;
      // do zrobienia: possibly build graph only once
      andres::graph::Graph<> originalGraph(noNodes);
      andres::graph::Graph<> liftedGraph(noNodes);
      std::vector<REAL> edgeValues;
      edgeValues.reserve(baseEdges_.size() + liftedEdges_.size());

      std::cout << "# base edges = " << baseEdges_.size() << ", # lifted edges = " << liftedEdges_.size() << "\n";
      // do zrobienia: initalize the graph structures only once
      for(const auto& e : baseEdges_) {
         originalGraph.insertEdge(e.i, e.j);
         liftedGraph.insertEdge(e.i,e.j);
         edgeValues.push_back(e.f->operator[](0));
      }
      for(const auto& e : liftedEdges_) {
         liftedGraph.insertEdge(e.i, e.j);
         edgeValues.push_back(e.f->operator[](0));
      }

      std::vector<char> labeling(baseEdges_.size() + liftedEdges_.size(),0);
      andres::graph::multicut_lifted::greedyAdditiveEdgeContraction(originalGraph,liftedGraph,edgeValues,labeling);
      andres::graph::multicut_lifted::kernighanLin(originalGraph,liftedGraph,edgeValues,labeling,labeling);

      // now write back primal solution and evaluate cost. 
      // Problem: we have to infer state of tightening edges. Note: Multicut is uniquely determined by baseEdges_, hence labeling of thightening edges can be inferred with union find datastructure.
      assert(false); // primals are different now!

      UnionFind uf(noNodes);
      for(INDEX e=0; e<baseEdges_.size(); ++e) {
         std::array<unsigned char,1> l { bool(labeling[e]) };
         baseEdges_[e].f->SetAndPropagatePrimal(primal, l.begin());
         if(labeling[e] == 0) {
            uf.merge(baseEdges_[e].i, baseEdges_[e].j);
         }
      }
      for(INDEX e=0; e<liftedEdges_.size(); ++e) {
         std::array<unsigned char,1> l { (unsigned char)(labeling[e + baseEdges_.size()]) };
         liftedEdges_[e].f->SetAndPropagatePrimal(primal, l.begin());
      }
      // now go over all edges: additionally tightening edges will not have received primal value, set if now.
      for(const auto& e : this->unaryFactors_) {
         const auto* f = e.second;
         if(primal[f->GetPrimalOffset()] == unknownState) { // tightening edge
            std::array<unsigned char,1> l;
            if(uf.connected(e.first[0],e.first[1])) {
               l[0] = 0;
            } else {
               l[0] = 1;
            }
            f->SetAndPropagatePrimal(primal, l.begin());
         }
      }
   }
   */


   private:
   //struct Edge {INDEX i; INDEX j; REAL w;}; // replace by WeightedEdge
   struct MulticutEdge {
      INDEX i; 
      INDEX j; 
      typename MULTICUT_CONSTRUCTOR::UnaryFactorContainer* f;
      REAL weight() const { return (*f->GetFactor()); }
   };
   bool addingTighteningEdges = false; // controls whether edges are added to baseEdges_
   std::vector<MulticutEdge> baseEdges_;
   std::vector<MulticutEdge> liftedEdges_;

   std::vector<std::vector<INDEX>> cutEdgesLiftedMulticutFactors_;
   std::vector<std::vector<INDEX>> liftedEdgesLiftedMulticutFactors_;

   std::map<CutId,std::pair<LiftedMulticutCutFactorContainer*,std::vector<Edge>>> liftedMulticutFactors_;
};

} // end namespace LP_MP

#endif // LP_MP_MULTICUT_CONSTRUCTOR_HXX

