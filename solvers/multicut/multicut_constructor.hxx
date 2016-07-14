#ifndef LP_MP_MULTICUT_CONSTRUCTOR_HXX
#define LP_MP_MULTICUT_CONSTRUCTOR_HXX

#include "LP_MP.h"
#include "multicut_unary_factor.hxx"
#include "multicut_triplet_factor.hxx"
#include "permutation.hxx"
#include "multicut_odd_wheel.hxx"

#include "union_find.hxx"
#include "max_flow.hxx"

#include <unordered_map>
#include <unordered_set>
#include <queue>

#include "andres/graph/graph.hxx"
#include "andres/graph/grid-graph.hxx"
#include "andres/graph/multicut-lifted/kernighan-lin.hxx"
#include "andres/graph/multicut-lifted/greedy-additive.hxx"


namespace LP_MP {


// hash function for maps used in constructors. Do zrobienia: define hash functions used somewhere globally in config.hxx
static auto hf2 = [](const std::array<INDEX,2> x) { return std::hash<INDEX>()(x[0])^std::hash<INDEX>()(x[1]); };
static auto hf3 = [](const std::array<INDEX,3> x) { return std::hash<INDEX>()(x[0])^std::hash<INDEX>()(x[1])^std::hash<INDEX>()(x[2]); };
static auto hf4 = [](const std::array<INDEX,4> x) { return std::hash<INDEX>()(x[0])^std::hash<INDEX>()(x[1])^std::hash<INDEX>()(x[2])^std::hash<INDEX>()(x[3]); };

template<class FACTOR_MESSAGE_CONNECTION, INDEX UNARY_FACTOR_NO, INDEX TRIPLET_FACTOR_NO, INDEX UNARY_TRIPLET_MESSAGE_NO>
class MulticutConstructor {
protected:
   using FMC = FACTOR_MESSAGE_CONNECTION;

   using UnaryFactorContainer = meta::at_c<typename FMC::FactorList, UNARY_FACTOR_NO>;
   using TripletFactorContainer = meta::at_c<typename FMC::FactorList, TRIPLET_FACTOR_NO>;
   //using GlobalFactorContainer = meta::at_c<typename FMC::FactorList, GLOBAL_FACTOR_NO>;
   using UnaryTripletMessageContainer = typename meta::at_c<typename FMC::MessageList, UNARY_TRIPLET_MESSAGE_NO>::MessageContainerType;
   using UnaryTripletMessageType = typename UnaryTripletMessageContainer::MessageType;
   //using UnaryGlobalMessageContainer = typename meta::at_c<typename FMC::MessageList, UNARY_GLOBAL_MESSAGE_NO>::MessageContainerType;
   //using UnaryGlobalMessageType = typename UnaryGlobalMessageContainer::MessageType;

   // graph for edges with positive cost
   struct Arc;
   struct Node {
      Arc* first; // first outgoing arc from node
   };
   struct Arc {
		Node* head;
		Arc* next;
      REAL cost; // not needed for Bfs.
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
      std::tuple<REAL,std::vector<INDEX>> FindPath(const INDEX startNode, const INDEX endNode, const Graph& g) 
      {
         Reset();
         std::queue<INDEX> visit; // do zrobienia: do not allocate each time, make visit a member
         visit.push(startNode);
         Label1(startNode);
         Parent(startNode) = startNode;
         visit.push(endNode);
         Label2(endNode);
         Parent(endNode) = endNode;

         while(!visit.empty()) {
            const INDEX i=visit.front();
            visit.pop();

            if(Labelled1(i)) {
               for(Arc* a=g[i].first; a!=nullptr; a = a->next) { // do zrobienia: make iteration out of this?
                  Node* head = a->head;
                  const INDEX j = g[head];
                  if(!Labelled(j)) {
                     visit.push(j);
                     Parent(j) = i;
                     Label1(j);
                  } else if(Labelled2(j)) { // shortest path found
                     // trace bacj path from j to endNode and from i to startNode
                     std::vector<INDEX> startPath = TracePath(i);
                     std::vector<INDEX> endPath = TracePath(j);
                     std::reverse(endPath.begin(), endPath.end());
                     startPath.insert( startPath.end(), endPath.begin(), endPath.end());
                     return std::make_tuple(0.0, startPath);
                  }
               }
            } else {
               assert(Labelled2(i));
               for(Arc* a=g[i].first; a!=nullptr; a = a->next) { // do zrobienia: make iteration out of this?
                  Node* head = a->head;
                  const INDEX j = g[head];
                  if(!Labelled(j)) {
                     visit.push(j);
                     Parent(j) = i;
                     Label2(j);
                  } else if(Labelled1(j)) { // shortest path found
                     // trace bacj path from j to endNode and from i to startNode
                     std::vector<INDEX> startPath = TracePath(j);
                     std::vector<INDEX> endPath = TracePath(i);
                     std::reverse(endPath.begin(), endPath.end());
                     startPath.insert( startPath.end(), endPath.begin(), endPath.end());
                     return std::make_tuple(0.0, startPath);
                  }
               }
            }
         }
         return std::make_tuple(0.0,std::vector<INDEX>(0));
   }

private:
   std::vector<Item> d;
   INDEX flag1, flag2;
};

public:
MulticutConstructor(Solver<FMC>& pd) 
   : pd_(pd),
   //unaryFactors_(100,hf2),
   tripletFactors_(100,hf3)
   {
      //unaryFactors_.max_load_factor(0.7);
      tripletFactors_.max_load_factor(0.7);
   }
   ~MulticutConstructor()
   {
      static_assert(std::is_same<typename UnaryFactorContainer::FactorType, MulticutUnaryFactor>::value,"");
      static_assert(std::is_same<typename TripletFactorContainer::FactorType, MulticutTripletFactor>::value,"");
      //static_assert(std::is_same<typename MessageContainer::MessageType, MulticutUnaryTripletMessage<MessageSending::SRMP>>::value,"");
   }

   virtual UnaryFactorContainer* AddUnaryFactor(const INDEX i1, const INDEX i2, const REAL cost) // declared virtual so that derived class notices when unary factor is added
   {
      //if(globalFactor_ == nullptr) {
         //globalFactor_ = new GlobalFactorContainer(MulticutGlobalFactor(), 1); // we have one element currently
      //   pd_.GetLP().AddFactor(globalFactor_);
      //} else {
      //   globalFactor_->ResizeRepam(unaryFactors_.size()+1);
      //}
      assert(i1 < i2);
      assert(!HasUnaryFactor(i1,i2));
      
      auto* u = new UnaryFactorContainer(MulticutUnaryFactor(cost), std::vector<REAL>{cost});
      pd_.GetLP().AddFactor(u);
      auto it = unaryFactors_.insert(std::make_pair(std::array<INDEX,2>{i1,i2}, u)).first;
      //std::cout << "current edge: (" << i1 << "," << i2 << ")";
      if(it != unaryFactors_.begin()) {
         auto prevIt = it;
         --prevIt;
         //std::cout << ", prev edge: (" << prevIt->first.operator[](0) << "," << prevIt->first.operator[](1) << ")";
         assert(prevIt->second != u);
         pd_.GetLP().AddFactorRelation(prevIt->second, u);
      }
      auto nextIt = it;
      ++nextIt;
      if(nextIt != unaryFactors_.end()) {
         assert(nextIt->second != u);
         //std::cout << ", next edge: (" << nextIt->first.operator[](0) << "," << nextIt->first.operator[](1) << ")";
         pd_.GetLP().AddFactorRelation(u, nextIt->second);
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
      auto* m = new UnaryTripletMessageContainer(UnaryTripletMessageType(i), u, t, UnaryTripletMessageType::size());
      pd_.GetLP().AddMessage(m);
      return m;
   }
   //UnaryGlobalMessageContainer* LinkUnaryGlobal(UnaryFactorContainer* u, GlobalFactorContainer* g, const INDEX i1, const INDEX i2)
   //{
   //   auto* m = new UnaryGlobalMessageContainer(UnaryGlobalMessageType( g->GetFactor()->AddEdge(i1,i2) ), u, g, 0);
   //   pd_.GetLP().AddMessage(m);
   //   return m;
   //}
   virtual TripletFactorContainer* AddTripletFactor(const INDEX i1, const INDEX i2, const INDEX i3) // declared virtual so that derived constructor notices when triplet factor is added
   {
      assert(i1 < i2 && i2 < i3);
      assert(!HasTripletFactor(i1,i2,i3));
      assert(HasUnaryFactor(i1,i2) && HasUnaryFactor(i1,i3) && HasUnaryFactor(i2,i3));
      auto* t = new TripletFactorContainer(MulticutTripletFactor(), std::vector<REAL>(MulticutTripletFactor::size(),0.0));
      pd_.GetLP().AddFactor(t);
      tripletFactors_.insert(std::make_pair( std::array<INDEX,3>{i1,i2,i3}, t ));
      // use following ordering of unary and triplet factors: triplet comes after edge factor (i1,i2) and before (i2,i3)
      auto* before = GetUnaryFactor(i1,i2);
      pd_.GetLP().AddFactorRelation(before,t);
      auto* middle = GetUnaryFactor(i1,i3);
      pd_.GetLP().AddFactorRelation(middle,t);
      auto* after = GetUnaryFactor(i2,i3);
      pd_.GetLP().AddFactorRelation(t,after);
      // get immediate predeccessor and successor and place new triplet in between
      //auto succ = tripletFactors_.upper_bound(std::make_tuple(i1,i2,i3));
      //if(succ != tripletFactors_.end()) {
      //   assert(t != succ->second);
      //   pd_.GetLP().AddFactorRelation(t,succ->second);
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
         //logger->info() << "Edge present:  " << i << ", " << (i+1)%cycle.size() << " ; " << cycle[i] << "," << cycle[(i+1)%cycle.size()] << " ; " << std::get<0>(GetEdge(cycle[i], cycle[(i+1)%cycle.size()])) << "," << std::get<1>(GetEdge(cycle[i], cycle[(i+1)%cycle.size()]));
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
         //logger->info() << "Add triplet (" << minNode << "," << secondNode << "," << thirdNode << ") -- ";
         if(!HasTripletFactor(minNode, secondNode, thirdNode)) {
            //logger->info() << "do so"; 
            AddTripletFactor(minNode, secondNode, thirdNode);
            ++noTripletsAdded;
         } else { 
            //logger->info() << "already added";
         }
      }
      //logger->info() << "Added " << noTripletsAdded << " triplet(s)";
      return noTripletsAdded;
   }

   // search for cycles to add such that coordinate ascent will be possible
   INDEX Tighten(const REAL minDualIncrease, const INDEX maxCuttingPlanesToAdd)
   {
      const INDEX tripletsAdded = FindNegativeCycles(minDualIncrease,maxCuttingPlanesToAdd);
      spdlog::get("logger")->info() << "Added " << tripletsAdded << " triplet(s)";
      return tripletsAdded;
   }

   // returns number of triplets added
   INDEX FindNegativeCycles(const REAL minDualIncrease, const INDEX maxTripletsToAdd)
   {
      std::vector<std::tuple<INDEX,INDEX,REAL> > negativeEdges;
      // we can speed up compution by skipping path searches for node pairs which lie in different connected components. Connectedness is stored in a union find structure
      UnionFind uf(noNodes_);
      Graph posEdgesGraph(noNodes_,unaryFactors_.size()); // graph consisting of positive edges
      for(auto& it : unaryFactors_) {
         const REAL v = it.second->operator[](0);
         const INDEX i = std::get<0>(it.first);
         const INDEX j = std::get<1>(it.first);
         if(v > minDualIncrease) {
            posEdgesGraph.AddEdge(i,j,v);
            uf.merge(i,j);
         }
      }
      for(auto& it : unaryFactors_) {
         const REAL v = it.second->operator[](0);
         const INDEX i = std::get<0>(it.first);
         const INDEX j = std::get<1>(it.first);
         if(v <= -minDualIncrease && uf.connected(i,j)) {
            negativeEdges.push_back(std::make_tuple(i,j,v));
         }
      }
      // do zrobienia: possibly add reparametrization of triplet factors additionally

      //logger->info() << "Found " << negativeEdges.size() << " negative edges." << std::endl;

      std::sort(negativeEdges.begin(), negativeEdges.end(), [](const std::tuple<INDEX,INDEX,REAL>& e1, const std::tuple<INDEX,INDEX,REAL>& e2)->bool {
            return std::get<2>(e1) > std::get<2>(e2);
            });

      // now search for every negative edge for most negative path from end point to starting point. Do zrobienia: do this in parallel
      // here, longest path is sought after only the edges with positive taken into account
      // the cost of the path is the minimum of the costs of its edges. The guaranteed increase in the dual objective is v_min > -v ? -v : v_min

      //MostViolatedPathData mp(posEdgesGraph);
      BfsData mp(posEdgesGraph);

      // better: collect all cycles, then sort them according to violation, and only then add them






      INDEX tripletsAdded = 0;
      using CycleType = std::tuple<REAL, std::vector<INDEX>>;
      std::vector< CycleType > cycles;
      for(auto& it : negativeEdges) {
         const INDEX i = std::get<0>(it);
         const INDEX j = std::get<1>(it);
         const REAL v = std::get<2>(it);
         if(-v < minDualIncrease) break;
         auto cycle = mp.FindPath(i,j,posEdgesGraph);
         const REAL dualIncrease = std::min(-v, std::get<0>(cycle));
         assert(std::get<1>(cycle).size() > 0);
         if(std::get<1>(cycle).size() > 0) {
            //cycles.push_back( std::make_tuple(dualIncrease, std::move(std::get<1>(cycle))) );
            tripletsAdded += AddCycle(std::get<1>(cycle));
            if(tripletsAdded > maxTripletsToAdd) {
               return tripletsAdded;
            }
         } else {
            throw std::runtime_error("No path found although there should be one"); 
         }
      }
      // sort by guaranteed increase in decreasing order
      /*
      std::sort(cycles.begin(), cycles.end(), [](const CycleType& i, const CycleType& j) { return std::get<0>(i) > std::get<0>(j); });
      for(auto& cycle : cycles) {
         const REAL cycleDualIncrease = std::get<0>(cycle);
         tripletsAdded += AddCycle(std::get<1>(cycle));
         if(tripletsAdded > maxTripletsToAdd) {
            return tripletsAdded;
         }
      }
      */

      return tripletsAdded;
   }

   bool CheckPrimalConsistency(PrimalSolutionStorage::Element primal) const
   {
      //std::cout << "checking primal feasibility for multicut\n";
      UnionFind uf(noNodes_);
      for(const auto& e : unaryFactors_) {
         UnaryFactorContainer* f = e.second; 
         assert(primal[f->GetPrimalOffset()] != unknownState);
         if(primal[f->GetPrimalOffset()] == false) {
            // connect components 
            const INDEX i = std::get<0>(e.first);
            const INDEX j = std::get<1>(e.first);
            uf.merge(i,j);
         }
      }
      for(const auto& e : unaryFactors_) {
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

protected:
   //GlobalFactorContainer* globalFactor_;
   // do zrobienia: replace this by unordered_map, provide hash function.
   // possibly dont do this, but use sorting to provide ordering for LP
   std::map<std::array<INDEX,2>, UnaryFactorContainer*> unaryFactors_; // actually unary factors in multicut are defined on edges. assume first index < second one
   //std::unordered_map<std::array<INDEX,2>, UnaryFactorContainer*, decltype(hf2)> unaryFactors_; // actually unary factors in multicut are defined on edges. assume first index < second one
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
   std::unordered_map<std::array<INDEX,3>, TripletFactorContainer*, decltype(hf3)> tripletFactors_; // triplet factors are defined on cycles of length three
   INDEX noNodes_ = 0;

   Solver<FMC>& pd_;
};



template<class FACTOR_MESSAGE_CONNECTION, INDEX UNARY_FACTOR_NO, INDEX TRIPLET_FACTOR_NO, INDEX UNARY_TRIPLET_MESSAGE_NO,
   INDEX TRIPLET_PLUS_SPOKE_FACTOR_NO, INDEX TRIPLET_PLUS_SPOKE_MESSAGE_NO, INDEX TRIPLET_PLUS_SPOKE_COVER_MESSAGE_NO>
class MulticutOddWheelConstructor : public MulticutConstructor<FACTOR_MESSAGE_CONNECTION, UNARY_FACTOR_NO, TRIPLET_FACTOR_NO, UNARY_TRIPLET_MESSAGE_NO> {
   using FMC = FACTOR_MESSAGE_CONNECTION;
   using BaseConstructor = MulticutConstructor<FACTOR_MESSAGE_CONNECTION, UNARY_FACTOR_NO, TRIPLET_FACTOR_NO, UNARY_TRIPLET_MESSAGE_NO>;

   using TripletPlusSpokeFactorContainer = meta::at_c<typename FMC::FactorList, TRIPLET_PLUS_SPOKE_FACTOR_NO>;
   using TripletPlusSpokeMessageContainer = typename meta::at_c<typename FMC::MessageList, TRIPLET_PLUS_SPOKE_MESSAGE_NO>::MessageContainerType;
   using TripletPlusSpokeCoverMessageContainer = typename meta::at_c<typename FMC::MessageList, TRIPLET_PLUS_SPOKE_COVER_MESSAGE_NO>::MessageContainerType;
public:
   MulticutOddWheelConstructor(Solver<FMC>& pd) : BaseConstructor(pd), tripletPlusSpokeFactors_(100,hf4) { 
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
      auto* tps = new TripletPlusSpokeFactorContainer(MulticutTripletPlusSpokeFactor(),std::vector<REAL>(MulticutTripletPlusSpokeFactor::size(),0.0));
      BaseConstructor::pd_.GetLP().AddFactor(tps);
      assert(n1<n2);
      tripletPlusSpokeFactors_.insert(std::make_pair(std::array<INDEX,4>({n1,n2,centerNode,spokeNode}),tps));
      std::array<INDEX,3> tripletIndices{n1,n2,centerNode};
      std::sort(tripletIndices.begin(), tripletIndices.end());
      auto* t = BaseConstructor::GetTripletFactor(tripletIndices[0], tripletIndices[1], tripletIndices[2]);
      auto* m = new TripletPlusSpokeCoverMessageContainer(MulticutTripletPlusSpokeCoverMessage(n1,n2,centerNode,spokeNode), t, tps, MulticutTripletPlusSpokeCoverMessage::size());
      BaseConstructor::pd_.GetLP().AddMessage(m);
      //BaseConstructor::pd_.GetLP().AddFactorRelation(t,tps);
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

      auto* m = new TripletPlusSpokeMessageContainer(MulticutTripletPlusSpokeMessage(n1,n2,centerNode,spokeNode,tripletIndices[0],tripletIndices[1],tripletIndices[2]),t,tps,MulticutTripletPlusSpokeMessage::size());
      //auto* m = new TripletPlusSpokeMessageContainer(MulticutTripletPlusSpokeMessage(sharedTripletEdgeTriplet,spokeEdgeTriplet,sharedTripletEdgeTripletPlusSpoke),t,tps,MulticutTripletPlusSpokeMessage::size());
      BaseConstructor::pd_.GetLP().AddMessage(m);
      if(tripletNode == n1) {
         //BaseConstructor::pd_.GetLP().AddFactorRelation(t,tps);
      } else {
         assert(tripletNode == n2);
         //BaseConstructor::pd_.GetLP().AddFactorRelation(tps,t);
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

   INDEX FindOddWheels(const REAL minDualIncrease, const INDEX maxCuttingPlanesToAdd)
   {
      INDEX oddWheelsAdded = 0;

      // search for all triangles present in the graph, also in places where no triplet factor has been added
      // Construct adjacency list
      std::vector<std::vector<INDEX>> adjacencyList(BaseConstructor::noNodes_);
      for(auto& e : BaseConstructor::unaryFactors_) {
         const INDEX i = std::get<0>(e)[0];
         const INDEX j = std::get<0>(e)[1];
         assert(i<j);
         adjacencyList[i].push_back(j);
         adjacencyList[j].push_back(i);
      }
      // sort it for fast intersection
      for(INDEX i=0; i<adjacencyList.size(); ++i) {
         std::sort(adjacencyList[i].begin(), adjacencyList[i].end());
      }
      /*
      // count triangles
      INDEX noTriangles = 0;
      for(auto& e : BaseConstructor::unaryFactors_) {
         const INDEX i = std::get<0>(e)[0];
         const INDEX j = std::get<0>(e)[1];
         assert(i<j);
         auto commonNodesEnd = std::set_intersection(adjacencyList[i].begin(), adjacencyList[i].end(), adjacencyList[j].begin(), adjacencyList[j].end(), commonNodes.begin());
         for(auto it=commonNodes.begin(); it!=commonNodesEnd; ++it) {
            if(j < *it) { // triangles are counted multiple times
               ++noTriangles;
            }
         }
      }
      std::cout << "no triangles = " << noTriangles << "\n";
      exit(1);
      // construct all triangles
      std::unordered_set<std::array<INDEX,3>,decltype(hf3)> triangles {noTriangles,hf3};
      for(auto& e : BaseConstructor::unaryFactors_) {
         const INDEX i = std::get<0>(e)[0];
         const INDEX j = std::get<0>(e)[1];
         assert(i<j);
         auto commonNodesEnd = std::set_intersection(adjacencyList[i].begin(), adjacencyList[i].end(), adjacencyList[j].begin(), adjacencyList[j].end(), commonNodes.begin());
         for(auto it=commonNodes.begin(); it!=commonNodesEnd; ++it) {
            if(j < *it) { // triangles are counted multiple times
               triangles.insert(std::array<INDEX,3>{i,j,*it});
            }
         }
      }
      */


      std::vector<INDEX> commonNodes(BaseConstructor::noNodes_); // for detecting triangles
      for(INDEX i=0; i<BaseConstructor::noNodes_; ++i) {
         // compressed nodes for later usage in bfs
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
                  if(BaseConstructor::HasTripletFactor(triplet[0],triplet[1],triplet[2])) {
                     auto* t = BaseConstructor::GetTripletFactor(triplet[0],triplet[1],triplet[2]);
                     INDEX l1, l2;
                     INDEX l3, l4; // the other labelings
                     assert(t->GetFactor()->size() == 4);
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
                     const REAL ij = BaseConstructor::GetUnaryFactor(std::min(i,j), std::max(i,j))->operator[](0);
                     const REAL ik = BaseConstructor::GetUnaryFactor(std::min(i,k), std::max(i,k))->operator[](0);
                     const REAL jk = BaseConstructor::GetUnaryFactor(j,k)->operator[](0);
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
            typename BaseConstructor::TripletFactorContainer*  t = std::get<2>(tripletByIndices_[i][j]);
            const INDEX i1 = std::get<0>(tripletByIndices_[i][j]);
            const INDEX i2 = std::get<1>(tripletByIndices_[i][j]);
            // check if 110 or 101 are among the minimal labelings, where the first edge is the outer cycle edge
            // the corresponding labeling numbers are
            INDEX l1, l2;
            INDEX l3, l4; // the other labelings
            assert(t->GetFactor()->size() == 4);
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
   }

   INDEX Tighten(const REAL minDualIncrease, const INDEX maxCuttingPlanesToAdd)
   {
      const INDEX tripletsAdded = BaseConstructor::Tighten(minDualIncrease, maxCuttingPlanesToAdd);
      if(tripletsAdded >= maxCuttingPlanesToAdd ) {
         return tripletsAdded;
      } else {
         // require the odd wheels to have larger impact relative to violated cycles. This ensures that odd wheels are only added late in the optimziation, when no good violated cycles are present any more
         const INDEX oddWheelsAdded = FindOddWheels(1e10*minDualIncrease, maxCuttingPlanesToAdd - tripletsAdded);
         spdlog::get("logger")->info() << "Added " << oddWheelsAdded << " factors for odd wheel constraints";
         return tripletsAdded + oddWheelsAdded;
      }
      /*
      if(minDualIncrease <= 1e-14) {
         const INDEX oddWheelsAdded = FindOddWheels(minDualIncrease*1e12, maxCuttingPlanesToAdd);
         logger->info() << "Added " << oddWheelsAdded << " factors for odd wheel constraints";
         return tripletsAdded + oddWheelsAdded;
      }
      return tripletsAdded;
      */

      assert(false);
   }


private:
   std::vector<std::vector<std::tuple<INDEX,INDEX,typename BaseConstructor::TripletFactorContainer*>>> tripletByIndices_; // of triplet factor with indices (i1,i2,i3) exists, then (i1,i2,i3) will be in the vector of index i1, i2 and i3
   // the format for TripletPlusSpoke is (node1,node2, centerNode, spokeNode) and we assume n1<n2
   // hash for std::array<INDEX,4>
   std::unordered_map<std::array<INDEX,4>,TripletPlusSpokeFactorContainer*,decltype(hf4)> tripletPlusSpokeFactors_;
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

   LiftedMulticutConstructor(Solver<FMC>& pd) : MULTICUT_CONSTRUCTOR(pd) {}

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

   // do zrobienia: probide AddCutFactor(const CutId& cut, const INDEX i1, const INDEX i2) as well
   LiftedMulticutCutFactorContainer* AddCutFactor(const CutId& cut)
   {
      assert(!HasCutFactor(cut));
      auto* f = new LiftedMulticutCutFactorContainer(LiftedMulticutCutFactor(cut.size()),std::vector<REAL>(cut.size(),0));
      MULTICUT_CONSTRUCTOR::pd_.GetLP().AddFactor(f);
      // connect the cut edges
      for(INDEX e=0; e<cut.size(); ++e) {
         auto* unaryFactor = MULTICUT_CONSTRUCTOR::GetUnaryFactor(cut[e][0],cut[e][1]);
         auto* m = new CutEdgeLiftedMulticutFactorMessageContainer(CutEdgeLiftedMulticutFactorMessage(e),unaryFactor,f,1); // do zrobienia: remove 1
         MULTICUT_CONSTRUCTOR::pd_.GetLP().AddMessage(m);
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
      assert(i1<i2);
      auto& c = liftedMulticutFactors_[cut];
      auto* f = c.first;
      auto& edgeList = c.second;
      auto* unaryFactor = MULTICUT_CONSTRUCTOR::GetUnaryFactor(i1,i2);
      auto* m = new LiftedEdgeLiftedMulticutFactorMessageContainer(LiftedEdgeLiftedMulticutFactorMessage(edgeList.size() + cut.size()), unaryFactor, f, 1); // do zrobienia: remove 1
      MULTICUT_CONSTRUCTOR::pd_.GetLP().AddMessage(m);
      c.second.push_back(Edge({i1,i2}));
      f->resize(f->size()+1,0.0);
      f->GetFactor()->IncreaseLifted();
   }



   INDEX Tighten(const REAL minDualIncrease, const INDEX maxCuttingPlanesToAdd)
   {
      const bool prevMode = addingTighteningEdges;
      addingTighteningEdges = true;
      assert(maxCuttingPlanesToAdd > 5); //otherwise the below arrangement makes little sense.
      const INDEX noBaseConstraints = MULTICUT_CONSTRUCTOR::Tighten(minDualIncrease, 0.8*maxCuttingPlanesToAdd);
      INDEX noLiftingConstraints = 0;
      spdlog::get("logger")->info() << "number of cut constraints: " << liftedMulticutFactors_.size();
      //for(INDEX i=0; i<baseEdges_.size(); ++i) {
      //   std::cout << baseEdges_[i].weight() << ", ";
      //}
      //std::cout << "\n";
      //for(INDEX i=0; i<liftedEdges_.size(); ++i) {
      //   std::cout << liftedEdges_[i].weight() << ", ";
      //}
      std::cout << "\n";
      if(noBaseConstraints < maxCuttingPlanesToAdd) {
         noLiftingConstraints = FindViolatedCuts(minDualIncrease, maxCuttingPlanesToAdd - noBaseConstraints);
         spdlog::get("logger")->info() << "added " << noLiftingConstraints << " lifted cut factors";
      }
      addingTighteningEdges = prevMode;
      return noBaseConstraints + noLiftingConstraints;
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

   void ComputePrimal(PrimalSolutionStorage::Element primal) const
   {
      // do zrobienia: templatize for correct graph type
      // do zrobienia: put original graph into multicut constructor and let it be constant, i.e. not reallocate it every time for primal computation. Problem: When adding additional edges, we may not add them to the lifted multicut solver, as extra edges must not participate in cut inequalities

      // use GAEC and Kernighan&Lin algorithm of andres graph package to compute primal solution
      const INDEX noNodes = MULTICUT_CONSTRUCTOR::noNodes_;
      andres::graph::Graph<> originalGraph(noNodes);
      andres::graph::Graph<> liftedGraph(noNodes);
      std::vector<REAL> edgeValues;
      edgeValues.reserve(baseEdges_.size() + liftedEdges_.size());

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
      UnionFind uf(noNodes);
      for(INDEX e=0; e<baseEdges_.size(); ++e) {
         std::array<bool,1> l { bool(labeling[e]) };
         baseEdges_[e].f->SetAndPropagatePrimal(primal, l.begin());
         uf.merge(baseEdges_[e].i, baseEdges_[e].j);
      }
      for(INDEX e=0; e<liftedEdges_.size(); ++e) {
         std::array<bool,1> l { bool(labeling[e + baseEdges_.size()]) };
         liftedEdges_[e].f->SetAndPropagatePrimal(primal, l.begin());
      }
      // now go over all edges: additionally tightening edges will not have received primal value, set if now.
      for(const auto& e : MULTICUT_CONSTRUCTOR::unaryFactors_) {
         const auto* f = e.second;
         if(primal[f->GetPrimalOffset()] != unknownState) {
            std::array<bool,1> l;
            if(uf.connected(e.first[0],e.first[1])) {
               l[0] = 0;
            } else {
               l[0] = 1;
            }
            f->SetAndPropagatePrimal(primal, l.begin());
         }
      }
   }


   private:
   //struct Edge {INDEX i; INDEX j; REAL w;}; // replace by WeightedEdge
   struct MulticutEdge {
      INDEX i; 
      INDEX j; 
      typename MULTICUT_CONSTRUCTOR::UnaryFactorContainer* f;
      REAL weight() const { return (*f)[0]; }
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

