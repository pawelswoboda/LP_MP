#ifndef LP_MP_MULTICUT_CONSTRUCTOR_HXX
#define LP_MP_MULTICUT_CONSTRUCTOR_HXX

#include "multicut_unary_factor.hxx"
#include "multicut_triplet_factor.hxx"
#include "multicut_message.hxx"

#include <unordered_map>
#include <queue>

namespace LP_MP {

// do zrobienia: sort the unary factors such that the original multicut graph is iterated over in some sensible order, i.e. taking into account the node order
// to this end sort factors as follows: first come all unary factors, then all triplet ones.
// let unary factor 1 have edge indices (i1,i2) and factor 2 have edge indices (j1,j2).
// Then (i1,i2) < (j1,j2) if


template<class FACTOR_MESSAGE_CONNECTION, INDEX UNARY_FACTOR_NO, INDEX TRIPLET_FACTOR_NO, INDEX GLOBAL_FACTOR_NO, INDEX UNARY_TRIPLET_MESSAGE_NO, INDEX UNARY_GLOBAL_MESSAGE_NO>
class MulticutConstructor {
protected:
   using FMC = FACTOR_MESSAGE_CONNECTION;

   using UnaryFactorContainer = meta::at_c<typename FMC::FactorList, UNARY_FACTOR_NO>;
   using TripletFactorContainer = meta::at_c<typename FMC::FactorList, TRIPLET_FACTOR_NO>;
   using GlobalFactorContainer = meta::at_c<typename FMC::FactorList, GLOBAL_FACTOR_NO>;
   using UnaryTripletMessageContainer = typename meta::at_c<typename FMC::MessageList, UNARY_TRIPLET_MESSAGE_NO>::MessageContainerType;
   using UnaryTripletMessageType = typename UnaryTripletMessageContainer::MessageType;
   using UnaryGlobalMessageContainer = typename meta::at_c<typename FMC::MessageList, UNARY_GLOBAL_MESSAGE_NO>::MessageContainerType;
   using UnaryGlobalMessageType = typename UnaryGlobalMessageContainer::MessageType;

   // graph for edges with positive cost
   // note: it is not the nicest representation of nodes. head pointers are in the node structure, hence in Dijkstra,w e cannot make the graph const. Store heap pointer in separate structure. Also then concurrent access is possible, hence parallelization would be doable. Same for flag and parent. Each thread would get an additional structure with parent, heap pointer and flag
   // One more option would be to hold the full graph all the time and then making a view on its positive edges.
   struct Arc;
   struct Node {
      Arc* first; // first outgoing arc from node
   };
   struct Arc {
		Node* head;
		Arc* prev; // do zrobienia: not needed currently
		Arc* next;
      Arc* sister; // do zrobienia: not really needed
      REAL cost;
   };
   struct Graph {
      Graph(const INDEX noNodes, const INDEX noEdges) : nodes_(noNodes,{nullptr}), arcs_(2*noEdges), noArcs_(0) {}
      INDEX size() const { return nodes_.size(); }
      // note: this is possibly not very intuitive: we have two ways to index the graph: either by Node* structors or by node numbers, returning the other structure
      // also mixture of references and pointers is ugly
      Node& operator[](const INDEX i) { return nodes_[i]; }
      INDEX operator[](const Node* i) { assert(i - &*nodes_.begin() < nodes_.size()); return i - &*nodes_.begin(); } // do zrobienia: this is ugly pointer arithmetic
      const Node& operator[](const INDEX i) const { return nodes_[i]; }
      const INDEX operator[](const Node* i) const { assert(i - &*nodes_.begin() < nodes_.size()); return i - &*nodes_.begin(); } // do zrobienia: this is ugly pointer arithmetic
      void AddEdge(INDEX i, INDEX j, REAL cost) {
         //std::cout << "Add edge (" << i << "," << j << ") with cost = " << cost << "\n";
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
         //std::cout << "in finding most violated path from " << startNode << " to " << endNode << "\n";
         //PriorityQueue p;
         //p.Add(&g[startNode], std::numeric_limits<REAL>::max());
         AddPQ(startNode, std::numeric_limits<REAL>::max());
         LabelTemporarily(startNode);
         REAL curCost;
         while(!EmptyPQ()) {
            const INDEX i = RemoveMaxPQ(curCost);
            LabelPermanently(i);
            //std::cout << "Investigating node " << g[i] << " with cost " << curCost << "\n";

            // found path
            if(i == endNode) { // trace back path to startNode
               //std::cout << "Found shortest cycle with cost = " << curCost << "\n";
               std::vector<INDEX> path({i});
               INDEX j = i;
               while(Parent(j) != startNode) {
                  j = Parent(j);
                  path.push_back(j);
               }
               path.push_back(startNode);
               //std::cout << "Found cycle ";
               for(INDEX i=0;i<path.size();++i) {
                  //std::cout << path[i] << ",";
               }
               //std::cout << " with violation " << curCost << "\n";
               return std::make_tuple(curCost,path);;
            } 

            for(Arc* a=g[i].first; a!=nullptr; a = a->next) { // do zrobienia: make iteration out of this?
               const REAL edgeCost = a->cost;
               const REAL pathCost = std::min(curCost,edgeCost); // cost of path until head
               Node* head = a->head;
               //std::cout << " Investigating edge to " << g[head] << " with cost = " << edgeCost << "\n";
               if(LabelledPermanently(g[head])) { // permanently labelled
                  //std::cout << "         permanently labelled\n";
                  continue; 
               } else if(LabelledTemporarily(g[head])) { // temporarily labelled
                  //std::cout << "         temporarily labelled, cost = " << p.GetKey(head) << "\n";
                  if(GetKeyPQ(g[head]) < pathCost) {
                     //std::cout << "    Increasing value of node " << g[head] << " to value " << pathCost << "\n";
                     IncreaseKeyPQ(g[head], pathCost);
                     Parent(g[head]) = i;
                  }
               } else { // seen for first time
                  //std::cout << "    Inserting node " << g[head] << " with cost " << pathCost << "\n";
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

   // do zrobienia: one can reuse the parent structure in subsequent FindPath computations
   // However one then needs to store in union find data structures whether end node is in component. Similar can be done in most violated path search, more complicated, though
   struct BfsData {
      struct Item { INDEX parent; INDEX flag; };
      BfsData(const Graph& g) 
      {
         d.resize(g.size());
         for(INDEX i=0; i<d.size(); ++i) {
            d[i].flag = 0;
         }
         flag = 0;
      }
      void Reset() { ++flag; }
      Item& operator[](const INDEX i) { return d[i]; }
      void Label(const INDEX i) { d[i].flag = flag; }
      bool Labelled(const INDEX i) const { return d[i].flag == flag; }
      INDEX& Parent(const INDEX i) { return d[i].parent; }

      // do bfs with thresholded costs and iteratively lower threshold until enough cycles are found
      std::tuple<REAL,std::vector<INDEX>> FindPath(const INDEX startNode, const INDEX endNode, const REAL th, const INDEX maxLength, const Graph& g) 
      {
         Reset();
         std::queue<INDEX> visit; // do zrobienia: do not allocate each time, make visit a member
         visit.push(endNode);
         Label(endNode);
         while(!visit.empty()) {
            const INDEX i=visit.front();
            visit.pop();
            //std::cout << "Investigating node " << i << "\n";

            // found path
            if(i == startNode) { // trace back path to startNode
               //std::cout << "Found shortest cycle\n";
               REAL minPathCost = std::numeric_limits<REAL>::max();
               std::vector<INDEX> path({i});
               INDEX j = i;
               while(Parent(j) != endNode) {
                  //minPathCost = std::min(minPathCost);
                  //std::cout << j << ",";
                  j = Parent(j);
                  path.push_back(j);
               }
               //std::cout << "\n";
               path.push_back(endNode);
               //std::cout << "Found cycle ";
               for(INDEX i=0;i<path.size();++i) {
                  //std::cout << path[i] << ",";
               }
               return std::make_tuple(minPathCost,path);
            } 

            for(Arc* a=g[i].first; a!=nullptr; a = a->next) { // do zrobienia: make iteration out of this?
               if(a->cost >= th) {
                  Node* head = a->head;
                  if(!Labelled(g[head])) {
                     visit.push(g[head]);
                     Parent(g[head]) = i;
                     Label(g[head]);
                  }
               }
            }
         }
         return std::make_tuple(0.0,std::vector<INDEX>(0));
      }




      private:
      std::vector<Item> d;
      INDEX flag;
   };

public:
   MulticutConstructor(ProblemDecomposition<FMC>& pd) : pd_(pd) 
   {
      globalFactor_ = nullptr;
   }
   ~MulticutConstructor()
   {
      static_assert(std::is_same<typename UnaryFactorContainer::FactorType, MulticutUnaryFactor>::value,"");
      static_assert(std::is_same<typename TripletFactorContainer::FactorType, MulticutTripletFactor>::value,"");
      //static_assert(std::is_same<typename MessageContainer::MessageType, MulticutUnaryTripletMessage<MessageSending::SRMP>>::value,"");
   }

   UnaryFactorContainer* AddUnaryFactor(const INDEX i1, const INDEX i2, const REAL cost)
   {
      if(globalFactor_ == nullptr) {
         globalFactor_ = new GlobalFactorContainer(MulticutGlobalFactor(), std::vector<REAL>(0));
         pd_.GetLP()->AddFactor(globalFactor_);
      }
      assert(i1 < i2);
      assert(!HasUnaryFactor(i1,i2));
      
      auto* u = new UnaryFactorContainer(MulticutUnaryFactor(cost), std::vector<REAL>{cost});
      unaryFactors_.insert(std::make_pair( std::make_tuple(i1,i2), u ));
      pd_.GetLP()->AddFactor(u);
      noNodes_ = std::max(noNodes_,std::max(i1,i2)+1);

      LinkUnaryGlobal(u,globalFactor_,i1,i2);
      
      //std::cout << "Add unary factor (" << i1 << "," << i2 << ") with cost = " << cost << "\n";
      return u;
   }
   UnaryFactorContainer* GetUnaryFactor(const INDEX i1, const INDEX i2) const {
      assert(HasUnaryFactor(i1,i2));
      return unaryFactors_.find(std::make_tuple(i1,i2))->second;
   }
   UnaryTripletMessageContainer* LinkUnaryTriplet(UnaryFactorContainer* u, TripletFactorContainer* t, const INDEX i) // argument i denotes which edge the unary factor connects to
   {
      assert(i < 3);
      auto* m = new UnaryTripletMessageContainer(UnaryTripletMessageType(i), u, t, 1);
      pd_.GetLP()->AddMessage(m);
      return m;
   }
   UnaryGlobalMessageContainer* LinkUnaryGlobal(UnaryFactorContainer* u, GlobalFactorContainer* g, const INDEX i1, const INDEX i2)
   {
      auto* m = new UnaryGlobalMessageContainer(UnaryGlobalMessageType( g->GetFactor()->AddEdge(i1,i2) ), u, g, 0);
      pd_.GetLP()->AddMessage(m);
      return m;
   }
   TripletFactorContainer* AddTripletFactor(const INDEX i1, const INDEX i2, const INDEX i3)
   {
      assert(i1 < i2 && i2 < i3);
      assert(!HasTripletFactor(i1,i2,i3));
      assert(HasUnaryFactor(i1,i2) && HasUnaryFactor(i1,i3) && HasUnaryFactor(i2,i3));
      auto* t = new TripletFactorContainer(MulticutTripletFactor(), std::vector<REAL>{0.0,0.0,0.0});
      pd_.GetLP()->AddFactor(t);
      tripletFactors_.insert(std::make_pair( std::make_tuple(i1,i2,i3), t ));
      // get immediate predeccessor and successor and place new triplet in between
      auto succ = tripletFactors_.upper_bound(std::make_tuple(i1,i2,i3));
      if(succ != tripletFactors_.end()) {
         assert(t != succ->second);
         pd_.GetLP()->AddFactorRelation(t,succ->second);
      }
      // link with all three unary factors
      LinkUnaryTriplet(GetUnaryFactor(i1,i2), t, 0);
      LinkUnaryTriplet(GetUnaryFactor(i1,i3), t, 1);
      LinkUnaryTriplet(GetUnaryFactor(i2,i3), t, 2);
      return t;
   }
   bool HasUnaryFactor(const std::tuple<INDEX,INDEX> e) const 
   {
      assert(std::get<0>(e) < std::get<1>(e));
      return unaryFactors_.find(e) != unaryFactors_.end();
   }
   bool HasUnaryFactor(const INDEX i1, const INDEX i2) const 
   {
      //std::cout << "Has Unary factor with: " << i1 << ":" << i2 << "\n";
      assert(i1 < i2);
      return (unaryFactors_.find(std::make_tuple(i1,i2)) != unaryFactors_.end());
   }
   bool HasTripletFactor(const INDEX i1, const INDEX i2, const INDEX i3) const 
   {
      assert(i1 < i2 && i2 < i3);
      return (tripletFactors_.find(std::make_tuple(i1,i2,i3)) != tripletFactors_.end());
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
         //std::cout << "Edge present:  " << i << ", " << (i+1)%cycle.size() << " ; " << cycle[i] << "," << cycle[(i+1)%cycle.size()] << " ; " << std::get<0>(GetEdge(cycle[i], cycle[(i+1)%cycle.size()])) << "," << std::get<1>(GetEdge(cycle[i], cycle[(i+1)%cycle.size()])) << "\n";
         assert(HasUnaryFactor(GetEdge(cycle[i], cycle[(i+1)%cycle.size()])));
      }
      // now we add all triplets with triangulation edge
      INDEX noTripletsAdded = 0;
      for(INDEX i=2; i<cycle.size(); ++i) {
         if(!HasUnaryFactor(minNode, cycle[i])) {
            AddUnaryFactor(minNode,cycle[i],0.0);
         }
         const INDEX secondNode = std::min(cycle[i], cycle[i-1]);
         const INDEX thirdNode = std::max(cycle[i], cycle[i-1]);
         //std::cout << "Add triplet (" << minNode << "," << secondNode << "," << thirdNode << ") -- ";
         if(!HasTripletFactor(minNode, secondNode, thirdNode)) {
            //std::cout << "do so\n"; 
            AddTripletFactor(minNode, secondNode, thirdNode);
            ++noTripletsAdded;
         } else { 
            //std::cout << "already added\n";
         }
      }
      //std::cout << "Added " << noTripletsAdded << " triplet(s)\n";
      return noTripletsAdded;
   }

   // search for cycles to add such that coordinate ascent will be possible
   INDEX Tighten(const REAL minDualIncrease, const INDEX maxCuttingPlanesToAdd)
   {
      return FindNegativeCycles(minDualIncrease,maxCuttingPlanesToAdd);
   }

   // returns number of triplets added
   INDEX FindNegativeCycles(const REAL minDualIncrease, const INDEX maxTripletsToAdd)
   {
      std::vector<std::tuple<INDEX,INDEX,REAL> > negativeEdges;
      Graph posEdgesGraph(noNodes_,unaryFactors_.size()); // graph consisting of positive edges
      for(auto& it : unaryFactors_) {
         const REAL v = it.second->operator[](0);
         const INDEX i = std::get<0>(it.first);
         const INDEX j = std::get<1>(it.first);
         if(v < 0.0) {
            negativeEdges.push_back(std::make_tuple(i,j,v));
         } else {
            posEdgesGraph.AddEdge(i,j,v);
         }
      }

      //std::cout << "Found " << negativeEdges.size() << " negative edges." << std::endl;

      std::sort(negativeEdges.begin(), negativeEdges.end(), [](const std::tuple<INDEX,INDEX,REAL>& e1, const std::tuple<INDEX,INDEX,REAL>& e2)->bool {
            return std::get<2>(e1) > std::get<2>(e2);
            });

      // now search for every negative edge for most negative path from end point to starting point. Do zrobienia: do this in parallel
      // here, longest path is sought after only the edges with positive taken into account
      // the cost of the path is the minimum of the costs of its edges. The guaranteed increase in the dual objective is v_min > -v ? -v : v_min

      //MostViolatedPathData mp(posEdgesGraph);
      BfsData mp(posEdgesGraph);

      INDEX tripletsAdded = 0;
      using CycleType = std::tuple<REAL, std::vector<INDEX>>;
      std::vector< CycleType > cycles;
      for(auto& it : negativeEdges) {
         const INDEX i1 = std::get<0>(it);
         const INDEX i2 = std::get<1>(it);
         const REAL v = std::get<2>(it);
         if(-v < minDualIncrease) break;
         //auto cycle = mp.FindPath(i2,i1,posEdgesGraph);
         auto cycle = mp.FindPath(i2,i1,0.0001,0,posEdgesGraph);
         const REAL dualIncrease = std::min(-v, std::get<0>(cycle));
         if(std::get<1>(cycle).size() > 0) {
            cycles.push_back( std::make_tuple(dualIncrease, std::move(std::get<1>(cycle))) );
         }
      }
      // sort by guaranteed increase in decreasing order
      std::sort(cycles.begin(), cycles.end(), [](const CycleType& i, const CycleType& j) { return std::get<0>(i) > std::get<0>(j); });
      for(auto& cycle : cycles) {
         const REAL cycleDualIncrease = std::get<0>(cycle);
         tripletsAdded += AddCycle(std::get<1>(cycle));
         if(tripletsAdded > maxTripletsToAdd || cycleDualIncrease < minDualIncrease) {
            return tripletsAdded;
         }
      }

      return tripletsAdded;
   }


protected:
   GlobalFactorContainer* globalFactor_;
   // do zrobienia: replace this by unordered_map, provide hash function.
   // possibly dont do this, but use sorting to provide ordering for LP
   std::map<std::tuple<INDEX,INDEX>, UnaryFactorContainer*> unaryFactors_; // actually unary factors in multicut are defined on edges. assume first index < second one
   // sort triplet factors as follows: Let indices be i=(i1,i2,i3) and j=(j1,j2,j3). Then i<j iff i1+i2+i3 < j1+j2+j3 or for ties sort lexicographically
   struct tripletComp {
      bool operator()(const std::tuple<INDEX,INDEX,INDEX> i, const std::tuple<INDEX,INDEX,INDEX> j) const
      {
         const INDEX iSum = std::get<0>(i) + std::get<1>(i) + std::get<2>(i);
         const INDEX jSum = std::get<0>(j) + std::get<1>(j) + std::get<2>(j);
         //if(std::get<0>(i) > std::get<0>(j) && std::get<2>(i) < std::get<2>(j)) return false;
         //else if(std::get<0>(i) < std::get<0>(j) && std::get<2>(i) > std::get<2>(j)) return true;
         //else 
         if(iSum < jSum) return true;
         else if(iSum > jSum) return false;
         else return i<j; // lexicographic comparison
      }
   };
   std::map<std::tuple<INDEX,INDEX,INDEX>, TripletFactorContainer*, tripletComp> tripletFactors_; // triplet factors are defined on cycles of length three
   INDEX noNodes_ = 0;

   ProblemDecomposition<FMC>& pd_;
};

} // end namespace LP_MP

#endif // LP_MP_MULTICUT_CONSTRUCTOR_HXX

