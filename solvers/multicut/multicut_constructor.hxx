#ifndef LP_MP_MULTICUT_CONSTRUCTOR_HXX
#define LP_MP_MULTICUT_CONSTRUCTOR_HXX

#include "LP_MP.h"
#include "multicut_unary_factor.hxx"
#include "multicut_triplet_factor.hxx"
//#include "multicut_message.hxx"
#include "permutation.hxx"
#include "multicut_odd_wheel.hxx"

#include <unordered_map>
#include <queue>

namespace LP_MP {

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
		      //std::cout << j << ",";
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

   virtual UnaryFactorContainer* AddUnaryFactor(const INDEX i1, const INDEX i2, const REAL cost) // declared virtual so that derived class notices when unary factor is added
   {
      if(globalFactor_ == nullptr) {
         globalFactor_ = new GlobalFactorContainer(MulticutGlobalFactor(), 1); // we have one element currently
         pd_.GetLP()->AddFactor(globalFactor_);
      } else {
         globalFactor_->ResizeRepam(unaryFactors_.size()+1);
      }
      assert(i1 < i2);
      assert(!HasUnaryFactor(i1,i2));
      
      auto* u = new UnaryFactorContainer(MulticutUnaryFactor(cost), std::vector<REAL>{cost});
      unaryFactors_.insert(std::make_pair(std::make_tuple(i1,i2), u));
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
      auto* m = new UnaryTripletMessageContainer(UnaryTripletMessageType(i), u, t, UnaryTripletMessageType::size());
      pd_.GetLP()->AddMessage(m);
      return m;
   }
   UnaryGlobalMessageContainer* LinkUnaryGlobal(UnaryFactorContainer* u, GlobalFactorContainer* g, const INDEX i1, const INDEX i2)
   {
      auto* m = new UnaryGlobalMessageContainer(UnaryGlobalMessageType( g->GetFactor()->AddEdge(i1,i2) ), u, g, 0);
      pd_.GetLP()->AddMessage(m);
      return m;
   }
   virtual TripletFactorContainer* AddTripletFactor(const INDEX i1, const INDEX i2, const INDEX i3) // declared virtual so that derived constructor notices when triplet factor is added
   {
      assert(i1 < i2 && i2 < i3);
      assert(!HasTripletFactor(i1,i2,i3));
      assert(HasUnaryFactor(i1,i2) && HasUnaryFactor(i1,i3) && HasUnaryFactor(i2,i3));
      auto* t = new TripletFactorContainer(MulticutTripletFactor(), std::vector<REAL>(MulticutTripletFactor::size(),0.0));
      pd_.GetLP()->AddFactor(t);
      tripletFactors_.insert(std::make_pair( std::make_tuple(i1,i2,i3), t ));
      // get immediate predeccessor and successor and place new triplet in between
      auto succ = tripletFactors_.upper_bound(std::make_tuple(i1,i2,i3));
      if(succ != tripletFactors_.end()) {
         assert(t != succ->second);
         pd_.GetLP()->AddFactorRelation(t,succ->second);
      }
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

   TripletFactorContainer* GetTripletFactor(const INDEX i1, const INDEX i2, const INDEX i3) const 
   {
      assert(HasTripletFactor(i1,i2,i3));
      return tripletFactors_.find(std::make_tuple(i1,i2,i3))->second;
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
      // now we add all triplets with triangulation edge. Possibly, a better triangulation scheme would be possible
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
            throw std::runtime_error("No path found although there is one"); 
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
         if(iSum < jSum) return true;
         else if(iSum > jSum) return false;
         else return i<j; // lexicographic comparison
      }
   };
   std::map<std::tuple<INDEX,INDEX,INDEX>, TripletFactorContainer*, tripletComp> tripletFactors_; // triplet factors are defined on cycles of length three
   INDEX noNodes_ = 0;

   ProblemDecomposition<FMC>& pd_;
};

template<class FACTOR_MESSAGE_CONNECTION, INDEX UNARY_FACTOR_NO, INDEX TRIPLET_FACTOR_NO, INDEX GLOBAL_FACTOR_NO, INDEX UNARY_TRIPLET_MESSAGE_NO, INDEX UNARY_GLOBAL_MESSAGE_NO,
   //INDEX ODD_WHEEL_FACTOR_NO, INDEX TRIPLET_ODD_WHEEL_CYCLE_MESSAGE_NO>
   //INDEX CRITICAL_TRIPLET_FACTOR_NO, INDEX CRITICAL_TRIPLET_COUNTING_FACTOR_NO, INDEX CRITICAL_TRIPLET_MESSAGE_NO, INDEX CRITICAL_TRIPLET_COUNTING_MESSAGE_NO>
   INDEX TRIPLET_PLUS_SPOKE_FACTOR_NO, INDEX TRIPLET_PLUS_SPOKE_MESSAGE_NO, INDEX TRIPLET_PLUS_SPOKE_COVER_MESSAGE_NO>
class MulticutOddWheelConstructor : public MulticutConstructor<FACTOR_MESSAGE_CONNECTION, UNARY_FACTOR_NO, TRIPLET_FACTOR_NO, GLOBAL_FACTOR_NO, UNARY_TRIPLET_MESSAGE_NO, UNARY_GLOBAL_MESSAGE_NO> {
   using FMC = FACTOR_MESSAGE_CONNECTION;
   using BaseConstructor = MulticutConstructor<FACTOR_MESSAGE_CONNECTION, UNARY_FACTOR_NO, TRIPLET_FACTOR_NO, GLOBAL_FACTOR_NO, UNARY_TRIPLET_MESSAGE_NO, UNARY_GLOBAL_MESSAGE_NO>;
   //using OddWheelFactorContainer = meta::at_c<typename FMC::FactorList, ODD_WHEEL_FACTOR_NO>;
   //using TripletOddWheelMessageContainer = typename meta::at_c<typename FMC::MessageList, TRIPLET_ODD_WHEEL_CYCLE_MESSAGE_NO>::MessageContainerType;
   //using TripletOddWheelMessageType = typename TripletOddWheelMessageContainer::MessageType;

   //using CriticalTripletFactorContainer = meta::at_c<typename FMC::FactorList, CRITICAL_TRIPLET_FACTOR_NO>;
   using TripletPlusSpokeFactorContainer = meta::at_c<typename FMC::FactorList, TRIPLET_PLUS_SPOKE_FACTOR_NO>;
   //using CriticalTripletMessageContainer = typename meta::at_c<typename FMC::MessageList, CRITICAL_TRIPLET_MESSAGE_NO>::MessageContainerType;
   using TripletPlusSpokeMessageContainer = typename meta::at_c<typename FMC::MessageList, TRIPLET_PLUS_SPOKE_MESSAGE_NO>::MessageContainerType;
   using TripletPlusSpokeCoverMessageContainer = typename meta::at_c<typename FMC::MessageList, TRIPLET_PLUS_SPOKE_COVER_MESSAGE_NO>::MessageContainerType;
   /*
   using UnaryOddWheelCycleMessageContainer = typename meta::at_c<typename FMC::MessageList, UNARY_ODD_WHEEL_CYCLE_MESSAGE_NO>::MessageContainerType;
   using UnaryOddWheelCycleMessageType = typename UnaryOddWheelCycleMessageContainer::MessageType;
   using UnaryOddWheelCenterMessageContainer = typename meta::at_c<typename FMC::MessageList, UNARY_ODD_WHEEL_CENTER_MESSAGE_NO>::MessageContainerType;
   using UnaryOddWheelCenterMessageType = typename UnaryOddWheelCenterMessageContainer::MessageType;
   */
public:
   MulticutOddWheelConstructor(ProblemDecomposition<FACTOR_MESSAGE_CONNECTION>& pd) : BaseConstructor(pd) {}

   // add triplet indices additionally to tripletIndices_
   typename BaseConstructor::UnaryFactorContainer* AddUnaryFactor(const INDEX i1, const INDEX i2, const REAL cost) final
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
   typename BaseConstructor::TripletFactorContainer* AddTripletFactor(const INDEX i1, const INDEX i2, const INDEX i3) final
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

   /*
   bool HasOddWheelFactor(const INDEX centerNode, const std::vector<INDEX>& cycle) const
   {
      // expects cycle to be in normal form
      //assert(CycleNormalForm(
      return (oddWheelFactors_.find(std::make_tuple(centerNode,cycle)) != oddWheelFactors_.end());
   }
   */
   /*
   TripletOddWheelMessageContainer* LinkTripletOddWheel(typename BaseConstructor::TripletFactorContainer* t, OddWheelFactorContainer* o, std::array<INDEX,3> indices) 
   {
      auto* m = new TripletOddWheelMessageContainer(TripletOddWheelMessageType(indices), t, o, 3);
      BaseConstructor::pd_.GetLP()->AddMessage(m);
      return m;
   }
   */

   /*
   CriticalTripletCountingMessageContainer* LinkCountingTriplet(typename BaseConstructor::TripletFactorContainer* t, CriticalTripletCountingFactorContainer* c, const INDEX edgeNo, const INDEX whichEdgeIsCycle)
   {
      auto* m = new CriticalTripletCountingMessageContainer(MulticutCriticalTripletCountingMessage(edgeNo,whichEdgeIsCycle),t,c,1);
      BaseConstructor::pd_.GetLP()->AddMessage(m);
      return m;
   }
   */
   TripletPlusSpokeFactorContainer* AddTripletPlusSpokeFactor(const INDEX n1, const INDEX n2, const INDEX centerNode, const INDEX spokeNode)
   {
      assert(!HasTripletPlusSpokeFactor(n1,n2,centerNode,spokeNode));
      auto* tps = new TripletPlusSpokeFactorContainer(MulticutTripletPlusSpokeFactor(),std::vector<REAL>(MulticutTripletPlusSpokeFactor::size(),0.0));
      assert(n1<n2);
      tripletPlusSpokeFactors_.insert(std::make_pair(std::array<INDEX,4>({n1,n2,centerNode,spokeNode}),tps));
      BaseConstructor::pd_.GetLP()->AddFactor(tps);
      std::array<INDEX,3> tripletIndices{n1,n2,centerNode};
      std::sort(tripletIndices.begin(), tripletIndices.end());
      auto* t = BaseConstructor::GetTripletFactor(tripletIndices[0], tripletIndices[1], tripletIndices[2]);
      std::array<INDEX,3> ti = tripletIndices;
      std::array<std::array<INDEX,2>,3> te{{{ti[0],ti[1]},{ti[0],ti[2]},{ti[1],ti[2]}}};
      std::array<std::array<INDEX,2>,3> tspE{{{std::min(n1,centerNode), std::max(n1,centerNode)}, {std::min(n2,centerNode),std::max(n2,centerNode)}, {n1,n2}}};
      Permutation<3> p(te,tspE);
      p.Invert();
      std::cout << "triplet      edges: " << "(" << te[0][0] << "," << te[0][1] << "), " << "(" << te[1][0] << "," << te[1][1] << "), " << "(" << te[2][0] << "," << te[2][1] << ")" << "\n";
      std::cout << "tripletSpoke edges: " << "(" << tspE[0][0] << "," << tspE[0][1] << "), " << "(" << tspE[1][0] << "," << tspE[1][1] << "), " << "(" << tspE[2][0] << "," << tspE[2][1] << ")" << "\n";
      std::cout << "permutation for covering message: " << INDEX(p[0]) << INDEX(p[1]) << INDEX(p[2]) << "\n";
      std::cout << "attention: not using permutation currently\n";
      auto* m = new TripletPlusSpokeCoverMessageContainer(MulticutTripletPlusSpokeCoverMessage(Permutation<3>({0,1,2})), t, tps, MulticutTripletPlusSpokeCoverMessage::size());
      //auto* m = new TripletPlusSpokeCoverMessageContainer(MulticutTripletPlusSpokeCoverMessage(Permutation<3>({n1,n2,centerNode},{tripletIndices[0],tripletIndices[1],tripletIndices[2]})), 
      //                                                    t, tps, MulticutTripletPlusSpokeCoverMessage::size());
      BaseConstructor::pd_.GetLP()->AddMessage(m);
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

      auto* m = new TripletPlusSpokeMessageContainer(MulticutTripletPlusSpokeMessage(sharedTripletEdgeTriplet,spokeEdgeTriplet,sharedTripletEdgeTripletPlusSpoke),t,tps,MulticutTripletPlusSpokeMessage::size());
      BaseConstructor::pd_.GetLP()->AddMessage(m);
      return m;
   }

   // rename to EnforceOddWheel. No single odd wheel factor is added
   INDEX AddOddWheelFactor(const INDEX centerNode, std::vector<INDEX> cycle)
   {
      CycleNormalForm(cycle);
      std::cout << "Enforce odd wheel with center node " << centerNode << " and cycle nodes ";
      for(auto i : cycle) std::cout << i << ",";
      std::cout << "\n";
   
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

      /*
      assert(!HasOddWheelFactor(centerNode, cycle));
      auto* o = new CriticalTripletCountingFactorContainer(MulticutCriticalTripletCountingFactor(cycle.size()), std::vector<REAL>(cycle.size(),0.0));
      //auto* o = new OddWheelFactorContainer(MulticutOddWheelFactor(cycle.size()), std::vector<REAL>(2*cycle.size(),0.0));
      BaseConstructor::pd_.GetLP()->AddFactor(o);
      oddWheelFactors_.insert(std::make_pair( std::make_tuple(centerNode,cycle), o ));
        
      // Link with triplet edges 
      for(INDEX c=0; c<cycle.size(); ++c) {
         // first create an intermediate critical triplet factor
         //auto* critical = new CriticalTripletFactorContainer(MulticutCriticalTripletFactor(), std::vector<REAL>(2,0.0));
         //BaseConstructor::pd_.GetLP()->AddFactor(critical); 
         // link with triplet counting factor

         std::array<INDEX,3> nodeIndices{centerNode,cycle[c],cycle[(c+1)%cycle.size()]};
         std::sort(nodeIndices.begin(), nodeIndices.end());
         typename BaseConstructor::TripletFactorContainer* t = BaseConstructor::GetTripletFactor(nodeIndices[0],nodeIndices[1],nodeIndices[2]);

         // To which edges in odd wheel do the triplet edges correspond.
         std::array<INDEX,2> cycleEdge{cycle[c],cycle[(c+1)%cycle.size()]};
         std::sort(cycleEdge.begin(),cycleEdge.end());
         std::array<INDEX,2> centerEdge1{centerNode,cycle[c]};
         std::sort(centerEdge1.begin(),centerEdge1.end());
         std::array<INDEX,2> centerEdge2{centerNode,cycle[(c+1)%cycle.size()]};
         std::sort(centerEdge2.begin(),centerEdge2.end());

         // the edges in the triplet factor are ordered as (i1,i2), (i1,i3), (i2,i3), i1<i2<i3, i.e. lexicographically
         std::array<std::array<INDEX,2>,3> tripletEdges{cycleEdge,centerEdge1,centerEdge2};
         std::sort(tripletEdges.begin(),tripletEdges.end());

         // the corresponding edges in the odd wheel factor are ordered as follows
         std::array<std::array<INDEX,2>,3> oddWheelEdges;
         std::array<INDEX,3> oddWheelEdgeIndices;
         oddWheelEdges[0] = cycleEdge;   oddWheelEdgeIndices[0] = c;
         oddWheelEdges[1] = centerEdge1; oddWheelEdgeIndices[1] = cycle.size() + c;
         oddWheelEdges[2] = centerEdge2; oddWheelEdgeIndices[2] = cycle.size() + ((c+1)%cycle.size());


         // try all six possible permutations explicitly
         std::cout << "triplet nodes = " << nodeIndices[0] << "," << nodeIndices[1] << "," << nodeIndices[2] << "\n";
                if(tripletEdges[0] == oddWheelEdges[0] && tripletEdges[1] == oddWheelEdges[1] &&  tripletEdges[2] == oddWheelEdges[2]) {
            //LinkTripletOddWheel(t,o, {oddWheelEdgeIndices[0], oddWheelEdgeIndices[1], oddWheelEdgeIndices[2]});
            LinkCountingTriplet(t,o,c,0);
            //std::cout << "connect triplet with nodes " << 
         } else if(tripletEdges[0] == oddWheelEdges[0] && tripletEdges[1] == oddWheelEdges[2] &&  tripletEdges[2] == oddWheelEdges[1]) {
            //LinkTripletOddWheel(t,o, {oddWheelEdgeIndices[0], oddWheelEdgeIndices[2], oddWheelEdgeIndices[1]});
            LinkCountingTriplet(t,o,c,0);
         } else if(tripletEdges[0] == oddWheelEdges[1] && tripletEdges[1] == oddWheelEdges[0] &&  tripletEdges[2] == oddWheelEdges[2]) {
            //LinkTripletOddWheel(t,o, {oddWheelEdgeIndices[1], oddWheelEdgeIndices[0], oddWheelEdgeIndices[2]});
            LinkCountingTriplet(t,o,c,1);
         } else if(tripletEdges[0] == oddWheelEdges[1] && tripletEdges[1] == oddWheelEdges[2] &&  tripletEdges[2] == oddWheelEdges[0]) {
            //LinkTripletOddWheel(t,o, {oddWheelEdgeIndices[1], oddWheelEdgeIndices[2], oddWheelEdgeIndices[0]});
            LinkCountingTriplet(t,o,c,1);
         } else if(tripletEdges[0] == oddWheelEdges[2] && tripletEdges[1] == oddWheelEdges[0] &&  tripletEdges[2] == oddWheelEdges[1]) {
            //LinkTripletOddWheel(t,o, {oddWheelEdgeIndices[2], oddWheelEdgeIndices[0], oddWheelEdgeIndices[1]});
            LinkCountingTriplet(t,o,c,2);
         } else if(tripletEdges[0] == oddWheelEdges[2] && tripletEdges[1] == oddWheelEdges[1] &&  tripletEdges[2] == oddWheelEdges[0]) {
            //LinkTripletOddWheel(t,o, {oddWheelEdgeIndices[2], oddWheelEdgeIndices[1], oddWheelEdgeIndices[0]});
            LinkCountingTriplet(t,o,c,2);
         } else {
            assert(false);
         }
      }
      return o;
      */



      // Link with unary edges
      // first link cycle edges
      /*
      for(INDEX c=0; c<cycle.size(); ++c) {
         std::array<INDEX,2> e({cycle[c], cycle[(c+1)%cycle.size()]});
         std::sort(e.begin(),e.end());
         auto* u = BaseConstructor::GetUnaryFactor(e[0],e[1]);
         using UnaryOddWheelMessageContainer = typename meta::at_c<typename FMC::MessageList, 3>::MessageContainerType; // for now only
         using UnaryOddWheelMessageType = typename UnaryOddWheelMessageContainer::MessageType;
         auto m = new UnaryOddWheelMessageContainer(UnaryOddWheelMessageType(c), u, o, 1);
         BaseConstructor::pd_.GetLP()->AddMessage(m);
      }
      // now the edges linking to center node
      for(INDEX c=0; c<cycle.size(); ++c) {
         std::array<INDEX,2> e({cycle[c], centerNode});
         std::sort(e.begin(),e.end());
         auto* u = BaseConstructor::GetUnaryFactor(e[0],e[1]);
         using UnaryOddWheelMessageContainer = typename meta::at_c<typename FMC::MessageList, 3>::MessageContainerType; // for now only
         using UnaryOddWheelMessageType = typename UnaryOddWheelMessageContainer::MessageType;
         auto m = new UnaryOddWheelMessageContainer(UnaryOddWheelMessageType(cycle.size()+c), u, o, 1);
         BaseConstructor::pd_.GetLP()->AddMessage(m);
      }
      */

      // Link with triplet edges 
      /*
      for(INDEX c=0; c<cycle.size(); ++c) {
         std::array<INDEX,3> nodeIndices{centerNode,cycle[c],cycle[(c+1)%cycle.size()]};
         std::sort(nodeIndices.begin(), nodeIndices.end());
         typename BaseConstructor::TripletFactorContainer* t = BaseConstructor::GetTripletFactor(nodeIndices[0],nodeIndices[1],nodeIndices[2]);

         // To which edges in odd wheel do the triplet edges correspond.
         std::array<INDEX,2> cycleEdge{cycle[c],cycle[(c+1)%cycle.size()]};
         std::sort(cycleEdge.begin(),cycleEdge.end());
         std::array<INDEX,2> centerEdge1{centerNode,cycle[c]};
         std::sort(centerEdge1.begin(),centerEdge1.end());
         std::array<INDEX,2> centerEdge2{centerNode,cycle[(c+1)%cycle.size()]};
         std::sort(centerEdge2.begin(),centerEdge2.end());

         // the edges in the triplet factor are ordered as (i1,i2), (i1,i3), (i2,i3), i1<i2<i3, i.e. lexicographically
         std::array<std::array<INDEX,2>,3> tripletEdges{cycleEdge,centerEdge1,centerEdge2};
         std::sort(tripletEdges.begin(),tripletEdges.end());

         // the corresponding edges in the odd wheel factor are ordered as follows
         std::array<std::array<INDEX,2>,3> oddWheelEdges;
         std::array<INDEX,3> oddWheelEdgeIndices;
         oddWheelEdges[0] = cycleEdge;   oddWheelEdgeIndices[0] = c;
         oddWheelEdges[1] = centerEdge1; oddWheelEdgeIndices[1] = cycle.size() + c;
         oddWheelEdges[2] = centerEdge2; oddWheelEdgeIndices[2] = cycle.size() + ((c+1)%cycle.size());


         // try all six possible permutations explicitly
         std::cout << "triplet nodes = " << nodeIndices[0] << "," << nodeIndices[1] << "," << nodeIndices[2] << "\n";
                if(tripletEdges[0] == oddWheelEdges[0] && tripletEdges[1] == oddWheelEdges[1] &&  tripletEdges[2] == oddWheelEdges[2]) {
            LinkTripletOddWheel(t,o, {oddWheelEdgeIndices[0], oddWheelEdgeIndices[1], oddWheelEdgeIndices[2]});
            //std::cout << "connect triplet with nodes " << 
         } else if(tripletEdges[0] == oddWheelEdges[0] && tripletEdges[1] == oddWheelEdges[2] &&  tripletEdges[2] == oddWheelEdges[1]) {
            LinkTripletOddWheel(t,o, {oddWheelEdgeIndices[0], oddWheelEdgeIndices[2], oddWheelEdgeIndices[1]});
         } else if(tripletEdges[0] == oddWheelEdges[1] && tripletEdges[1] == oddWheelEdges[0] &&  tripletEdges[2] == oddWheelEdges[2]) {
            LinkTripletOddWheel(t,o, {oddWheelEdgeIndices[1], oddWheelEdgeIndices[0], oddWheelEdgeIndices[2]});
         } else if(tripletEdges[0] == oddWheelEdges[1] && tripletEdges[1] == oddWheelEdges[2] &&  tripletEdges[2] == oddWheelEdges[0]) {
            LinkTripletOddWheel(t,o, {oddWheelEdgeIndices[1], oddWheelEdgeIndices[2], oddWheelEdgeIndices[0]});
         } else if(tripletEdges[0] == oddWheelEdges[2] && tripletEdges[1] == oddWheelEdges[0] &&  tripletEdges[2] == oddWheelEdges[1]) {
            LinkTripletOddWheel(t,o, {oddWheelEdgeIndices[2], oddWheelEdgeIndices[0], oddWheelEdgeIndices[1]});
         } else if(tripletEdges[0] == oddWheelEdges[2] && tripletEdges[1] == oddWheelEdges[1] &&  tripletEdges[2] == oddWheelEdges[0]) {
            LinkTripletOddWheel(t,o, {oddWheelEdgeIndices[2], oddWheelEdgeIndices[1], oddWheelEdgeIndices[0]});
         } else {
            assert(false);
         }
      }
      return o;
      */
   }

   // do zrobienia: include minDualIncrease
   INDEX FindOddWheels(const REAL minDualIncrease, const INDEX maxCuttingPlanesToAdd)
   {
      INDEX oddWheelsAdded = 0;
      //std::cout << "find odd wheel, " << BaseConstructor::tripletFactors_.size() << "\n";
      // do zrobienia: use reparametrization of edge potentials or of triplet potentials or of both simultaneously?
      // currently we assume we have triplet edges, which may not hold true. Better use original edges
      for(INDEX i=0; i<BaseConstructor::noNodes_; ++i) {
         // get all triplet factors attached to node i and build bipartite subgraph with doubled edges as described in Nowozin's thesis
         // make a hash of this. Give hints to load etc.
         std::map<INDEX,INDEX> origToCompressedNode; // compresses node indices
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
               //compressedToOrigNode.insert(std::make_pair(origToCompressedNode.size()-1,i1));
            }
            if(origToCompressedNode.find(i2) == origToCompressedNode.end()) {
               origToCompressedNode.insert(std::make_pair(i2, origToCompressedNode.size()));
               compressedToOrigNode.push_back(i2);
               //compressedToOrigNode.insert(std::make_pair(origToCompressedNode.size()-1,i2));
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
            // check if 110 or 101 are among the minimal labelings. where the first edge is the outer cycle edge
            // the corresponding labeling numbers are
            INDEX l1, l2;
            if(i < i1 && i < i2) { // the cycle edge is the last one
               l1 = 0; l2 = 1;
            } else if (i > i1 && i > i2) { // the cycle edge is the first one
               l1 = 1; l2 = 2;
            } else { // i1 < i < i2, the cycle edge is the second one
               l1 = 0; l2 = 2;
            }
            if( std::min((*t)[l1], (*t)[l2]) <= t->GetFactor()->LowerBound(*t) + eps ) {
               const INDEX i1c = origToCompressedNode[i1];
               const INDEX i2c = origToCompressedNode[i2];
               assert(i1c != i2c);
               g.AddEdge(i1c,noCompressedNodes + i2c,0.0); // do zrobienia: possibly only add edges for labelings which have a margin >= minDualIncrease
               g.AddEdge(noCompressedNodes + i1c,i2c,0.0); // do zrobienia: possibly only add edges for labelings which have a margin >= minDualIncrease
               uf.merge(i1c,noCompressedNodes + i2c);
               uf.merge(noCompressedNodes + i1c,i2c);
            }
         }
         // now check whether path exists between any given edges on graph
         //std::cout << "built auxiliary graph for node " << i << " with " << origToCompressedNode.size() << " nodes\n";
         for(INDEX j=0; j<noCompressedNodes; ++j) { // not nice: this has to be original number of nodes and bipartiteNumberOfNodes
            // find path from node j to node noNodes+j in g
            if(uf.connected(j,noCompressedNodes+j)) {
               //std::cout << "find path from " << j << " to " << j+noNodes << " and add corresponding wheel\n";
               auto path = mp.FindPath(j,noCompressedNodes+j,g);
               auto pathNormalized = std::get<1>(path);
               std::cout << "found compressed path ";
               for(INDEX k=0; k<pathNormalized.size(); ++k) {
                  std::cout << pathNormalized[k] << ", ";
               }
               std::cout << "\n";
               for(INDEX k=0; k<pathNormalized.size()-1; ++k) {
                  //assert(compressedToOrigNode.find(pathNormalized[k]%noCompressedNodes) != compressedToOrigNode.end());
                  pathNormalized[k] = compressedToOrigNode[pathNormalized[k]%noCompressedNodes];
               }
               //for(auto& k : pathNormalized) {
               //   assert(compressedToOrigNode.find(k) != compressedToOrigNode.end());
               //   k = compressedToOrigNode[k%noNodes];
               //}
               pathNormalized.resize(pathNormalized.size()-1); // first and last node coincide
               if(HasUniqueValues(pathNormalized)) { // possibly already add the subpath that is unique and do not search for it later. Indicate this with a std::vector<bool>
                  //assert(HasUniqueValues(pathNormalized)); // if not, a shorter subpath has been found. This subpath will be detected or has been deteced and has been added
                  CycleNormalForm(pathNormalized);
                  //CycleNormalForm called unnecesarily in AddOddWheelFactor
                  oddWheelsAdded += AddOddWheelFactor(i,pathNormalized);
               }
            }
         }

      }
      
      return oddWheelsAdded;
   }

   INDEX Tighten(const REAL minDualIncrease, const INDEX maxCuttingPlanesToAdd)
   {
      //std::cout << "remove this\n\n\n";
      //AddOddWheelFactor(1, std::vector<INDEX>{0,2,3});
      //return 0;
      const INDEX tripletsAdded = BaseConstructor::Tighten(minDualIncrease, maxCuttingPlanesToAdd);
      if(tripletsAdded > 0 ) {
         std::cout << "Added " << tripletsAdded << " triplet(s)\n";
         return tripletsAdded;
      } else {
         const INDEX oddWheelsAdded = FindOddWheels(minDualIncrease, maxCuttingPlanesToAdd - tripletsAdded);
         std::cout << "Added " << oddWheelsAdded << " factor(s) for odd wheel constraints\n";
         return tripletsAdded + oddWheelsAdded;
      }
      return 0;
   }


private:
   //std::map<std::tuple<INDEX,std::vector<INDEX>>, OddWheelFactorContainer*> oddWheelFactors_; // first entry denotes center node, then come cycle nodes with smallest node first, as computed by CycleNormalForm
   //std::map<std::tuple<INDEX,std::vector<INDEX>>, CriticalTripletCountingFactorContainer*> oddWheelFactors_; // first entry denotes center node, then come cycle nodes with smallest node first, as computed by CycleNormalForm
   // do zrobienia: possibly also hold pointer to factor
   std::vector<std::vector<std::tuple<INDEX,INDEX,typename BaseConstructor::TripletFactorContainer*>>> tripletByIndices_; // of triplet factor with indices (i1,i2,i3) exists, then (i1,i2,i3) will be in the vector of index i1, i2 and i3
   // the format for TripletPlusSpoke is (node1,node2, centerNode, spokeNode) and we assume n1<n2
   std::map<std::array<INDEX,4>,TripletPlusSpokeFactorContainer*> tripletPlusSpokeFactors_;
   //std::vector<std::vector<std::tuple<INDEX,typename BaseConstructor::UnaryFactorContainer*>>> unaryIndices_;

};
} // end namespace LP_MP

#endif // LP_MP_MULTICUT_CONSTRUCTOR_HXX

