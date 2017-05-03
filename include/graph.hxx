#ifndef LP_MP_GRAPH
#define LP_MP_GRAPH

// graph and related utilities for e.g. shortest path based separation algorithms
namespace LP_MP {

   // Structure for storing a candidate triplet cluster for tightening
   struct triplet_candidate {
      triplet_candidate(const INDEX i1, const INDEX i2, const INDEX i3, const REAL c)
         : cost(c)
      {
         i = std::min({i1,i2,i3});
         j = std::max(std::min(i1,i2), std::min(std::max(i1,i2),i3));
         k = std::max({i1,i2,i3}); 
         assert(i < j && j < k);
      }

      REAL cost;
      int i,j,k;

   };

   bool operator<(const triplet_candidate& l, const triplet_candidate& r) {
      return l.cost > r.cost;
   }
   bool operator==(const triplet_candidate& l, const triplet_candidate& r) {
      return l.i == r.i && l.j == r.j && l.k == r.k;
   }


   // efficient graph structure for shortest path computation.
   // arcs are held contiguously w.r.t. tail node and descending w.r.t. cost.
   // Hence, iterating over all outgoing arcs from a given node is fast. 
   // Also early stopping, when we go over a thresholded graph is possible, as arcs are ordered w.r.t. cost when tail node is fixed.
   class Graph {
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
      Graph() {}
      template<typename VEC>
      Graph(const INDEX no_nodes, const INDEX no_arcs, const VEC& number_outgoing_arcs)
      : nodes_(no_nodes,{nullptr,nullptr}),
         arcs_(no_arcs) 
      {
         assert(no_nodes > 0);
         if(no_arcs > 0) {
            nodes_[0].first = &arcs_[0];
            nodes_[0].last = nodes_[0].first;
            for(INDEX i=1; i<no_nodes; ++i) {
               nodes_[i].first = nodes_[i-1].first + number_outgoing_arcs[i-1];
               nodes_[i].last = nodes_[i].first;
            }
         }
      }
      Graph(const Graph& o) = delete;
      //{
      //   assert(false);
      //}
      Graph(Graph&& o) = default;
      Graph& operator=(Graph&& o) = default;

      INDEX size() const { return nodes_.size(); }
      void add_arc(INDEX i, INDEX j, REAL cost) {
         assert(cost >= 0.0);
         assert(i<nodes_.size());
         assert(j<nodes_.size());
         *nodes_[i].last = arc({cost, &nodes_[j]});
         if(i+1<nodes_.size()) {
            assert(nodes_[i].last < nodes_[i+1].first);
         }
         nodes_[i].last += 1;
      }
      void add_edge(const INDEX i, const INDEX j, const REAL cost)
      {
         add_arc(i,j,cost);
         add_arc(j,i,cost);
      }

      INDEX operator[](node* n) const { return n - &nodes_[0]; }
      const node& operator[](INDEX i) const { return nodes_[i]; }

      void sort()
      {
#pragma omp parallel for
         for(INDEX i=0; i<nodes_.size(); ++i) {
            std::sort(nodes_[i].begin(), nodes_[i].end(), [](auto& a, auto& b) { return a.cost > b.cost; });
         }

         for(INDEX i=0; i<nodes_.size()-1; ++i) {
            assert(nodes_[i].last == nodes_[i+1].first);
         }
      }

      private:
         std::vector<node> nodes_;
         std::vector<arc> arcs_;
   };

   struct BfsData {
      struct Item { INDEX parent; REAL cost; INDEX flag; };
      BfsData(const INDEX no_nodes)
      {
         d.resize(no_nodes);
         for(INDEX i=0; i<d.size(); ++i) {
            d[i].flag = 0;
         }
         flag1 = 0;
         flag2 = 1; 
      }
      BfsData(const Graph& g) : BfsData(g.size())
      {}
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

      //static auto no_mask_op = [](const INDEX i, const INDEX j, const REAL weight) { return true; };
      static bool no_mask_op(const INDEX, const INDEX, const REAL) { return true;}
      
      // do bfs with thresholded costs and iteratively lower threshold until enough cycles are found
      // only consider edges that have cost equal or larger than th
      std::tuple<REAL, std::vector<INDEX>> FindPath(const INDEX startNode, const INDEX endNode, const Graph& g, const REAL th = 0.0)
      {
         return FindPath(startNode, endNode, g, th, no_mask_op);
      }

      template<typename MASK_OP>
      std::tuple<REAL,std::vector<INDEX>> 
      FindPath(
         const INDEX startNode, const INDEX endNode, const Graph& g, const REAL th,
         MASK_OP mask_op
         )
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


            if(Labelled1(i)) {
               for(auto* a=g[i].begin(); a->cost>=th && a!=g[i].end(); ++a) { 
                  auto* head = a->head;
                  const INDEX j = g[head];

                  if(mask_op(i,j,a->cost)) {

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
               }
            } else {
               assert(Labelled2(i));
               for(auto* a=g[i].begin(); a->cost>=th && a!=g[i].end(); ++a) { 
                  auto* head = a->head;
                  const INDEX j = g[head];

                  if(mask_op(i,j,a->cost)) {

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
         return std::make_tuple(-std::numeric_limits<REAL>::infinity(),std::vector<INDEX>(0));
      }

protected:
   std::vector<Item> d;
   std::deque<std::array<INDEX,2>> visit; // node number, distance from start or end 
   INDEX flag1, flag2;
};

// possible idea for faster search: order arcs so that those with cost 0 are inserted at the beginning, those with cost 1 at the end of each nodes arc list
class multicut_path_search : public BfsData {
public:
   using BfsData::BfsData;

   bool path_exists(const INDEX startNode, const INDEX endNode, Graph& g) // find path of length 1, where we assume that all edges have cost 0 or 1. We do not assume anymore that arcs are ordered with descending cost
   {
      Reset();
      visit.push_back({startNode, 0});
      Label1(startNode);
      Parent(startNode) = startNode;
      Cost(startNode) = 0;
      visit.push_back({endNode, 0});
      Label2(endNode);
      Parent(endNode) = endNode;
      Cost(endNode) = 0;

      while(!visit.empty()) {
         const INDEX i = visit.front()[0];
         const INDEX distance = visit.front()[1];
         visit.pop_front();

         if(Labelled1(i)) {
            for(auto* a=g[i].begin(); a!=g[i].end(); ++a) { 
               auto* head = a->head;
               const INDEX j = g[head];
               const REAL cost_to_j = a->cost + Cost(i);

               assert(cost_to_j == 0.0 || cost_to_j == 1.0 || cost_to_j == 2.0);
               if(cost_to_j <= 1.0) {

                  if(!Labelled(j)) {
                     visit.push_back({j, distance+1});
                     Parent(j) = i;
                     Cost(j) = cost_to_j;
                     Label1(j);
                  } else if(Labelled2(j)) { // shortest path found
                     if(cost_to_j + Cost(j) == 1.0) { return true; }
                  } else if(Labelled1(j)) { // new shortest path found
                     Cost(j) = std::min(Cost(j), cost_to_j); // should never occur: if there exists path of length 1, then a path of length 0 should not be possible. Larger values are not considered in any case.
                  }
               }
            }
         } else {
            assert(Labelled2(i));
            for(auto* a=g[i].begin(); a!=g[i].end(); ++a) { 
               auto* head = a->head;
               const INDEX j = g[head];
               const REAL cost_to_j = a->cost + Cost(i);

               assert(cost_to_j == 0.0 || cost_to_j == 1.0 || cost_to_j == 2.0);
               if(cost_to_j <= 1.0) {

                  if(!Labelled(j)) {
                     visit.push_back({j, distance+1});
                     Parent(j) = i;
                     Cost(j) = cost_to_j;
                     Label2(j);
                  } else if(Labelled1(j)) { // shortest path found
                     if(cost_to_j + Cost(j) == 1.0) { return true; }
                  } else if(Labelled2(j)) {
                     Cost(j) = std::min(Cost(j), cost_to_j); 
                  }

               }
            }
         }
      }

      return false;
   } 
};

};

#endif
