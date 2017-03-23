/*
 *  Clean and efficient reimplementation of separation for k-ary cycle inequalities, as proposed by David Sontag, Do Kook Choe and Yitao Li in UAI 2012 paper.
 */

#ifndef LP_MP_CYCLE_INEQUALITIES_HXX
#define LP_MP_CYCLE_INEQUALITIES_HXX

#include <stdlib.h>
#include <vector>
#include <algorithm>
#include <set>
#include <list>
#include <map>
#include <queue>

#include "config.hxx"
#include "vector.hxx"
#include "graph.hxx"
#include "union_find.hxx"

namespace LP_MP {

template<typename MRF_CONSTRUCTOR>
class triplet_search 
{
public:
   typedef std::vector<std::vector<std::pair<INDEX, REAL> > > adj_type;

   triplet_search(const MRF_CONSTRUCTOR& mrf)
      : gm_(mrf)
   {}
   ~triplet_search() {};

   std::vector<triplet_candidate> search()
   {
   // Initialize adjacency list (filled in later) TODO: only do this when needed
   std::vector<std::vector<INDEX> > adjacency_list(gm_.GetNumberOfVariables());

   // Construct adjacency list for the graph
   // Iterate over all of the edges (we do this by looking at the edge intersection sets)
   for(size_t factorId=0; factorId<gm_.GetNumberOfPairwiseFactors(); factorId++) {
      auto vars = gm_.GetPairwiseVariables(factorId);
      // Get the two nodes i & j
      const size_t i=std::get<0>(vars);
      const size_t j=std::get<1>(vars);
      assert(i<j);
      adjacency_list[i].push_back(j);
      adjacency_list[j].push_back(i);
   }

   // Sort the adjacency list, for fast intersections later
#pragma omp parallel for
   for(INDEX i=0; i < adjacency_list.size(); i++) {
      std::sort(adjacency_list[i].begin(), adjacency_list[i].end());
   }

  std::vector<triplet_candidate> triplet_candidates;
  INDEX index=0;
			
  // Iterate over all of the edge intersection sets
#pragma omp parallel 
  {
     std::vector<INDEX> commonNodes(gm_.GetNumberOfVariables()-1);
     std::vector<triplet_candidate> triplet_candidates_local;
#pragma omp for
     for(size_t factorId=0; factorId<gm_.GetNumberOfPairwiseFactors(); factorId++) {
        auto vars = gm_.GetPairwiseVariables(factorId);
        const INDEX i=std::get<0>(vars);
        const INDEX j=std::get<1>(vars);
        const auto& factor_ij = *gm_.GetPairwiseFactor(i,j)->GetFactor();
        const REAL lb_ij = factor_ij.LowerBound();

        // Now find all neighbors of both i and j to see where the triangles are
        auto intersects_iter_end = set_intersection(adjacency_list[i].begin(), adjacency_list[i].end(), adjacency_list[j].begin(), adjacency_list[j].end(), commonNodes.begin());

        for(std::vector<INDEX>::const_iterator n=commonNodes.begin(); n != intersects_iter_end; ++n) {
           INDEX k = *n;

           // Since a triplet shows up multiple times, we only consider it for the case when i<j<k 
           if(!(j<k))
              continue;

           const auto& factor_ik = *gm_.GetPairwiseFactor(i,k)->GetFactor();
           const auto& factor_jk = *gm_.GetPairwiseFactor(j,k)->GetFactor();
           const REAL boundIndep = lb_ij + factor_ik.LowerBound() + factor_jk.LowerBound();
           const REAL boundCycle = minimizeTriangle(factor_ij, factor_ik, factor_jk);

           const REAL bound = boundCycle - boundIndep; 
           assert(bound >=  - eps);
           if(bound > eps) {
              triplet_candidate t(i,j,k, bound);
              triplet_candidates_local.push_back(t);
           }
        }
     }
#pragma omp critical
     triplet_candidates.insert(triplet_candidates.end(), triplet_candidates_local.begin(), triplet_candidates_local.end());
  }

  std::sort(triplet_candidates.begin(), triplet_candidates.end());

  return std::move(triplet_candidates);
}

protected:
   template<typename PAIRWISE_REPAM>
   REAL minimizeTriangle(const PAIRWISE_REPAM& factor_ij, const PAIRWISE_REPAM& factor_ik, const PAIRWISE_REPAM& factor_jk) const
   {
      REAL max_val = std::numeric_limits<REAL>::infinity();

      // Fix value of the first variable
      const INDEX dim1 = factor_ij.dim1();
      const INDEX dim2 = factor_ij.dim2();
      const INDEX dim3 = factor_ik.dim2();
      for(INDEX i1=0; i1<dim1; ++i1) {
         for(INDEX i2=0; i2<dim2; ++i2) {
            for(INDEX i3=0; i3<dim3; ++i3) {
               max_val = std::min(max_val, factor_ij(i1,i2) + factor_ik(i1,i3) + factor_jk(i2,i3));
            }
         }
      }
      return max_val;
   }

   const MRF_CONSTRUCTOR& gm_;
};

// search for violated k-ary cycles either with k-projection graph or with expanded projection graph (indicated by EXTENDED template flag)
template<typename MRF_CONSTRUCTOR, bool EXTENDED=false>
class k_ary_cycle_inequalities_search
{
public:
   k_ary_cycle_inequalities_search(const MRF_CONSTRUCTOR& mrf) 
      : gm_(mrf)
   {}

   std::vector<triplet_candidate> search(const INDEX max_triplets = std::numeric_limits<INDEX>::max())
   {
      if(!EXTENDED) {
         construct_projection_graph();
      } else {
         auto partitions = compute_partitions(); 
         construct_projection_graph(std::move(partitions));
      }
      return find_cycles(max_triplets);
   }

   template<typename PAIRWISE_REPAM>
      static matrix<REAL> row_minima(const PAIRWISE_REPAM& f);
   template<typename PAIRWISE_REPAM>
      static matrix<REAL> column_minima(const PAIRWISE_REPAM& f);
   template<typename PAIRWISE_REPAM>
      static matrix<REAL> principal_minima(const PAIRWISE_REPAM& f, const matrix<REAL>& _column_minima);

protected:

   template<typename PAIRWISE_REPAM> 
   std::pair<std::vector<bool>,std::vector<bool>> compute_partitions(const PAIRWISE_REPAM& f);
   std::vector<std::vector<std::vector<bool>>> compute_partitions(); 

   void construct_projection_graph(const std::vector<std::vector<std::vector<bool>>> partitions = {});

   template<typename PAIRWISE_REPAM, typename V>
   REAL compute_projection_weight_singleton_2(const PAIRWISE_REPAM& f, const std::vector<bool>& part_i, const INDEX x2, const V& row_minima);
   template<typename PAIRWISE_REPAM, typename V>
   REAL compute_projection_weight_singleton_1(const PAIRWISE_REPAM& f, const INDEX x1, const std::vector<bool>& part_j, const V& column_minima);
   template<typename PAIRWISE_REPAM>
   REAL compute_projection_weight_on_partitions(const PAIRWISE_REPAM& f, const std::vector<bool>& part_i, const std::vector<bool>& part_j);

   std::vector<triplet_candidate> find_cycles(const INDEX max_triplets);
   void triangulate(std::vector<triplet_candidate>& triplet_candidates, std::tuple<REAL,std::vector<INDEX>>& path);

   std::vector<std::tuple<INDEX,INDEX,REAL>> projection_edges_;
   std::vector<INDEX> proj_graph_offsets_;
   std::vector<INDEX> proj_graph_to_gm_node_;
   std::vector<std::tuple<INDEX,INDEX,REAL>> proj_graph_edges_;
   Graph proj_graph_;

   const MRF_CONSTRUCTOR& gm_;
};


// possibly do the three claculations below in one pass for greater efficiency
template<typename MRF_CONSTRUCTOR, bool EXTENDED>
template<typename PAIRWISE_REPAM>
matrix<REAL> k_ary_cycle_inequalities_search<MRF_CONSTRUCTOR, EXTENDED>::row_minima(const PAIRWISE_REPAM& f)
{
   matrix<REAL> _row_minima(f.dim1(),2);
   for(INDEX x1=0; x1<f.dim1(); ++x1) {
      // find smallest and second smallest entry w.r.t. second index given current first index
      REAL smallest = std::numeric_limits<REAL>::infinity();
      REAL second_smallest = std::numeric_limits<REAL>::infinity();
      for(INDEX x2=0; x2<f.dim2(); ++x2) {
         const REAL val = f(x1,x2);
         const REAL min = std::min(smallest, val);
         const REAL max = std::max(smallest, val);
         smallest = min;
         second_smallest = std::min(max, second_smallest);
      }
      _row_minima(x1,0) = smallest;
      _row_minima(x1,1) = second_smallest;
   }
   return std::move(_row_minima);
}

template<typename MRF_CONSTRUCTOR, bool EXTENDED>
template<typename PAIRWISE_REPAM>
matrix<REAL> k_ary_cycle_inequalities_search<MRF_CONSTRUCTOR, EXTENDED>::column_minima(const PAIRWISE_REPAM& f)
{
   matrix<REAL> _column_minima(f.dim2(),2);
   std::fill(_column_minima.begin(), _column_minima.end(), std::numeric_limits<REAL>::infinity());
   for(INDEX x1=0; x1<f.dim1(); ++x1) {
      for(INDEX x2=0; x2<f.dim2(); ++x2) {
         const REAL val = f(x1,x2);
         const REAL min = std::min(_column_minima(x2,0), val);
         const REAL max = std::max(_column_minima(x2,0), val);
         _column_minima(x2,0) = min;
         _column_minima(x2,1) = std::min(max, _column_minima(x2,1));
      }
   }
   return std::move(_column_minima);
}

template<typename MRF_CONSTRUCTOR, bool EXTENDED>
template<typename PAIRWISE_REPAM>
matrix<REAL> k_ary_cycle_inequalities_search<MRF_CONSTRUCTOR, EXTENDED>::principal_minima(const PAIRWISE_REPAM& f, const matrix<REAL>& _column_minima)
{
   // possibly can be computed more efficiently. Note that at most 4 values need to be stored
   matrix<REAL> _principal_minima(f.dim1(), f.dim2());
   for(INDEX x1=0; x1<f.dim1(); ++x1) {
      REAL smallest = std::numeric_limits<REAL>::infinity();
      REAL second_smallest = std::numeric_limits<REAL>::infinity();
      INDEX smallest_ind = 0;
      for(INDEX x2=0; x2<f.dim2(); ++x2) {
         const REAL val_xij = f(x1,x2);
         const REAL x = _column_minima(x2,0) == val_xij ? _column_minima(x2,1) : _column_minima(x2,0);
         if(x<smallest) {
            smallest = x;
            smallest_ind = x2;
         }
      }
      for(INDEX x2=0; x2<f.dim2(); ++x2) {
         const REAL val_xij = f(x1,x2);
         const REAL x = _column_minima(x2,0) == val_xij ? _column_minima(x2,1) : _column_minima(x2,0);
         if(x < second_smallest && x2 != smallest_ind) {
            second_smallest = x;
         }
      }
      for(INDEX x2=0; x2<f.dim2(); ++x2) {
         _principal_minima(x1,x2) = smallest;
      }
      _principal_minima(x1,smallest_ind) = second_smallest; 
   }
   return std::move(_principal_minima);
}

// Given an undirected graph, finds odd-signed cycles.
template<typename MRF_CONSTRUCTOR, bool EXTENDED>
std::vector<triplet_candidate> 
k_ary_cycle_inequalities_search<MRF_CONSTRUCTOR, EXTENDED>::find_cycles(const INDEX max_triplets)
{
   REAL largest_th;
   UnionFind uf(proj_graph_.size());
   INDEX e=0;

   auto merge_edge = [&uf](const INDEX m, const INDEX n, const REAL s) {
      if(s < 0) {
        uf.merge(2*n,2*m);
        uf.merge(2*m,2*n);
        uf.merge(2*n+1,2*m+1);
        uf.merge(2*m+1,2*n+1); 
     } else {
        uf.merge(2*n,2*m+1);
        uf.merge(2*m+1,2*n);
        uf.merge(2*n+1,2*m);
        uf.merge(2*m,2*n+1); 
     }
   };

   for(; e<projection_edges_.size(); ++e) {
      const INDEX i = std::get<0>(projection_edges_[e]);
      const INDEX j = std::get<1>(projection_edges_[e]);
      const REAL s = std::get<2>(projection_edges_[e]);
      merge_edge(i,j,s);

      if(uf.connected(2*i,2*i+1) || uf.connected(2*j,2*j+1)) {
         largest_th = std::abs(s);
         break;
      } 
   }

   std::vector<triplet_candidate> triplet_candidates;
   std::vector<bool> already_searched(proj_graph_to_gm_node_.size(),false);
   //BfsData bfs(proj_graph_);
   // first update union find datastructure by merging additional edges with cost greater than th
   REAL th = 0.5*largest_th;
   for(INDEX iter=0; iter<8 && th>=eps; ++iter, th*=0.1) {
      // update connectivity information
      for(; e<projection_edges_.size(); ++e) {
         const INDEX i = std::get<0>(projection_edges_[e]);
         const INDEX j = std::get<1>(projection_edges_[e]);
         const REAL s = std::get<2>(projection_edges_[e]);
         if(std::abs(s) >= th) {
            merge_edge(i,j,s);
         } else {
            break;
         }
      }

      // now actually search for odd signed cycles
#pragma omp parallel
      {
         std::vector<triplet_candidate> triplet_candidates_local;
         BfsData bfs(proj_graph_);
#pragma omp for
         for(INDEX i=0; i<proj_graph_to_gm_node_.size(); ++i) {
            if(!already_searched[i] && uf.thread_safe_connected(2*i, 2*i+1)) {
               already_searched[i] = true;
               auto path = bfs.FindPath(2*i, 2*i+1, proj_graph_, th);
               assert(std::get<1>(path).size() >= 3);
               triangulate(triplet_candidates_local, path);
            }
         }
#pragma omp critical
         triplet_candidates.insert(triplet_candidates.end(), triplet_candidates_local.begin(), triplet_candidates_local.end()); 
      }
      if(triplet_candidates.size() > max_triplets) {
         break;
      }
   }

   std::sort(triplet_candidates.begin(), triplet_candidates.end());
   if(triplet_candidates.size() > 0) {
      assert(triplet_candidates[0].cost >= triplet_candidates.back().cost);
   }
   triplet_candidates.erase( unique( triplet_candidates.begin(), triplet_candidates.end() ), triplet_candidates.end() );
   return std::move(triplet_candidates);
}

 
template<typename MRF_CONSTRUCTOR, bool EXTENDED>
void
k_ary_cycle_inequalities_search<MRF_CONSTRUCTOR, EXTENDED>::triangulate(std::vector<triplet_candidate>& triplet_candidates, std::tuple<REAL,std::vector<INDEX>>& path)
{
   const REAL th = std::get<0>(path);
   auto cycle = std::get<1>(path);
   // halve indices, so we get node indices in projection graph, then project back onto indices of graphical model
   for(INDEX& i : cycle) {
      i = proj_graph_to_gm_node_[i/2];
   }
   assert(cycle[0] == cycle.back());
   cycle.resize(cycle.size()-1);

   // possibly cycle has a sub-cycle, or even uses an edge twice consecutively, still add it

   assert(cycle.size() >= 3);

   for(INDEX i=1; i<cycle.size()-1; ++i) {
      const INDEX j = cycle[i];
      const INDEX k = cycle[i+1];
      assert(j!=k);
      if(j != cycle[0] && k != cycle[0]) {
         triplet_candidates.push_back(triplet_candidate(cycle[0], j,k, th));
      }
   }
   /*
   auto min_node = std::min_element(cycle.begin(), cycle.end());
   for(auto it=cycle.begin(); it<min_node-2; ++it) {
      triplet_candidates.push_back(triplet_candidate(*it,*(it+1),*min_node, th));
   }
   for(auto it=min_node+1; it<cycle.end()-2; ++it) {
      triplet_candidates.push_back(triplet_candidate(*it,*(it+1),*min_node, th));
   }
   if(min_node != cycle.begin() && min_node != cycle.end()-1) { // when min_node is somewhere in the middle.
      triplet_candidates.push_back(triplet_candidate(*cycle.begin(), *cycle.rbegin(), *min_node, th));
   }
   */
}


template<typename MRF_CONSTRUCTOR, bool EXTENDED>
template<typename PAIRWISE_REPAM> 
std::pair<std::vector<bool>,std::vector<bool>> 
k_ary_cycle_inequalities_search<MRF_CONSTRUCTOR, EXTENDED>::compute_partitions(const PAIRWISE_REPAM& f)
{
   if(f.dim1() <= 3 && f.dim2() <= 3) {
      return std::make_pair(std::vector<bool>{}, std::vector<bool>{});
   }
   std::vector<std::tuple<INDEX,INDEX,REAL>> sorted_factor(f.dim1()*f.dim2());
   for(INDEX x1=0; x1<f.dim1(); ++x1) {
      for(INDEX x2=0; x2<f.dim2(); ++x2) {
         sorted_factor[x1*f.dim2() + x2] = std::make_tuple(x1,x2,f(x1,x2));
      }
   }
   std::sort(sorted_factor.begin(), sorted_factor.end(), [](auto a, auto b) { return std::get<2>(a) > std::get<2>(b); });
   UnionFind uf(f.dim1() + f.dim2());
   std::vector<bool> part_i(f.dim1());
   std::vector<bool> part_j(f.dim2());
   for(auto x : sorted_factor) {
      const INDEX x1 = std::get<0>(x);
      const INDEX x2 = std::get<1>(x);
      const REAL cost = std::get<2>(x);

      if(!uf.connected(x1, f.dim1() + x2)) {
         // check if merging c1 and c2 would result in one partition
         if(uf.count() > 2) {
            uf.merge(x1, f.dim1() + x2);
         } else {
            assert(uf.count() == 2);
            const INDEX c1 = uf.find(x1);
            const INDEX c2 = uf.find(f.dim1() + x2);
            // record labels of partition c1
            for(INDEX y1=0; y1<f.dim1(); ++y1) {
               if(uf.find(y1) == c1) {
                  part_i[y1] = true; 
               } else {
                  part_i[y1] = false; 
               }
            }
            for(INDEX y2=0; y2<f.dim2(); ++y2) {
               if(uf.find(f.dim1() + y2) == c1) {
                  part_j[y2] = true; 
               } else {
                  part_j[y2] = false; 
               }
            }
            break;
         }
      }
   }
   // normalize partitions
   if(part_i[0] == false) {
      part_i.flip();
   }
   if(part_j[0] == false) {
      part_j.flip();
   }
   // if partition consists of a singleton on any side, do not add it
   const INDEX sum_1 = std::count(part_i.begin(), part_i.end(), true);
   if(sum_1 == 1 || sum_1 == f.dim1()-1) {
      part_i.clear();
   }
   const INDEX sum_2 = std::count(part_j.begin(), part_j.end(), true);
   if(sum_2 == 1 || sum_2 == f.dim2()-1) {
      part_j.clear();
   }
   if(f.dim1() <= 3) {
      part_i.clear();
   }
   if(f.dim2() <= 3) {
      part_j.clear();
   }

   return std::move(std::make_pair(std::move(part_i), std::move(part_j)));
}

template<typename MRF_CONSTRUCTOR, bool EXTENDED>
std::vector<std::vector<std::vector<bool>>> 
k_ary_cycle_inequalities_search<MRF_CONSTRUCTOR, EXTENDED>::compute_partitions()
{
   std::vector<std::vector<std::vector<bool>>> partitions(gm_.GetNumberOfVariables());
#pragma omp parallel
   {
      auto partitions_local = partitions;
#pragma omp for nowait
      for(size_t factorId=0; factorId<gm_.GetNumberOfPairwiseFactors(); factorId++) {
         auto vars = gm_.GetPairwiseVariables(factorId);
         const auto& factor = *gm_.GetPairwiseFactor(factorId)->GetFactor();
         const size_t i=std::get<0>(vars);
         const size_t j=std::get<1>(vars);
         assert(i<j);
         const auto part_ij = compute_partitions(factor);
         const auto& part_i = std::get<0>(part_ij);
         if(part_i.size() > 0) {
            partitions_local[i].push_back(std::move(part_i));
         }
         const auto& part_j = std::get<1>(part_ij);
         if(part_j.size() > 0) {
            partitions_local[j].push_back(std::move(part_j));
         }
      }
#pragma omp critical
      for(INDEX i=0; i<partitions.size(); ++i) {
         partitions[i].insert(partitions[i].end(), partitions_local[i].begin(), partitions_local[i].end());
      }
   }

   // remove duplicate partitions
#pragma omp parallel for
   for(INDEX i=0; i<partitions.size(); ++i) {
      std::sort(partitions[i].begin(), partitions[i].end());
      partitions[i].erase( unique( partitions[i].begin(), partitions[i].end() ), partitions[i].end() );
   }

   return std::move(partitions);
}

template<typename MRF_CONSTRUCTOR, bool EXTENDED>
template<typename PAIRWISE_REPAM, typename V>
REAL 
k_ary_cycle_inequalities_search<MRF_CONSTRUCTOR, EXTENDED>::compute_projection_weight_singleton_2(const PAIRWISE_REPAM& f, const std::vector<bool>& part_i, const INDEX x2, const V& _row_minima)
{
   // same dimension if part_i[x1] and x2
   assert(f.dim1() == part_i.size());
   assert(f.dim1() == _row_minima.dim1() && 2 == _row_minima.dim2());
   REAL min_same_part = std::numeric_limits<REAL>::infinity();
   REAL min_different_part = std::numeric_limits<REAL>::infinity();
   for(INDEX x1=0; x1<f.dim1(); ++x1) {
      if(part_i[x1]) {
         min_same_part = std::min(min_same_part, f(x1,x2)); 
      } else {
         min_different_part = std::min(min_different_part, f(x1,x2));
      }
   }

   for(INDEX x1=0; x1<f.dim1(); ++x1) {
      const REAL val_xij = f(x1,x2);
      const REAL row_min = _row_minima(x1,0) == val_xij ? _row_minima(x1,1) : _row_minima(x1,0);
      if(part_i[x1]) {
         min_different_part = std::min(min_different_part, row_min);
      } else {
         min_same_part = std::min(min_same_part, row_min);
      }
   }
   return min_same_part - min_different_part;
}

template<typename MRF_CONSTRUCTOR, bool EXTENDED>
template<typename PAIRWISE_REPAM, typename V>
REAL 
k_ary_cycle_inequalities_search<MRF_CONSTRUCTOR, EXTENDED>::compute_projection_weight_singleton_1(const PAIRWISE_REPAM& f, const INDEX x1, const std::vector<bool>& part_j, const V& _column_minima)
{
   assert(f.dim2() == part_j.size());
   assert(f.dim2() == _column_minima.dim1() && 2 == _column_minima.dim2());
   REAL min_same_part = std::numeric_limits<REAL>::infinity();
   REAL min_different_part = std::numeric_limits<REAL>::infinity();
   for(INDEX x2=0; x2<f.dim2(); ++x2) {
      if(part_j[x2]) {
         min_same_part = std::min(min_same_part, f(x1,x2)); 
      } else {
         min_different_part = std::min(min_different_part, f(x1,x2));
      }
   }

   for(INDEX x2=0; x2<f.dim2(); ++x2) {
      const REAL val_xij = f(x1,x2);
      const REAL column_min = _column_minima(x2,0) == val_xij ? _column_minima(x2,1) : _column_minima(x2,0);
      if(part_j[x2]) {
         min_different_part = std::min(min_different_part, column_min);
      } else {
         min_same_part = std::min(min_same_part, column_min);
      }
   }
   return min_same_part - min_different_part;
}

template<typename MRF_CONSTRUCTOR, bool EXTENDED>
template<typename PAIRWISE_REPAM>
REAL 
k_ary_cycle_inequalities_search<MRF_CONSTRUCTOR, EXTENDED>::compute_projection_weight_on_partitions(const PAIRWISE_REPAM& f, const std::vector<bool>& part_i, const std::vector<bool>& part_j)
{
   assert(f.dim1() == part_i.size());
   assert(f.dim2() == part_j.size());
   REAL min_same_part = std::numeric_limits<REAL>::infinity();
   REAL min_different_part = std::numeric_limits<REAL>::infinity();

   for(INDEX x1=0; x1<f.dim1(); ++x1) {
      for(INDEX x2=0; x2<f.dim2(); ++x2) {
         const REAL val_xij = f(x1,x2);
         if(part_i[x1] == part_j[x2]) {
            min_same_part = std::min(min_same_part, val_xij);
         } else {
            min_different_part = std::min(min_different_part, val_xij); 
         }
      }
   }
   return min_same_part - min_different_part;
}

template<typename MRF_CONSTRUCTOR, bool EXTENDED>
void 
k_ary_cycle_inequalities_search<MRF_CONSTRUCTOR, EXTENDED>::construct_projection_graph(const std::vector<std::vector<std::vector<bool>>> partitions)
{
   proj_graph_offsets_ = std::vector<INDEX>(gm_.GetNumberOfVariables());
   proj_graph_offsets_[0] = 0;
   for(INDEX i=1; i<gm_.GetNumberOfVariables(); i++) {
      proj_graph_offsets_[i] = gm_.GetNumberOfLabels(i-1);
      if(EXTENDED) { 
         proj_graph_offsets_[i] += partitions[i-1].size();
      }
   }

   std::partial_sum(proj_graph_offsets_.begin(), proj_graph_offsets_.end(), proj_graph_offsets_.begin());
   INDEX proj_graph_nodes = proj_graph_offsets_.back() + gm_.GetNumberOfLabels(gm_.GetNumberOfVariables()-1);
   if(EXTENDED) {
      proj_graph_nodes += partitions.back().size();
   }

   proj_graph_to_gm_node_ = std::vector<INDEX>(proj_graph_nodes);
   {
      INDEX c=0;
      for(INDEX i=0; i<gm_.GetNumberOfVariables(); i++) {
         INDEX no_proj_graph_nodes_for_label = gm_.GetNumberOfLabels(i);
         if(EXTENDED) {
            no_proj_graph_nodes_for_label += partitions[i].size();
         }
         for(INDEX x=0; x<no_proj_graph_nodes_for_label; ++x) {
            proj_graph_to_gm_node_[c] = i;
            ++c;
         }
      }
      assert(c == proj_graph_nodes);
   }

   std::vector<INDEX> no_outgoing_arcs(2*proj_graph_nodes,0);
   projection_edges_.clear();

   auto add_to_projection_edges = [&no_outgoing_arcs](auto& projection_edges, const INDEX n, const INDEX m, const REAL val) {
      if(std::abs(val) >= eps && !std::isnan(val)) {          
         projection_edges.push_back(std::make_tuple(m,n,val));
         no_outgoing_arcs[2*m]++;
         no_outgoing_arcs[2*n]++;
         no_outgoing_arcs[2*m+1]++;
         no_outgoing_arcs[2*n+1]++;
      } 
   };

   for(size_t factorId=0; factorId<gm_.GetNumberOfPairwiseFactors(); factorId++) {
      // Get the two nodes i & j and the edge intersection set. Put in right order.
      const INDEX i = std::get<0>(gm_.GetPairwiseVariables(factorId));
      const INDEX j = std::get<1>(gm_.GetPairwiseVariables(factorId));

      // Check to see if i and j have at least two states each -- otherwise, cannot be part of any frustrated edge
      if(gm_.GetNumberOfLabels(i) <= 1 || gm_.GetNumberOfLabels(j) <= 1)
         continue;

      // For each of their singleton states efficiently compute edge weights
      const auto& factor_ij = *gm_.GetPairwiseFactor(i,j)->GetFactor(); // better retrieve by factor id

      const auto row_min = row_minima(factor_ij);
      const auto col_min = column_minima(factor_ij);
      const auto principal_min = principal_minima(factor_ij, col_min);

      assert(i<j);
      for(INDEX xi=0; xi<factor_ij.dim1(); xi++) {
         const INDEX m = proj_graph_offsets_[i] + xi;
         assert(i == proj_graph_to_gm_node_[m]);

         for(INDEX xj=0; xj<factor_ij.dim2(); xj++) {
            const INDEX n = proj_graph_offsets_[j] + xj;
            assert(j == proj_graph_to_gm_node_[n]);

            const REAL val_xij = factor_ij(xi,xj);

            const REAL val_not_xi = row_min(xi,0) == val_xij ? row_min(xi,1) : row_min(xi,0);
            const REAL val_not_xj = col_min(xj,0) == val_xij ? col_min(xj,1) : col_min(xj,0);

            // val_s < 0 means same projection < different projection, > 0 the opposite
            // Hence we search for a cycle with an odd number of entries > 0      
            const REAL cost_projection_same = std::min(val_xij, principal_min(xi,xj));
            const REAL cost_projection_different = std::min(val_not_xi, val_not_xj);
            const REAL val_s = cost_projection_same - cost_projection_different;

            // TODO: use threshold here, to make next stage faster
            add_to_projection_edges(projection_edges_,n,m,val_s);
         }
      }

      if(EXTENDED) {
         // add edge weights between each general projection and each singleton state 
         for(INDEX x1=0; x1<factor_ij.dim1(); ++x1) {
            const INDEX m = proj_graph_offsets_[i] + x1;
            for(INDEX p2=0; p2<partitions[j].size(); ++p2) {
               const INDEX n = proj_graph_offsets_[j] + gm_.GetNumberOfLabels(j) + p2;
               const REAL val_s = compute_projection_weight_singleton_1(factor_ij, x1, partitions[j][p2], col_min);
               add_to_projection_edges(projection_edges_,n,m,val_s);
            }
         }
         for(INDEX x2=0; x2<factor_ij.dim2(); ++x2) {
            const INDEX m = proj_graph_offsets_[j] + x2;
            for(INDEX p1=0; p1<partitions[i].size(); ++p1) {
               const INDEX n = proj_graph_offsets_[i] + gm_.GetNumberOfLabels(i) + p1;
               const REAL val_s = compute_projection_weight_singleton_2(factor_ij, partitions[i][p1], x2, row_min);
               add_to_projection_edges(projection_edges_,n,m,val_s);
            }
         }

         // compute edge weights between general projections
         for(INDEX p1=0; p1<partitions[i].size(); ++p1) {
            const INDEX n = proj_graph_offsets_[i] + gm_.GetNumberOfLabels(i) + p1;
            for(INDEX p2=0; p2<partitions[j].size(); ++p2) {
               const INDEX m = proj_graph_offsets_[j] + gm_.GetNumberOfLabels(j) + p2;
                  const REAL val_s = compute_projection_weight_on_partitions(factor_ij, partitions[i][p1], partitions[j][p2]);
                  add_to_projection_edges(projection_edges_,n,m,val_s);
               }
            }
         }
      }


      proj_graph_ = Graph(2*proj_graph_nodes, 4*projection_edges_.size(), no_outgoing_arcs);
      for(const auto edge : projection_edges_) {
         const INDEX m = std::get<0>(edge);
         const INDEX n = std::get<1>(edge);
         const REAL s = std::get<2>(edge);
         if(s < 0) {
            proj_graph_.add_arc(2*n,2*m,-s);
            proj_graph_.add_arc(2*m,2*n,-s);
            proj_graph_.add_arc(2*n+1,2*m+1,-s);
            proj_graph_.add_arc(2*m+1,2*n+1,-s); 
         } else {
            proj_graph_.add_arc(2*n,2*m+1,s);
            proj_graph_.add_arc(2*m+1,2*n,s);
            proj_graph_.add_arc(2*n+1,2*m,s);
            proj_graph_.add_arc(2*m,2*n+1,s); 
         }
      }

      proj_graph_.sort();

      std::sort(projection_edges_.begin(), projection_edges_.end(), [](auto a, auto b) { return std::get<2>(a) > std::get<2>(b); }); 
   }
} // end namespace LP_MP

#endif // LP_MP_CYCLE_INEQUALITIES_HXX

