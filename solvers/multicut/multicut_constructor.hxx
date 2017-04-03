#ifndef LP_MP_MULTICUT_CONSTRUCTOR_HXX
#define LP_MP_MULTICUT_CONSTRUCTOR_HXX

#include "LP_MP.h"
#include "multicut_unary_factor.hxx"
#include "multicut_triplet_factor.hxx"
#include "multicut_odd_wheel.hxx"
#include "lifted_multicut_factors_messages.hxx"

#include "union_find.hxx"
#include "graph.hxx"
#include "max_flow.hxx"

#include <unordered_map>
#include <unordered_set>
#include <queue>
#include <list>

#ifdef LP_MP_PARALLEL
#include <omp.h>
#endif

#include "andres/graph/graph.hxx"
#include "andres/graph/grid-graph.hxx"
#include "andres/graph/multicut/kernighan-lin.hxx"
#include "andres/graph/multicut/greedy-additive.hxx"
#include "andres/graph/multicut-lifted/kernighan-lin.hxx"
#include "andres/graph/multicut-lifted/greedy-additive.hxx"


namespace LP_MP {


// hash function for maps used in constructors. Do zrobienia: define hash functions used somewhere globally in config.hxx

enum class cut_type { multicut, maxcut };

template<class FACTOR_MESSAGE_CONNECTION, INDEX UNARY_FACTOR_NO, INDEX TRIPLET_FACTOR_NO,
   INDEX UNARY_TRIPLET_MESSAGE_0_NO, INDEX UNARY_TRIPLET_MESSAGE_1_NO, INDEX UNARY_TRIPLET_MESSAGE_2_NO,
   INDEX CONSTANT_FACTOR_NO, cut_type CUT_TYPE = cut_type::multicut>
class MulticutConstructor {
public:
   using MulticutConstructorType = MulticutConstructor<FACTOR_MESSAGE_CONNECTION, UNARY_FACTOR_NO, TRIPLET_FACTOR_NO, UNARY_TRIPLET_MESSAGE_0_NO, UNARY_TRIPLET_MESSAGE_1_NO, UNARY_TRIPLET_MESSAGE_2_NO, CONSTANT_FACTOR_NO, CUT_TYPE>;
   using FMC = FACTOR_MESSAGE_CONNECTION;

   using UnaryFactorContainer = meta::at_c<typename FMC::FactorList, UNARY_FACTOR_NO>;
   using TripletFactorContainer = meta::at_c<typename FMC::FactorList, TRIPLET_FACTOR_NO>;
   using ConstantFactorContainer = meta::at_c<typename FMC::FactorList, CONSTANT_FACTOR_NO>;
   using edge_triplet_message_0_container = typename meta::at_c<typename FMC::MessageList, UNARY_TRIPLET_MESSAGE_0_NO>::MessageContainerType;
   using edge_triplet_message_1_container = typename meta::at_c<typename FMC::MessageList, UNARY_TRIPLET_MESSAGE_1_NO>::MessageContainerType;
   using edge_triplet_message_2_container = typename meta::at_c<typename FMC::MessageList, UNARY_TRIPLET_MESSAGE_2_NO>::MessageContainerType;

   template<typename SOLVER>
   MulticutConstructor(SOLVER& pd)
   : lp_(&pd.GetLP()),
   //unaryFactors_(100,hash::array2),
   tripletFactors_(100,hash::array3)
   {
      //unaryFactors_.max_load_factor(0.7);
      tripletFactors_.max_load_factor(0.7);
      constant_factor_ = new ConstantFactorContainer(0.0);
      pd.GetLP().AddFactor(constant_factor_);
   }
   MulticutConstructor(const MulticutConstructorType& o)
      : unaryFactors_(o.unaryFactors_),
      unaryFactorsVector_(o.unaryFactorsVector_),
      tripletFactors_(o.tripletFactors_),
      noNodes_(o.noNodes_),
      constant_factor_(o.constant_factor_), 
      lp_(o.lp_) 
   {}

   ~MulticutConstructor()
   {
      //static_assert(std::is_same<typename UnaryFactorContainer::FactorType, MulticutUnaryFactor>::value,"");
      //static_assert(std::is_same<typename TripletFactorContainer::FactorType, MulticutTripletFactor>::value,"");
      //static_assert(std::is_same<typename MessageContainer::MessageType, MulticutUnaryTripletMessage<MessageSending::SRMP>>::value,"");
   }

   void End() 
   {
      // wait for the primal rounding to finish.
      std::cout << "wait for primal computation to end\n";
      if(primal_handle_.valid()) {
         auto labeling = primal_handle_.get();
         write_labeling_into_factors(labeling);
      }
   }

   void AddToConstant(const REAL delta) { constant_factor_->GetFactor()->AddToOffset(delta); }

   virtual UnaryFactorContainer* AddUnaryFactor(const INDEX i1, const INDEX i2, const REAL cost) // declared virtual so that derived class notices when unary factor is added
   {
      assert(i1 < i2);
      assert(!HasUnaryFactor(i1,i2));
      
      auto* u = new UnaryFactorContainer();
      (*u->GetFactor())[0] = cost;
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
   template<typename MESSAGE_CONTAINER>
   MESSAGE_CONTAINER* LinkUnaryTriplet(UnaryFactorContainer* u, TripletFactorContainer* t) 
   {
      auto* m = new MESSAGE_CONTAINER(u, t);
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
      if(!HasUnaryFactor(i1,i2)) {
         AddUnaryFactor(i1,i2,0.0);
      }
      if(!HasUnaryFactor(i1,i3)) {
         AddUnaryFactor(i1,i3,0.0);
      }
      if(!HasUnaryFactor(i2,i3)) {
         AddUnaryFactor(i2,i3,0.0);
      }
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
      LinkUnaryTriplet<edge_triplet_message_0_container>(GetUnaryFactor(tripletEdges[0][0],tripletEdges[0][1]), t);
      LinkUnaryTriplet<edge_triplet_message_1_container>(GetUnaryFactor(tripletEdges[1][0],tripletEdges[1][1]), t);
      LinkUnaryTriplet<edge_triplet_message_2_container>(GetUnaryFactor(tripletEdges[2][0],tripletEdges[2][1]), t);
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
      if(tripletsAdded < 0.6*maxCuttingPlanesToAdd) {
         std::cout << "Additionally search via shortest paths for violated constraints\n";
         if(CUT_TYPE == cut_type::multicut) {
            tripletsAdded += find_violated_cycles_multicut(maxCuttingPlanesToAdd - tripletsAdded);
         } else if(CUT_TYPE == cut_type::maxcut) {
            tripletsAdded += find_violated_cycles_maxcut(maxCuttingPlanesToAdd - tripletsAdded);
         }
         std::cout << "Added " << tripletsAdded << " triplet(s) out of " <<  maxCuttingPlanesToAdd << " in total\n";
      }
      return tripletsAdded;
   }

   template<
   class InputIt1, class InputIt2, class OutputIt, class Compare, class Merge >
   static OutputIt set_intersection_merge
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
         const REAL cost_ij = (*it.second->GetFactor())[0];
         assert(i<j);
         adjacency_list[i][adjacency_list_count[i]] = std::make_tuple(j,cost_ij);
         adjacency_list_count[i]++;
         adjacency_list[j][adjacency_list_count[j]] = std::make_tuple(i,cost_ij);
         adjacency_list_count[j]++;
      }

      // Sort the adjacency list, for fast intersections later
      auto adj_sort = [](const auto a, const auto b) { return std::get<0>(a) < std::get<0>(b); };

#pragma omp parallel for schedule(guided)
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
      std::vector<std::tuple<INDEX,INDEX,INDEX,REAL>> triplet_candidates;
#pragma omp parallel
      {
         std::vector<intersection_type> commonNodes(noNodes_);
         std::vector<std::tuple<INDEX,INDEX,INDEX,REAL>> triplet_candidates_per_thread;
#pragma omp for schedule(guided)
         for(INDEX c=0; c<unaryFactorsVector_.size(); ++c) {
            const REAL cost_ij = (*unaryFactorsVector_[c].second->GetFactor())[0];
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

               if(CUT_TYPE == cut_type::multicut) {
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
                     triplet_candidates_per_thread.push_back(std::make_tuple(i,j,k,guaranteed_dual_increase));
                  } 
               } else if(CUT_TYPE == cut_type::maxcut) {
                  const REAL lb = std::min(0.0, cost_ij) + std::min(0.0, cost_ik) + std::min(0.0, cost_jk);
                  const REAL best_labeling = std::min({0.0, cost_ij+cost_ik, cost_ij+cost_jk, cost_ij+cost_jk});
                  assert(lb <= best_labeling+eps);
                  const REAL guaranteed_dual_increase = best_labeling - lb;
                  if(guaranteed_dual_increase > 0.0) {
                     triplet_candidates_per_thread.push_back(std::make_tuple(i,j,k,guaranteed_dual_increase));
                  } 
               }
            }
         }
#pragma omp critical
         {
            triplet_candidates.insert(triplet_candidates.end(), triplet_candidates_per_thread.begin(), triplet_candidates_per_thread.end()); 
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

   INDEX find_violated_cycles_maxcut(const INDEX max_triplets_to_add)
   {
      assert(false);
      // build doubled graph for separation. We search for a path with an odd number of negative ( = cut) edges. 
      // Hence negative values are between two components, while positive ones are inside each component.
      std::vector<std::tuple<INDEX,INDEX,REAL>> edges;
      edges.reserve(unaryFactorsVector_.size());
      std::vector<INDEX> number_outgoing_arcs(2*noNodes_,0); // number of arcs outgoing arcs of each node
      INDEX number_arcs_total = 0;
      for(auto& it : unaryFactorsVector_) {
         const REAL v = (*it.second->GetFactor())[0];
         const INDEX i = std::get<0>(it.first);
         const INDEX j = std::get<1>(it.first);
         if(v != 0.0) {
            number_outgoing_arcs[2*i]++;
            number_outgoing_arcs[2*i+1]++;
            number_outgoing_arcs[2*j]++;
            number_outgoing_arcs[2*j+1]++;
            edges.push_back(std::make_tuple(i,j,v));
         }
         number_arcs_total += 4;
      }

      Graph g(noNodes_, number_arcs_total, number_outgoing_arcs); // graph consisting of positive edges
      for(const auto& it : edges) {
         const INDEX i = std::get<0>(it);
         const INDEX j = std::get<1>(it);
         const REAL v = std::get<2>(it);
         assert(i<j);
         if(v > 0.0) {
            g.add_edge(2*i,2*j,v);
            g.add_edge(2*i+1,2*j+1,v);
         } else if(v < 0.0) {
            g.add_edge(2*i,2*j+1,-v);
            g.add_edge(2*i+1,2*j,-v);
         }
      }
      g.sort();

      std::sort(edges.begin(), edges.end(), [](auto& a, auto& b) { return std::abs(std::get<2>(a)) > std::abs(std::get<2>(b)); }); 

      UnionFind uf(2*noNodes_);
      auto merge_edges = [&uf](const INDEX i, const INDEX j, const REAL v) {
         assert(v != 0.0);
         if(v > 0.0) {
            uf.merge(2*i,2*j);
            uf.merge(2*i+1,2*j+1);
         } else if(v < 0.0) {
            uf.merge(2*i,2*j+1);
            uf.merge(2*i+1,2*j); 
         }
      };

      INDEX tripletsAdded = 0;

      INDEX e=0;
      REAL initial_th;
      for(; e<edges.size(); ++e) {
         const INDEX i = std::get<0>(edges[e]);
         const INDEX j = std::get<1>(edges[e]);
         const REAL v = std::get<2>(edges[e]);
         merge_edges(i,j,v);
         if(uf.connected(2*i,2*i+1) || uf.connected(2*j, 2*j+1)) {
            initial_th = std::abs(v);
            break;
         } 
      }

      bool zero_th_iteration = true;
      std::vector<bool> already_searched(noNodes_, false);
      for(REAL th=0.5*initial_th; th>=eps || zero_th_iteration; th*=0.1) {
         if(th < eps) {
            if(tripletsAdded <= 0.01*max_triplets_to_add) {
               std::cout << "additional separation with no guaranteed dual increase, i.e. th = 0\n";
               th = 0.0;
               zero_th_iteration = false;
            } else {
               break;
            }
         }
         // first update union find datastructure
         for(;e<edges.size(); ++e) {
            const INDEX i = std::get<0>(edges[e]);
            const INDEX j = std::get<1>(edges[e]);
            const REAL v = std::get<2>(edges[e]);
            if(std::abs(v) >= th) {
               merge_edges(i,j,v); 
            }
         }
         using CycleType = std::tuple<REAL, std::vector<INDEX>>;
         std::vector<CycleType > cycles;
#pragma omp parallel 
         {
            std::vector<CycleType > cycles_local;
            BfsData mp2(g);
#pragma for schedule(guided)
            for(INDEX i=0; i<noNodes_; ++i) {
               if(!already_searched[i] && uf.thread_safe_connected(2*i, 2*i+1)) {
                  already_searched[i] = true;
                  auto cycle = mp2.FindPath(2*i, 2*i+1, g, th);
                  const REAL dualIncrease = std::get<0>(cycle);
                  assert(std::get<1>(cycle).size() > 0);
                  if(std::get<1>(cycle).size() > 0) {
                     cycles_local.push_back( std::make_tuple(dualIncrease, std::move(std::get<1>(cycle))) );
                  } else {
                     throw std::runtime_error("No path found although there should be one"); 
                  } 
               } 
            }
#pragma omp critical
            {
               cycles.insert(cycles.end(), cycles_local.begin(), cycles_local.end());
            }
         }
         // sort by guaranteed increase in decreasing order
         std::sort(cycles.begin(), cycles.end(), [](const CycleType& i, const CycleType& j) { return std::get<0>(i) > std::get<0>(j); });
         for(auto& cycle : cycles) {
            if(std::get<1>(cycle).size() > 2) {
               tripletsAdded += AddCycle(std::move(std::get<1>(cycle)));
               if(tripletsAdded > max_triplets_to_add) {
                  return tripletsAdded;
               }
            }
         }
      }

      return tripletsAdded;
   }

   INDEX find_violated_cycles_multicut(const INDEX maxTripletsToAdd)
   {
      std::vector<std::tuple<INDEX,INDEX,REAL,bool> > negative_edges; // endpoints, edge cost, searched positive path with given endpoints?
      // we can speed up compution by skipping path searches for node pairs which lie in different connected components. Connectedness is stored in a union find structure
      REAL pos_th = 0.0;
      std::vector<INDEX> number_outgoing_arcs(noNodes_,0); // number of arcs outgoing arcs of each node
      INDEX number_outgoing_arcs_total = 0;
      for(auto& it : unaryFactorsVector_) {
         const REAL v = (*it.second->GetFactor())[0];
         const INDEX i = std::get<0>(it.first);
         const INDEX j = std::get<1>(it.first);
         if(v >= 0.0) {
            pos_th = std::max(pos_th, v);
            number_outgoing_arcs[i]++;
            number_outgoing_arcs[j]++;
            number_outgoing_arcs_total += 2;
         } else {
            negative_edges.push_back(std::make_tuple(i,j,v,false));
         }
      }
      if(negative_edges.size() == 0 || negative_edges.size() == unaryFactorsVector_.size()) { return 0; }

      Graph posEdgesGraph(noNodes_, number_outgoing_arcs_total, number_outgoing_arcs); // graph consisting of positive edges
      for(auto& it : unaryFactorsVector_) {
         const REAL v = (*it.second->GetFactor())[0];
         const INDEX i = std::get<0>(it.first);
         const INDEX j = std::get<1>(it.first);
         assert(i<j);
         if(v >= 0.0) {
            posEdgesGraph.add_arc(i,j,v);
            posEdgesGraph.add_arc(j,i,v);
         }
      }
      posEdgesGraph.sort();

      // do zrobienia: possibly add reparametrization of triplet factors additionally

      std::sort(negative_edges.begin(), negative_edges.end(), [](const auto& e1, const auto& e2)->bool {
            return std::get<2>(e1) < std::get<2>(e2);
            });


      // now search for every negative edge for most negative path from end point to starting point. Do zrobienia: do this in parallel
      // here, longest path is sought after only the edges with positive taken into account
      // the cost of the path is the minimum of the costs of its edges. The guaranteed increase in the dual objective is v_min > -v ? -v : v_min

      //MostViolatedPathData mp(posEdgesGraph);
      //BfsData mp(posEdgesGraph);

      UnionFind uf(noNodes_);
      INDEX tripletsAdded = 0;
      const REAL initial_th = 0.6*std::min(-std::get<2>(negative_edges[0]), pos_th);
      bool zero_th_iteration = true;
      for(REAL th=initial_th; th>=eps || zero_th_iteration; th*=0.1) {
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
            const REAL v = (*it.second->GetFactor())[0];
            const INDEX i = std::get<0>(it.first);
            const INDEX j = std::get<1>(it.first);
            if(v >= th) {
               uf.merge(i,j);   
            }
         }
         using CycleType = std::tuple<REAL, std::vector<INDEX>>;
         std::vector<CycleType > cycles;
#pragma omp parallel 
         {
            std::vector<CycleType > cycles_local;
            BfsData mp2(posEdgesGraph);
#pragma for schedule(guided)
            for(INDEX c=0; c<negative_edges.size(); ++c) {
               const INDEX i = std::get<0>(negative_edges[c]);
               const INDEX j = std::get<1>(negative_edges[c]);
               const REAL v = std::get<2>(negative_edges[c]);
               const bool already_used_for_path_search = std::get<3>(negative_edges[c]);
               //if(-v <= th) break;
               //if(already_used_for_path_search) continue;
               if(-v > th && !already_used_for_path_search && uf.thread_safe_connected(i,j)) {
                  //auto cycle = mp.FindPath(i,j,posEdgesGraph);
                  auto cycle = mp2.FindPath(i,j,posEdgesGraph, th);
                  const REAL dualIncrease = std::min(-v, std::get<0>(cycle));
                  assert(std::get<1>(cycle).size() > 0);
                  if(std::get<1>(cycle).size() > 0) {
                     cycles_local.push_back( std::make_tuple(dualIncrease, std::move(std::get<1>(cycle))) );
                     //tripletsAdded += AddCycle(std::get<1>(cycle));
                     //if(tripletsAdded > maxTripletsToAdd) {
                     //   return tripletsAdded;
                     //}
                  } else {
                     throw std::runtime_error("No path found although there should be one"); 
                  }
               }
            }
#pragma omp critical
            {
               cycles.insert(cycles.end(), cycles_local.begin(), cycles_local.end());
            }
         }
         // sort by guaranteed increase in decreasing order
         std::sort(cycles.begin(), cycles.end(), [](const CycleType& i, const CycleType& j) { return std::get<0>(i) > std::get<0>(j); });
         for(auto& cycle : cycles) {
            if(std::get<1>(cycle).size() > 2) {
               tripletsAdded += AddCycle(std::move(std::get<1>(cycle)));
               if(tripletsAdded > maxTripletsToAdd) {
                  return tripletsAdded;
               }
            }
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



   bool CheckPrimalConsistency() const
   {
      if(CUT_TYPE == cut_type::multicut) {
         std::cout << "checking primal feasibility for multicut ... ";
         UnionFind uf(noNodes_);
         for(const auto& e : unaryFactorsVector_) {
            UnaryFactorContainer* f = e.second; 
            if(f->GetFactor()->primal()[0] == false) {
               // connect components 
               const INDEX i = std::get<0>(e.first);
               const INDEX j = std::get<1>(e.first);
               uf.merge(i,j);
            }
         }
         for(const auto& e : unaryFactorsVector_) {
            UnaryFactorContainer* f = e.second; 
            if(f->GetFactor()->primal()[0] == true) {
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
      } else if(CUT_TYPE == cut_type::maxcut) {
         return false;
      } else {
         throw std::runtime_error("");
      }
   }

   void round()
   {
      std::cout << "compute multicut primal with GAEC + KLj\n";
      andres::graph::Graph<> graph(noNodes_);
      std::vector<REAL> edgeValues;
      edgeValues.reserve(unaryFactorsVector_.size());

      for(const auto& e : unaryFactorsVector_) {
         graph.insertEdge(e.first[0], e.first[1]);
         edgeValues.push_back((*e.second->GetFactor())[0]);
      }

      primal_handle_ = std::async(std::launch::async, gaec_klj, std::move(graph), std::move(edgeValues));
   }

   static std::vector<char> gaec_klj(andres::graph::Graph<> g, std::vector<REAL> edge_values)
   {
      std::vector<char> labeling(g.numberOfEdges(), 0);
      andres::graph::multicut::greedyAdditiveEdgeContraction(g, edge_values, labeling);
      andres::graph::multicut::kernighanLin(g, edge_values, labeling, labeling);
      return labeling;
   }

   void write_labeling_into_factors(const std::vector<char>& labeling) {
      assert(labeling.size() <= unaryFactorsVector_.size());
      for(INDEX c=0; c<labeling.size(); ++c) {
         auto* f = unaryFactorsVector_[c].second;
         f->GetFactor()->primal()[0] = labeling[c];
         f->ComputePrimalThroughMessages();
      }

      // possibly, additional edges have been added because of tightening. infer labeling of those from union find datastructure
      if(labeling.size() < unaryFactorsVector_.size()) {
         UnionFind uf(noNodes_);
         for(INDEX c=0; c<labeling.size(); ++c) {
            UnaryFactorContainer* f = unaryFactorsVector_[c].second; 
            if(f->GetFactor()->primal()[0] == false) {
               // connect components 
               const INDEX i = std::get<0>(unaryFactorsVector_[c].first);
               const INDEX j = std::get<1>(unaryFactorsVector_[c].first);
               uf.merge(i,j);
            }
         }
         std::cout << "built union find structure, propagate information now\n";
         for(INDEX c=labeling.size(); c<unaryFactorsVector_.size(); ++c) {
            UnaryFactorContainer* f = unaryFactorsVector_[c].second; 
            const INDEX i = std::get<0>(unaryFactorsVector_[c].first);
            const INDEX j = std::get<1>(unaryFactorsVector_[c].first);
            if(uf.connected(i,j)) {
               f->GetFactor()->primal()[0] = false;
            } else {
               f->GetFactor()->primal()[0] = true;
            }
            f->ComputePrimalThroughMessages();
         }
      }
   }

   // use GAEC and Kernighan&Lin algorithm of andres graph package to compute primal solution
   void ComputePrimal()
   {
      if(CUT_TYPE == cut_type::multicut) {
         if(!primal_handle_.valid()) { 
            round();
            return;
         }

         const auto primal_state = primal_handle_.wait_for(std::chrono::seconds(0));

         if(primal_state == std::future_status::deferred) {
            assert(false); // this should not happen, we launch immediately!
            throw std::runtime_error("asynchronuous primal multicut rounding was deferred, but this should not happen");
         } else if(primal_state == std::future_status::ready) {

            std::cout << "collect multicut rounding result\n";
            auto labeling = primal_handle_.get();
            write_labeling_into_factors(labeling);

            std::cout << "restart primal rounding\n";
            round();

         } else {
            std::cout << "multicut rounding is currently running.\n";
         }
      } else if(CUT_TYPE == cut_type::maxcut) {
         assert(false);
      } else {
         throw std::runtime_error("");
      } 
   }

protected:

   decltype(std::async(std::launch::async, gaec_klj, andres::graph::Graph<>(0), std::vector<REAL>{})) primal_handle_;

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



template<
   class MULTICUT_CONSTRUCTOR,
   INDEX ODD_3_WHEEL_FACTOR_NO,
   INDEX TRIPLET_ODD_3_WHEEL_MESSAGE_012_NO, INDEX TRIPLET_ODD_3_WHEEL_MESSAGE_013_NO, INDEX TRIPLET_ODD_3_WHEEL_MESSAGE_023_NO, INDEX TRIPLET_ODD_3_WHEEL_MESSAGE_123_NO
>
class MulticutOddWheelConstructor : public MULTICUT_CONSTRUCTOR {
public:
   using FMC = typename MULTICUT_CONSTRUCTOR::FMC;
   using BaseConstructor = MULTICUT_CONSTRUCTOR;

   using odd_3_wheel_factor_container = meta::at_c<typename FMC::FactorList, ODD_3_WHEEL_FACTOR_NO>;
   using triplet_odd_3_wheel_message_012_container = meta::at_c<typename FMC::MessageList, TRIPLET_ODD_3_WHEEL_MESSAGE_012_NO>;
   using triplet_odd_3_wheel_message_013_container = meta::at_c<typename FMC::MessageList, TRIPLET_ODD_3_WHEEL_MESSAGE_013_NO>;
   using triplet_odd_3_wheel_message_023_container = meta::at_c<typename FMC::MessageList, TRIPLET_ODD_3_WHEEL_MESSAGE_023_NO>;
   using triplet_odd_3_wheel_message_123_container = meta::at_c<typename FMC::MessageList, TRIPLET_ODD_3_WHEEL_MESSAGE_123_NO>;

   //using TripletPlusSpokeFactorContainer = meta::at_c<typename FMC::FactorList, TRIPLET_PLUS_SPOKE_FACTOR_NO>;
   //using TripletPlusSpokeMessageContainer = typename meta::at_c<typename FMC::MessageList, TRIPLET_PLUS_SPOKE_MESSAGE_NO>::MessageContainerType;
   //using TripletPlusSpokeCoverMessageContainer = typename meta::at_c<typename FMC::MessageList, TRIPLET_PLUS_SPOKE_COVER_MESSAGE_NO>::MessageContainerType;

   template<typename SOLVER>
   MulticutOddWheelConstructor(SOLVER& pd) :
      BaseConstructor(pd),
      odd_3_wheel_factors_(100, hash::array4)
      //, tripletPlusSpokeFactors_(100,hash::array4) { 
   {
      //tripletPlusSpokeFactors_.max_load_factor(0.7); 
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

   odd_3_wheel_factor_container* add_odd_3_wheel_factor(const INDEX i0, const INDEX i1, const INDEX i2, const INDEX i3)
   {
      assert(!has_odd_3_wheel_factor(i0,i1,i2,i3));

      auto* f = new odd_3_wheel_factor_container();
      this->lp_->AddFactor(f);
      odd_3_wheel_factors_.insert(std::make_pair(std::array<INDEX,4>({i0,i1,i2,i3}), f));

      if(!this->HasTripletFactor(i0,i1,i2)) {
         this->AddTripletFactor(i0,i1,i2);
      }
      auto* t012 = this->GetTripletFactor(i0,i1,i2);
      auto* m012 = new triplet_odd_3_wheel_message_012_container(t012, f);
      this->lp_->AddMessage(m012);

      if(!this->HasTripletFactor(i0,i1,i3)) {
         this->AddTripletFactor(i0,i1,i3);
      }
      auto* t013 = this->GetTripletFactor(i0,i1,i3);
      auto* m013 = new triplet_odd_3_wheel_message_013_container(t013, f);
      this->lp_->AddMessage(m013);

      if(!this->HasTripletFactor(i0,i2,i3)) {
         this->AddTripletFactor(i0,i2,i3);
      }
      auto* t023 = this->GetTripletFactor(i0,i2,i3);
      auto* m023 = new triplet_odd_3_wheel_message_023_container(t023, f);
      this->lp_->AddMessage(m023);

      if(!this->HasTripletFactor(i1,i2,i3)) {
         this->AddTripletFactor(i1,i2,i3);
      }
      auto* t123 = this->GetTripletFactor(i1,i2,i3);
      auto* m123 = new triplet_odd_3_wheel_message_123_container(t123, f);
      this->lp_->AddMessage(m123);

      return f;
   }
   bool has_odd_3_wheel_factor(const INDEX i0, const INDEX i1, const INDEX i2, const INDEX i3) const
   {
      assert(i0 < i1 && i1 < i2 && i2 < i3);
      return odd_3_wheel_factors_.find(std::array<INDEX,4>({i0,i1,i2,i3})) != odd_3_wheel_factors_.end();
   }
   odd_3_wheel_factor_container* get_odd_3_wheel_factor(const INDEX i0, const INDEX i1, const INDEX i2, const INDEX i3) const
   {
      return odd_3_wheel_factors_.find(std::array<INDEX,4>({i0,i1,i2,i3}))->second; 
   }
   /*
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
   */
   /*
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
   */

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

         std::array<INDEX,4> nodes({n1,n2,centerNode, cycle[0]});
         std::sort(nodes.begin(), nodes.end());
         if(!has_odd_3_wheel_factor(nodes[0], nodes[1], nodes[2], nodes[3])) {
            add_odd_3_wheel_factor(nodes[0], nodes[1], nodes[2], nodes[3]);
            ++tripletPlusSpokesAdded;
         }

         /*
         if(n1>n2) {
            std::swap(n1,n2); 
         }
         if(!HasTripletPlusSpokeFactor(n1,n2,centerNode, cycle[0])) {
            AddTripletPlusSpokeFactor(n1,n2,centerNode,cycle[0]);
            LinkTripletPlusSpokeFactor(n1,n2,centerNode,cycle[0],n1);
            LinkTripletPlusSpokeFactor(n1,n2,centerNode,cycle[0],n2);
            ++tripletPlusSpokesAdded;
         }
         */
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
      const REAL c01 = (*this->GetUnaryFactor(triplet[0], triplet[1])->GetFactor())[0];
      const REAL c02 = (*this->GetUnaryFactor(triplet[0], triplet[2])->GetFactor())[0];
      const REAL c12 = (*this->GetUnaryFactor(triplet[1], triplet[2])->GetFactor())[0];
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
      Graph g(noBipartiteCompressedNodes,4*compressedEdges.size(), no_outgoing_arcs);
      BfsData mp(g);
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

   }

   INDEX Tighten(const INDEX maxCuttingPlanesToAdd)
   {
      const INDEX tripletsAdded = BaseConstructor::Tighten(maxCuttingPlanesToAdd);
      if(tripletsAdded > 0.1*maxCuttingPlanesToAdd) {
         return tripletsAdded;
      } else {
         const INDEX oddWheelsAdded = FindOddWheels(maxCuttingPlanesToAdd);
         std::cout << "Added " << oddWheelsAdded << " factors for odd wheel constraints\n";
         return oddWheelsAdded;
      }
      assert(false);
   }

   

private:
   std::atomic<int> kwas_delete;

   std::unordered_map<std::array<INDEX,4>, odd_3_wheel_factor_container*,decltype(hash::array4)> odd_3_wheel_factors_;

   std::vector<std::vector<std::tuple<INDEX,INDEX,typename BaseConstructor::TripletFactorContainer*>>> tripletByIndices_; // of triplet factor with indices (i1,i2,i3) exists, then (i1,i2,i3) will be in the vector of index i1, i2 and i3
   // the format for TripletPlusSpoke is (node1,node2, centerNode, spokeNode) and we assume n1<n2
   // hash for std::array<INDEX,4>
   //std::unordered_map<std::array<INDEX,4>,TripletPlusSpokeFactorContainer*,decltype(hash::array4)> tripletPlusSpokeFactors_;
};

template<
   class MULTICUT_CONSTRUCTOR,
   INDEX LIFTED_MULTICUT_CUT_FACTOR_NO,
   INDEX CUT_EDGE_LIFTED_MULTICUT_FACTOR_NO, INDEX LIFTED_EDGE_LIFTED_MULTICUT_FACTOR_NO
   >
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
   bool CheckPrimalConsistency() const
   {
      const bool multicutConsistent = MULTICUT_CONSTRUCTOR::CheckPrimalConsistency();
      if(!multicutConsistent) {
         return false;
      }
      
      //collect connectivity information with union find w.r.t. base edges
      UnionFind uf(MULTICUT_CONSTRUCTOR::noNodes_);
      for(const auto& e : baseEdges_) {
         if(e.f->GetFactor()->primal()[0] == false) {
            uf.merge(e.i,e.j);
         }
      }
      for(const auto& e : liftedEdges_) {
         if(e.f->GetFactor()->primal()[0] == false) {
            if(!uf.connected(e.i,e.j)) {
              return false;
           }
        }
      }
      return true;
   }

   void ComputePrimal()
   {}

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
      REAL weight() const { return (*f->GetFactor())[0]; }
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

