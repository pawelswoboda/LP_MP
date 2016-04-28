#ifndef LP_MP_MIN_COST_FLOW_FACTOR_HXX
#define LP_MP_MIN_COST_FLOW_FACTOR_HXX

#include "MinCost/MinCost.h"

namespace LP_MP {

// build factor consisting of min cost network flow graph
// use as implementation the MinCost solver from Kolmogorov
// the reparametrization storage holds the min cost flow solver

// do zrobienia: a pointer to the min cost flow solver is held in the factor and in the ReparametrizationStorage, which is redundant. Possibly hold it just in one place via unique_ptr. Possibly do it via shared_ptr
class MinCostFlowFactor {
public:
   struct Edge {
      const INDEX start_node, end_node, cap; // we assume that the reverse capacity is 0
      const REAL cost;
   };

   MinCostFlowFactor(const std::vector<Edge>& edges, const std::vector<SIGNED_INDEX>& excess) 
   {
      const INDEX no_edges = edges.size();
      const INDEX no_nodes = excess.size();
      minCostFlow_ = new MinCost<SIGNED_INDEX,REAL>(no_nodes, no_edges);
      for(INDEX i=0; i<excess.size(); ++i) {
         minCostFlow_->AddNodeExcess(i, excess[i]);
      }
      assert(std::accumulate(excess.begin(), excess.end(),0) == 0);
      for(INDEX e=0; e<edges.size(); ++e) {
         assert(edges[e].start_node < no_nodes && edges[e].end_node < no_nodes);
         minCostFlow_->AddEdge(edges[e].start_node, edges[e].end_node, 0, edges[e].cap, edges[e].cost);
      }
   }

   ~MinCostFlowFactor()
   {}

   template<typename REPAM_ARRAY>
   void MaximizePotential(const REPAM_ARRAY& repam) { minCostFlow_->Solve();}

   template<typename REPAM_ARRAY>
   REAL LowerBound(const REPAM_ARRAY& repamPot) const
   {
      return repamPot.GetMinCostFlowSolver()->Solve(); // here we always assume that LowerBound is called with the current reparametrization
   }

   // do zrobienia: if shortest path is found, we can make provisions such that a subsequent call to UpdateCost with the Breakpoint value is computed fast.
   // we need only to call augment afterwards with minor modifications, all parent pointers are valid.
   // here we assume that capacity is [0,1]
   template<typename REPAM_ARRAY>
   const REAL GetBreakpoint(const REPAM_ARRAY& repamPot, const INDEX minCostFlowIdx)
   {
      const INDEX start_node = minCostFlow_->GetTailNodeId(minCostFlowIdx);
      const INDEX end_node = minCostFlow_->GetHeadNodeId(minCostFlowIdx);
      assert(minCostFlow_ == repamPot.GetMinCostFlowSolver());
      //const REAL cost = minCostFlow_->Solve(); // should not be needed
      const INDEX flow = 1-minCostFlow_->GetRCap(minCostFlowIdx);
      assert(flow == minCostFlow_->GetReverseRCap(minCostFlowIdx));
      assert(flow == 0 || flow == 1);
      if(flow == 0) { // compute cost of a shortest path from edge head to edge tail wrt. residual graph
         //std::cout << "Get breakpoint, flow is 0\n";
         const REAL breakpoint = -minCostFlow_->ShortestPath(end_node, start_node) - minCostFlow_->GetCost(minCostFlowIdx); // kwaskwas edge cost
         /*
         minCostFlow_->UpdateCost(minCostFlowIdx, 1, breakpoint + 1e-7);
         const REAL cost_1 = minCostFlow_->Solve();
         assert(cost_1 == cost);
         const INDEX flow_1 = 1-minCostFlow_->GetRCap(minCostFlowIdx);
         assert(flow_1 == flow);
         minCostFlow_->UpdateCost(minCostFlowIdx, 1, -2.0*1e-7);
         minCostFlow_->Solve();
         const INDEX flow_2 = 1-minCostFlow_->GetRCap(minCostFlowIdx);
         assert(flow_2 != flow);
         minCostFlow_->UpdateCost(minCostFlowIdx, 1, -breakpoint + 1e-7);
         minCostFlow_->Solve();
         */
         return breakpoint;
      } else { // compute cost of a shortest path from edge tail to edge head wrt. residual graph
         //std::cout << "Get breakpoint, flow is 1, edge cost is " << minCostFlow_->GetCost(minCostFlowIdx) << "\n";
         const REAL breakpoint = minCostFlow_->ShortestPath(start_node, end_node) - minCostFlow_->GetCost(minCostFlowIdx); // kwaskwas edge cost;
         /*
         //std::cout << "breakpoint is " << breakpoint << ", edge cost is " << minCostFlow_->GetCost(minCostFlowIdx) << "\n";
         minCostFlow_->UpdateCost(minCostFlowIdx, 1, breakpoint - 1e-7);
         const REAL cost_1 = minCostFlow_->Solve();
         assert(std::abs(cost_1 - (cost+breakpoint)) < 1e-6);
         const INDEX flow_1 = 1-minCostFlow_->GetRCap(minCostFlowIdx);
         assert(flow_1 == flow);
         minCostFlow_->UpdateCost(minCostFlowIdx, 1, +2.0*1e-7);
         minCostFlow_->Solve();
         const INDEX flow_2 = 1-minCostFlow_->GetRCap(minCostFlowIdx);
         assert(flow_2 != flow);
         minCostFlow_->UpdateCost(minCostFlowIdx, 1, -breakpoint - 1e-7);
         minCostFlow_->Solve();
         */
         return breakpoint;
      }
   }

   MinCost<SIGNED_INDEX,REAL>* GetMinCostFlowSolver() const
   {
      return minCostFlow_;
   }

   const INDEX size() const 
   { 
      throw std::runtime_error("not supperted"); return 0; 
      return minCostFlow_->GetEdgeNum();
   }
   const REAL operator[](const INDEX i) const 
   {
      throw std::runtime_error("not supperted"); return 0; 
      return minCostFlow_->GetCost(i);    
   }


   // for the SendMessages update step, compute maximal cost change such that current labeling stays optimal
   // do zrobienia: is possibly an array with different values from which to compute the maximal perturbation
   template<class REPAM_ARRAY>
   std::vector<REAL> MaximallyPerturbCosts(const REPAM_ARRAY& repam, const std::vector<bool>& active_edges) const
   {
      assert(active_edges.size() == minCostFlow_->GetEdgeNum());

      // do zrobienia: do not allocate anew each time but only once
      std::unique_ptr<MinCost<SIGNED_INDEX, REAL>> minCostFlowRepamUpdate_( new MinCost<SIGNED_INDEX,REAL>(minCostFlow_->GetNodeNum(), minCostFlow_->GetEdgeNum()) );
      
      const SIGNED_INDEX max_cap = 1000* minCostFlowRepamUpdate_->GetNodeNum();
      for(INDEX e=0; e<active_edges.size(); ++e) {
         assert(repam[e] == minCostFlow_->GetCost(e)); // for now: more general case requires storing primal solution
         if(active_edges[e] == true) {
            SIGNED_INDEX cur_min_cap, cur_max_cap;
            assert( minCostFlow_->GetRCap(e) == 1 || minCostFlow_->GetReverseRCap(e) == 1 );
            assert( minCostFlow_->GetRCap(e) == 0 || minCostFlow_->GetReverseRCap(e) == 0 );
            if(minCostFlow_->GetRCap(e) == 0) { // flow on edge is 1
               assert(minCostFlow_->GetFlow(e) == 1);
               cur_min_cap = -max_cap; 
               cur_max_cap = -1;
               minCostFlowRepamUpdate_->AddEdge(minCostFlow_->GetTailNodeId(e), minCostFlow_->GetHeadNodeId(e), cur_min_cap, cur_max_cap, minCostFlow_->GetCost(e));
            } else if(minCostFlow_->GetReverseRCap(e) == 0) { // flow on edge is 0
               assert(minCostFlow_->GetFlow(e) == 0);
               cur_min_cap = 1; 
               cur_max_cap = max_cap;
               minCostFlowRepamUpdate_->AddEdge(minCostFlow_->GetTailNodeId(e), minCostFlow_->GetHeadNodeId(e), cur_min_cap, cur_max_cap, minCostFlow_->GetCost(e));
               //cur_min_cap = 1; 
               //cur_max_cap = max_cap;
               //minCostFlowRepamUpdate_->AddEdge(minCostFlow_->GetHeadNodeId(e), minCostFlow_->GetTailNodeId(e), cur_min_cap, cur_max_cap, -minCostFlow_->GetCost(e));
            } else {
               throw std::runtime_error("Edge in min cost flow problem is not 0/1");
            }
         } else {
            const SIGNED_INDEX flow = minCostFlow_->GetFlow(e);
            if(flow == minCostFlow_->GetLowerBound(e)) {
               minCostFlowRepamUpdate_->AddEdge(minCostFlow_->GetTailNodeId(e), minCostFlow_->GetHeadNodeId(e), 0, max_cap, minCostFlow_->GetCost(e));
            } else if(flow == minCostFlow_->GetUpperBound(e)) {
               minCostFlowRepamUpdate_->AddEdge(minCostFlow_->GetTailNodeId(e), minCostFlow_->GetHeadNodeId(e), -max_cap, 0, minCostFlow_->GetCost(e));
            } else {
               minCostFlowRepamUpdate_->AddEdge(minCostFlow_->GetTailNodeId(e), minCostFlow_->GetHeadNodeId(e), -max_cap, max_cap, minCostFlow_->GetCost(e));
            }
         }
      }

      assert(minCostFlowRepamUpdate_->ExcessSum() == 0);
      std::cout << "Solving maximal reparametrization problem\n";
      minCostFlowRepamUpdate_->Solve();
      std::cout << "done\n";
      std::vector<REAL> repam_cost;
      repam_cost.reserve(std::count(active_edges.begin(),active_edges.end(),true));
      for(INDEX e=0; e<active_edges.size(); ++e) {
         if(active_edges[e] == true) {
            if(minCostFlow_->GetFlow(e) == 1) {
               repam_cost.push_back( std::max(minCostFlowRepamUpdate_->GetReducedCost(e), 0.0) );
               //repam_cost.push_back( -minCostFlowRepamUpdate_->GetReducedCost(e) );
               assert(minCostFlowRepamUpdate_->GetReducedCost(e) > -eps);
            } else if(minCostFlow_->GetFlow(e) == 0) {
               repam_cost.push_back( std::min(minCostFlowRepamUpdate_->GetReducedCost(e), 0.0) );
               //repam_cost.push_back( +minCostFlowRepamUpdate_->GetReducedCost(e) );
               assert(minCostFlowRepamUpdate_->GetReducedCost(e) < eps);
            } else {
               assert(false); // flow is not 0/1
            }
         }
      }
      // do zrobienia: construct test problem with perturbed cost and see whether solution is still optimal and just at edge
      return repam_cost;
   }

private:
   // note: this is also given to the reparametrization storage and hence points to the reparametrized potential. Use shared_ptr?
   MinCost<SIGNED_INDEX,REAL>* minCostFlow_;
};

template<typename FACTOR_CONTAINER>
class MinCostFlowReparametrizationStorage {
public:
   MinCostFlowReparametrizationStorage(const MinCostFlowFactor& mfc, const std::vector<REAL>& cost) 
      : minCostFlow_(mfc.GetMinCostFlowSolver())
   {
      assert(minCostFlow_ != nullptr);
   }

   ~MinCostFlowReparametrizationStorage() {
      static_assert( std::is_same<typename FACTOR_CONTAINER::FactorType, MinCostFlowFactor>::value, "");
      delete minCostFlow_;
   }

   MinCost<SIGNED_INDEX,REAL>* GetMinCostFlowSolver() const
   {
      return minCostFlow_;
   }

   const REAL operator[](const INDEX i) const 
   {
      return minCostFlow_->GetCost(i); // this could also be GetReducedCost(i) (?). But then this has to be done everywhere
   }

   // note: calling Solve() only necessary, if shortest path computation was done for updating min cost flow problem. But then again, we could solve mcf before calling shortest path computation.
   // calling Solve() is detrimental to performance for "reverse inverse linear program" reparametrization update.
   class MinCostWriteBack {
      public:
         MinCostWriteBack(MinCostFlowReparametrizationStorage* m, const INDEX dim) 
            : 
               minCostFlow_(m->GetMinCostFlowSolver()),
               dim_(dim)
         {}
         MinCostWriteBack& operator=(const REAL x)
         {
            const REAL diff = x - minCostFlow_->GetCost(dim_);
            assert(minCostFlow_->GetRCap(dim_) - minCostFlow_->GetReverseRCap(dim_) == 1);
            minCostFlow_->UpdateCost(dim_, 1, diff); // second parameter is the original capacity. We assume it is 1
            //minCostFlow_->Solve(); // do zrobienia: more efficient update might be possible after ShortestPath computation
            return *this;
         }
         MinCostWriteBack& operator+=(const REAL x)
         {
            assert(std::abs(minCostFlow_->GetRCap(dim_) - minCostFlow_->GetReverseRCap(dim_)) == 1);
            minCostFlow_->UpdateCost(dim_, 1, x); // second parameter is the original capacity. We assume it is 1
            //minCostFlow_->Solve(); // do zrobienia: more efficient update might be possible after ShortestPath computation
            return *this;

         }
         MinCostWriteBack& operator-=(const REAL x)
         {
            assert(std::abs(minCostFlow_->GetRCap(dim_) - minCostFlow_->GetReverseRCap(dim_)) == 1);
            minCostFlow_->UpdateCost(dim_, 1, -x); // second parameter is the original capacity. We assume it is 1
            //minCostFlow_->Solve(); // do zrobienia: more efficient update might be possible after ShortestPath computation
            return *this;
         }
      operator REAL() const { return minCostFlow_->GetCost(dim_); }
      private:
         MinCost<SIGNED_INDEX,REAL>* minCostFlow_;
         const INDEX dim_;
   };

   MinCostWriteBack operator[](const INDEX i) 
   {
      // build proxy object which automatically updates the cost in the flow graph
      return MinCostWriteBack(this,i);
   }

   const INDEX size() const 
   {
      return minCostFlow_->GetEdgeNum();
   }

private:
   // do zrobienia: what about shared_ptr?
   MinCost<SIGNED_INDEX,REAL>* minCostFlow_;
};

} // end namespace LP_MP

#endif // LP_MP_MIN_COST_FLOW_FACTOR_HXX
