#ifndef LP_MP_MIN_COST_FLOW_FACTOR_CS2_HXX
#define LP_MP_MIN_COST_FLOW_FACTOR_CS2_HXX

#include "mcmf.hxx"
#include "lib/MinCost/MinCost.h"
#include "config.hxx"

namespace LP_MP {


// mcf holds reverse copies of edges. arcs are ordered lexicographically by (tail,node).
// this creates problems for some applications, where only one direction of arcs is needed, but it may be advantageous for other applications.
class MinCostFlowFactorCS2 {
public:
   //using MinCostFlowSolverType = CS2_CPP::MCMF_CS2<>;
   using MinCostFlowSolverType = MCF::SSP<SIGNED_INDEX,REAL>;
   /*
   constexpr static REAL scalingFactor_ = 1000000000;
   // note: the scaling is brittle. Unless long unsigned int = 64 bit is used, monotone improvement often fails for eps in tolerance.hxx
   // do zrobienia: adaptive scaling?
   LONG_SIGNED_INDEX Scale(const REAL x) const { assert(std::abs(round(x)) < std::numeric_limits<LONG_SIGNED_INDEX>::max()/scalingFactor_ ); return std::round(scalingFactor_*x); }
   REAL Descale(const LONG_SIGNED_INDEX x) const { return REAL(x)/REAL(scalingFactor_); }
   */
   struct Edge {
      const INDEX start_node, end_node, lower_bound, upper_bound; // we assume that the minimum capacity is 0
      const REAL cost;
   };
   
   // provide method primal_size(), as not all arcs are [0,1]-flow and actually used as variables.

   MinCostFlowFactorCS2(const std::vector<Edge>& edges, const std::vector<SIGNED_INDEX>& supply_demand) 
   {
      INDEX noNodes_ = supply_demand.size();
      INDEX noEdges_ = edges.size(); // arc and reverse arc

      minCostFlow_ = new MinCostFlowSolverType(noNodes_, noEdges_);
      repamUpdateFlow_ = new MinCostFlowSolverType(noNodes_, noEdges_);
      for(auto& e : edges) {
         //minCostFlow_->set_arc(e.start_node, e.end_node, e.lower_bound, e.upper_bound, Scale(e.cost));
         //repamUpdateFlow_->set_arc(e.start_node, e.end_node, 0, 1, Scale(e.cost));
         minCostFlow_->AddEdge(e.start_node, e.end_node, e.lower_bound, e.upper_bound, e.cost);
         repamUpdateFlow_->AddEdge(e.start_node, e.end_node, 0, std::numeric_limits<SIGNED_INDEX>::max(), 0);
      }
      for(INDEX i=0; i<supply_demand.size(); ++i) {
         minCostFlow_->AddNodeExcess(i, supply_demand[i]);
         minCostFlow_->AddNodeExcess(i, 0);
         //minCostFlow_->set_supply_demand_of_node(i, supply_demand[i]);
         //repamUpdateFlow_->set_supply_demand_of_node(i,0);
      }
      minCostFlow_->SortArcs();
      minCostFlow_->Solve(); // to initialize data structures
      repamUpdateFlow_->SortArcs();
      //repamUpdateFlow_->run_cs2(); // to initialize data structures
   }

   ~MinCostFlowFactorCS2()
   {
      // delete minCostFlow_; deleted in repam storage: not very clean!
      // do zrobienia: memory leak
      //delete graph_; // one has to delete the graph at the end! However add copy constructor
   }

   template<typename REPAM_ARRAY>
   REAL EvaluatePrimal(const REPAM_ARRAY& repam, const PrimalSolutionStorage::Element primal) const
   {
      REAL cost = 0.0;
      //assert(repam.size() == arcs_.size());
      // check whether primal indeed belongs to a flow
      // the below is wrong whenever primal and repam have different sizes.
      for(INDEX e=0; e<repam.size(); e++) {
         if(minCostFlow_->GetCap(e) == 1) {
            assert(primal[e] == primal[ minCostFlow_->GetReverseArcId(e) ]);
            cost += repam[e]*primal[e]; 
         }
      }
      return cost;
   }

   void WritePrimal(const std::vector<INDEX>& flow, std::ofstream& fs) const
   {
      for(INDEX i=0; i<flow.size()-1; ++i) {
         fs << flow[i] << ", ";
      }
      fs << flow.back() << "\n";
   }

   template<typename REPAM_ARRAY>
   void MaximizePotential(const REPAM_ARRAY& repam)
   { 
      //for(INDEX e; e<repam.size(); ++e) {
      //   minCostFlow_->set_cost(e,repam[e]);
      //}
      //minCostFlow_->cs2_cost_restart();
      minCostFlow_->Solve();
      //for(INDEX e=0; e<arcs_.size(); ++e) {
      //   cost[arcs_[e]] = Scale(repam[e]);
      //}
      //auto ret = minCostFlow_->costMap(cost).run(); 
      //assert(ret == MinCostFlowSolverType::OPTIMAL);
   }

   template<typename REPAM_ARRAY>
   void MaximizePotentialAndComputePrimal(const REPAM_ARRAY& repam, typename PrimalSolutionStorage::Element primal)
   { 
      std::cout << "round mcf\n";
      // we assume here that repam is the current one as well
      //MaximizePotential(repam);
      //for(INDEX e; e<repam.size(); ++e) {
      //   minCostFlow_->set_cost(e,repam[e]);
      //}
      minCostFlow_->Solve();
      for(INDEX e=0; e<repam.size(); ++e) {
         // here we assume that flow is 0/1
         // but there may be edges that are not 0/1, e.g. slack edges, which are auxiliary ones
         //assert(-1 <= minCostFlow_->GetFlow(e) && minCostFlow_->GetFlow(e) <= 1);
         primal[e] = std::abs(minCostFlow_->GetFlow(e));
      }
   }

   // sometimes reverse edges have to be set as well
   void PropagatePrimal(PrimalSolutionStorage::Element primal)
   {
      for(INDEX a=0; a<size(); ++a) {
         if(primal[a] != unknownState) {
            const INDEX a_rev = minCostFlow_->GetReverseArcId(a);
            assert(primal[a_rev] == primal[a] || primal[a_rev] == unknownState);
            primal[a_rev] = primal[a];
         }
      }
   }

   template<typename REPAM_ARRAY>
   REAL LowerBound(const REPAM_ARRAY& repamPot) const
   {
      assert(repamPot.size() == size());
      return minCostFlow_->Solve();
      //return Descale(minCostFlow_->cs2_cost_restart()); // we assume here that current reparametrization is also repamPot
      //GraphType::ArcMap<LONG_SIGNED_INDEX> repam(*graph_);
      //for(INDEX e=0; e<arcs_.size(); ++e) {
      //   repam[arcs_[e]] = Scale(repamPot[e]);
      //}

      //auto ret = minCostFlow_->costMap(repam).run();
      //assert(ret == MinCostFlowSolverType::OPTIMAL);
      
      //return Descale(minCostFlow_->totalCost());
   }

   MinCostFlowSolverType* GetMinCostFlowSolver() const
   {
      return minCostFlow_;
   }
   MinCostFlowSolverType* GetRepamFlowSolver() const
   {
      return repamUpdateFlow_;
   }

   const INDEX size() const 
   { 
      //throw std::runtime_error("not supported"); return 0; 
      return minCostFlow_->GetArcNum();
   }
   const REAL operator[](const INDEX i) const 
   {
      return minCostFlow_->GetCost(i);
      throw std::runtime_error("not supported"); return 0; 
      return 0.0;
   }

   // for the SendMessages update step, compute maximal cost change such that current labeling stays optimal
   // do zrobienia: is possibly an array with different values from which to compute the maximal perturbation
   // do zrobienia: rename active_arcs
   template<class REPAM_ARRAY, typename ACTIVE_EDGES_ARRAY>
   void MaximallyPerturbCosts(const REPAM_ARRAY& repam, const ACTIVE_EDGES_ARRAY& active_edges) const
   {
      // first read reduced cost of mcf into repam update
      for(INDEX i=0; i<minCostFlow_->GetArcNum(); ++i) {
         if(minCostFlow_->GetTailNodeId(i) < minCostFlow_->GetHeadNodeId(i)) { // avoid copying twice
            repamUpdateFlow_->SetCost(i,minCostFlow_->GetCost(i));
         }
      }
      // set capacities
      assert(active_edges.size() == minCostFlow_->GetArcNum());
      for(INDEX i=0; i<minCostFlow_->GetArcNum(); ++i) {
         if(active_edges[i]) {
            if(minCostFlow_->GetRCap(i) == 0) {
               repamUpdateFlow_->SetCap(i,-1);
            } else {
               assert(minCostFlow_->GetRCap(i) == 1);
               // beware of overflows!
               repamUpdateFlow_->SetCap(i,std::numeric_limits<SIGNED_INDEX>::max()/16);
            }
         } else {
            repamUpdateFlow_->SetCap(i,std::numeric_limits<SIGNED_INDEX>::max()/16);
         }
      }
      repamUpdateFlow_->Solve();
      /*
      // we assume here that repam is the current one
      assert(active_edges.size() == size());
      const SIGNED_INDEX max_cap = 100000;
      for(INDEX a=0; a<minCostFlow_->no_arcs(); ++a) {
         // do zrobienia: make checks in terms of residual capacity
         const SIGNED_INDEX flow = minCostFlow_->get_flow(a);
         if(active_edges[a] == true) {
            if(flow == 0) { 
               repamUpdateFlow_->set_cap(a,max_cap);
            } else if(flow == 1) {
               repamUpdateFlow_->set_cap(a,-1);
            } else if(flow == -1) { // reverse edge, also treat it
               repamUpdateFlow_->set_cap(a,1);
            } else {
               assert(false);
            }
         } else {
            repamUpdateFlow_->set_cap(a,max_cap);
         }
      }
      // copy cost
      repamUpdateFlow_->copy_cost(minCostFlow_);
      // run
      repamUpdateFlow_->cs2_cost_restart();
      // read out dual variables which will go into reparametrization update
      std::vector<REAL> repam_cost;
      INDEX noActiveEdges = 0;
      for(INDEX e=0; e<active_edges.size(); ++e) {
         noActiveEdges += active_edges[e];
      }
      repam_cost.reserve(noActiveEdges);
      for(INDEX e=0; e<active_edges.size(); ++e) {
         if(active_edges[e] == true) {
            // do zrobienia: why should this be true?
            //const REAL reducedCost = Descale(Scale(minCostFlowSolver_->get_cost(e) + repamUpdateFlow_->.potential(graph_->source(arcs_[e])) - minCostFlowRepamUpdate.potential(graph_->target(arcs_[e])));
            const REAL reducedCost = Descale(repamUpdateFlow_->get_reduced_cost(e));
            if(minCostFlow_->get_flow(e) == 1) {
               repam_cost.push_back( -reducedCost );
            } else if(minCostFlow_->get_flow(e) == 0) {
               repam_cost.push_back( -reducedCost );
            } else if(minCostFlow_->get_flow(e) == -1) {
               repam_cost.push_back( -reducedCost );
            } else {
               assert(false); // flow is not 0/1
            }
         }
      }
      return repam_cost;
      */

      /*
      GraphType::ArcMap<LONG_SIGNED_INDEX> cost(*graph_);
      for(INDEX e=0; e<active_edges.size(); ++e) {
         const SIGNED_INDEX flow = minCostFlow_->flow(arcs_[e]);
         if(active_edges[e] == true) {
            assert( flow == 0 || flow == 1); 
            if(flow == 1) {
               lower[arcs_[e]] = -max_cap;
               upper[arcs_[e]] = -1;
            } else if(flow == 0) {
               lower[arcs_[e]] = 1;
               upper[arcs_[e]] = max_cap;
            }
         } else {
            if(flow == 0) {
               lower[arcs_[e]] = 0;
               upper[arcs_[e]] = +max_cap;
            } else {
               lower[arcs_[e]] = -max_cap;
               upper[arcs_[e]] = +max_cap;
            }
         }
         assert(lower[arcs_[e]] < upper[arcs_[e]]);
         cost[arcs_[e]] = Scale(repam[e]);
      }
      auto ret = minCostFlowRepamUpdate.lowerMap(lower).upperMap(upper).costMap(cost).run();
      assert(ret == MinCostFlowSolverType::OPTIMAL);
      if(ret != MinCostFlowSolverType::OPTIMAL) {
         throw std::runtime_error("Reparametrization problem could not be solved to optimality");
      }

      std::vector<REAL> repam_cost;
      INDEX noActiveEdges = 0;
      for(INDEX e=0; e<active_edges.size(); ++e) {
         noActiveEdges += active_edges[e];
      }
      repam_cost.reserve(noActiveEdges);
      for(INDEX e=0; e<active_edges.size(); ++e) {
         if(active_edges[e] == true) {
            // do zrobienia: why should this be true?
            const REAL reducedCost = Descale(Scale(repam[e]) + minCostFlowRepamUpdate.potential(graph_->source(arcs_[e])) - minCostFlowRepamUpdate.potential(graph_->target(arcs_[e])));
            if(minCostFlow_->flow(arcs_[e]) == 1) {
               repam_cost.push_back( -reducedCost );
            } else if(minCostFlow_->flow(arcs_[e]) == 0) {
               repam_cost.push_back( -reducedCost );
            } else {
               assert(false); // flow is not 0/1
            }
         }
      }
      // do zrobienia: construct test problem with perturbed cost and see whether solution is still optimal and just at edge
      return repam_cost;
      */
   }

  void CreateConstraints(LpInterfaceAdapter* lp) const { 
     // go through all edges and link let forward be 1-reverse edge
     for(INDEX a=0; a<minCostFlow_->GetArcNum(); ++a) {
        const INDEX i = minCostFlow_->GetTailNodeId(a);
        const INDEX j = minCostFlow_->GetHeadNodeId(a);
        assert(i != j);
        if(i < j) {
           const INDEX a_rev = minCostFlow_->GetReverseArcId(a);
           assert( minCostFlow_->GetCap(a) == 1 && minCostFlow_->GetCap(a_rev) == 0
                 || minCostFlow_->GetCap(a) == 0 && minCostFlow_->GetCap(a_rev) == 1);
           LinExpr lhs = lp->CreateLinExpr();
           lhs += lp->GetVariable(a);
           lhs -= lp->GetVariable(a_rev); // this is not the flow formulation, but one based on "undirected flow"
           LinExpr rhs = lp->CreateLinExpr();
           lp->addLinearEquality(lhs,rhs);
        }
     }
     // go through all nodes and add flow conservation constraints
     for(INDEX i=0; i<minCostFlow_->GetNodeNum(); ++i) {
        LinExpr lhs = lp->CreateLinExpr();
        for(INDEX a=minCostFlow_->StartingArc(i); a<minCostFlow_->StartingArc(i) + minCostFlow_->NoArcs(i); ++a) {
           lhs += lp->GetVariable(a);
        }
        LinExpr rhs = lp->CreateLinExpr();
        rhs += minCostFlow_->GetDemand(i); // do zrobienia: store this directly and not in min cost flow problem
        lp->addLinearEquality(lhs,rhs);
     }
  }

private:
   // note: this is also given to the reparametrization storage and hence points to the reparametrized potential. Use shared_ptr?
   MinCostFlowSolverType* minCostFlow_;
   MinCostFlowSolverType* repamUpdateFlow_;
};

template<typename FACTOR_CONTAINER>
class MinCostFlowReparametrizationStorageCS2 {
public:
   /*
   constexpr static REAL scalingFactor_ = MinCostFlowFactorCS2::scalingFactor_;
   // note: the scaling is brittle. Unless long unsigned int = 64 bit is used, monotone improvement often fails for eps in tolerance.hxx
   // do zrobienia: adaptive scaling?
   static LONG_SIGNED_INDEX Scale(const REAL x) { assert(std::abs(round(x)) < std::numeric_limits<LONG_SIGNED_INDEX>::max()/scalingFactor_ ); return std::round(scalingFactor_*x); }
   static REAL Descale(const LONG_SIGNED_INDEX x) { return REAL(x)/REAL(scalingFactor_); }
   */

   using MinCostFlowSolverType = MinCostFlowFactorCS2::MinCostFlowSolverType;
   MinCostFlowReparametrizationStorageCS2(const MinCostFlowFactorCS2& mfc) 
      : minCostFlow_(mfc.GetMinCostFlowSolver()),
      repamUpdateFlow_(mfc.GetRepamFlowSolver())
   {
      assert(minCostFlow_ != nullptr);
   }

   ~MinCostFlowReparametrizationStorageCS2() {
      static_assert( std::is_same<typename FACTOR_CONTAINER::FactorType, MinCostFlowFactorCS2>::value, "");
      delete minCostFlow_;
   }
   
   //possibly not needed here.
   MinCostFlowSolverType* GetMinCostFlowSolver() const
   {
      return minCostFlow_;
   }

   const REAL operator[](const INDEX i) const 
   {
      return minCostFlow_->GetCost(i);
   }

   class WriteBackProxy {
      public:
      WriteBackProxy(const INDEX i, MinCostFlowSolverType* mcf) 
         : i_(i),
         minCostFlow_(mcf)
      {}
      WriteBackProxy operator=(const REAL x) {
         minCostFlow_->SetCost(i_,x);
         return *this;
      }
      operator REAL() const { 
        return minCostFlow_->GetCost(i_);
      }

      // when we write to the mcf factor, the cost of the reverse arc is also modified
      // do zrobienia: introduce function for changing cost by constant
      WriteBackProxy operator+=(const REAL x) {
         // sign is necessary, as we may write cost to reverse edges -> reverse costs
         const REAL sign = minCostFlow_->GetCap(i_) == 1 ? 1.0 : -1.0;
         minCostFlow_->UpdateCost(i_, sign*x);
         return *this;
      }

      WriteBackProxy operator-=(const REAL x) {
         const REAL sign = minCostFlow_->GetCap(i_) == 1 ? 1.0 : -1.0;
         minCostFlow_->UpdateCost(i_, sign*(-x)); // do zrobienia: capacity is possibly wrong! reverse has cap 0
         return *this;
      }

      private:
      const INDEX i_;
      MinCostFlowSolverType* minCostFlow_;
   };

   WriteBackProxy operator[](const INDEX i) 
   {
      return WriteBackProxy(i, minCostFlow_);
   }

   const INDEX size() const 
   {
      return minCostFlow_->GetArcNum();
   }

private:
   // do zrobienia: what about shared_ptr?
   MinCostFlowSolverType* minCostFlow_; //possibly not needed here.
   MinCostFlowSolverType* repamUpdateFlow_;
};


// used in graph matching via mcf. 
// do zrobienia: rename
template<INDEX COVERING_FACTOR, Chirality PRIMAL_ROUNDING_DIRECTION> // COVERING_FACTOR specifies how often an mcf edge is covered at most. Usually, either one or two
class UnaryToAssignmentMessageCS2 {
// right factor is MinCostFlowFactor describing an 1:1 assignment as constructed above
// left factors are the left and right simplex corresponding to it. 
// Those are ordered as the left and right nodes
// templatize edgeIndex_ to allow for more compact representation. Possibly hold reference to original edgeId structure instead of full vector
public:
   UnaryToAssignmentMessageCS2(const INDEX start_arc, const INDEX no_arcs)
      : 
         start_arc_(start_arc),
         no_arcs_(no_arcs)
   {
      //assert(no_arcs >= 2);
   }

   template<typename LEFT_POT, typename MSG_ARRAY>
   void MakeLeftFactorUniform(const LEFT_POT& leftPot, MSG_ARRAY& msg, const REAL omega = 1.0) 
   {
      assert(leftPot.size() == no_arcs_ && msg.size() == no_arcs_);

      for(INDEX i=0; i<leftPot.size(); ++i) {
         msg[i] -= omega*(leftPot[i]);
      }
   }

   template<typename RIGHT_FACTOR, typename RIGHT_POT, typename MSG_ARRAY>
   void MakeRightFactorUniform(const RIGHT_FACTOR& rightFactor, const RIGHT_POT& perturb, MSG_ARRAY& msg, const REAL omega = 1.0) const 
   {
      assert(msg.size() == no_arcs_);
      for(INDEX i=0; i<no_arcs_; ++i) {
         const INDEX curEdgeId = start_arc_ + i;
         msg[i] -= omega*perturb[curEdgeId]; 
      }
   }

   // performance is enhanced when not using this. 
   template<typename LEFT_FACTOR, typename G1, typename G2>
   void ReceiveMessageFromLeft(LEFT_FACTOR* l, const G1& leftPot, G2& msg)
   {
      MakeLeftFactorUniform(leftPot, msg);
   }

   // for hungarian BP
   template<typename LEFT_FACTOR, typename G1, typename G2>
   void
   ReceiveRestrictedMessageFromLeft(LEFT_FACTOR* l, const G1& leftPot, G2& msg, typename PrimalSolutionStorage::Element leftPrimal)
   { 
      // we assume that leftPrimal is all unknownState
      for(INDEX i=0; i<no_arcs_; ++i) {
         assert(leftPrimal[i] == unknownState);
      }
      MakeLeftFactorUniform(leftPot, msg);
   }

   // for rounding based on unaries
   template<typename RIGHT_FACTOR, typename G1, typename G2>
   void
   ReceiveRestrictedMessageFromRight(RIGHT_FACTOR* r, const G1& rightPot, G2& msg, typename PrimalSolutionStorage::Element rightPrimal)
   { 
      auto* mcf = r->GetMinCostFlowSolver();
      // if an arc is true, enforce it. If an arc is false, forbid it. Otherwise, get reduced costs
      for(INDEX l=0; l<no_arcs_; ++l) {
         //msgs[i][l] -= omega_sum*1.0/REAL(COVERING_FACTOR)*(-mcf->GetReducedCost(start_arc + l) + mcf->GetCost(start_arc + l));
         const INDEX e = start_arc_ + l;
         assert(mcf->GetCap(e) == 0 || mcf->GetCap(e) == 1); // only {0,1}-flow accepted here!
         const REAL sign = mcf->GetCap(e) == 1 ? 1.0 : -1.0;
         if(rightPrimal[e] == unknownState) {
            msg[l] -= (sign*mcf->GetReducedCost(e) );  // do zrobienia: this only works for left side models, not for right hand side ones!
         } else if(rightPrimal[e] == true) {
            msg[l] -= -std::numeric_limits<REAL>::max();
         } else {
            assert(rightPrimal[e] == false);
            msg[l] -= std::numeric_limits<REAL>::max();
         }
      }
   }
   /*
   template<typename LEFT_FACTOR, typename RIGHT_FACTOR, typename G1, typename G2, typename G3>
   void SendMessageToRight(LEFT_FACTOR* const l, RIGHT_FACTOR* const r, const G1& leftPot, const G2& rightPot, G3& msg, const REAL omega)
   {
      MakeLeftFactorUniform(leftPot, msg, omega);
   }
   */

   template<typename RIGHT_FACTOR, typename MSG_ARRAY, typename RIGHT_REPAM, typename ITERATOR>
   static void
   SendMessagesToLeft(const RIGHT_FACTOR& rightFactor, const RIGHT_REPAM& rightRepam, const MSG_ARRAY& msgs, ITERATOR omegaIt)
   {
      auto* mcf = rightFactor.GetMinCostFlowSolver();
      assert(rightRepam.size() == rightFactor.size());
      for(INDEX i=0; i<rightFactor.size(); ++i) {
         assert(std::abs(rightRepam[i] - mcf->GetCost(i)) <= eps);
      }
      const REAL omega_sum = std::accumulate(omegaIt, omegaIt+msgs.size(),0.0); 
      assert(0.0 <= omega_sum && omega_sum < 1.0 + eps);
      /*
      std::vector<bool> active_edges(rightFactor.size(),false);
      // iterate over all edges covered by arcs and set to true
      for(INDEX m=0; m<msgs.size(); ++m) {
         auto& msg = msgs[m].GetMessageOp();
         // note: start_arc_ must be ordered
         for(INDEX i=0; i<msg.size(); ++i) {
            active_edges[msg.start_arc_+i] = true;
         }
      }
      std::vector<REAL> repam = rightFactor.MaximallyPerturbCosts(rightRepam, active_edges);
      */

      // the technique applied in Max-Weight Bipartite Matching ... CVPR16
      for(INDEX i=0; i<msgs.size(); ++i, ++omegaIt) {
         INDEX start_arc = msgs[i].GetMessageOp().start_arc_;
         INDEX no_arcs = msgs[i].GetMessageOp().no_arcs_;
         for(INDEX l=0; l<no_arcs; ++l) {
            //msgs[i][l] -= omega_sum*1.0/REAL(COVERING_FACTOR)*(-mcf->GetReducedCost(start_arc + l) + mcf->GetCost(start_arc + l));
            const INDEX e = start_arc + l;
            assert(mcf->GetCap(e) == 0 || mcf->GetCap(e) == 1); // only {0,1}-flow accepted here!
            const REAL sign = mcf->GetCap(e) == 1 ? 1.0 : -1.0;
            msgs[i][l] -= omega_sum*1.0/REAL(COVERING_FACTOR)*(sign*mcf->GetReducedCost(e) );  // do zrobienia: this only works for left side models, not for right hand side ones!
         }
         //msgs[i].GetMessageOp().MakeRightFactorUniform(rightFactor, repam, msgs[i], omega_sum*1.0/REAL(COVERING_FACTOR));
      }
   }

   // do zrobienia: possibly automatically query solution type of min cost flow factor and set type of right argument automatically

   template<bool ENABLE = (PRIMAL_ROUNDING_DIRECTION == Chirality::right), typename LEFT_FACTOR, typename RIGHT_FACTOR>
   typename std::enable_if<ENABLE,void>::type
   ComputeLeftFromRightPrimal(const typename PrimalSolutionStorage::Element left, LEFT_FACTOR* l, typename PrimalSolutionStorage::Element right, RIGHT_FACTOR* r)
   {
      for(INDEX e=0; e<no_arcs_; ++e) {
         left[e] = right[start_arc_ + e];
      }
   }

   template<bool ENABLE = (PRIMAL_ROUNDING_DIRECTION == Chirality::left), typename LEFT_FACTOR, typename RIGHT_FACTOR>
   typename std::enable_if<ENABLE,void>::type
   ComputeRightFromLeftPrimal(typename PrimalSolutionStorage::Element left, LEFT_FACTOR* l, const typename PrimalSolutionStorage::Element right, RIGHT_FACTOR* r)
   {
      auto* mcf = r->GetMinCostFlowSolver();
      for(INDEX e=0; e<no_arcs_; ++e) {
         right[start_arc_ + e] = left[e];
         const INDEX a_rev = mcf->GetReverseArcId(start_arc_ + e);
         right[a_rev] = left[e];
      }
   }

   INDEX size() const { return no_arcs_; }

   template<typename G>
   void RepamLeft(G& leftRepamPot, const REAL msg, const INDEX dim) 
   {
      leftRepamPot[dim] += msg;
   }

   template<typename G>
   void RepamRight(G& rightRepamPot, const REAL msg, const INDEX dim) 
   {
      rightRepamPot[start_arc_ + dim] += msg;
   }

   template<class LEFT_FACTOR_TYPE,class RIGHT_FACTOR_TYPE>
   void CreateConstraints(LpInterfaceAdapter* lp,LEFT_FACTOR_TYPE* LeftFactor,RIGHT_FACTOR_TYPE* RightFactor) const
   { 
      for(INDEX i=0; i<no_arcs_; ++i) {
         LinExpr lhs = lp->CreateLinExpr();
         LinExpr rhs = lp->CreateLinExpr();
         lhs += lp->GetLeftVariable(i);
         rhs += lp->GetRightVariable(i + start_arc_);
         lp->addLinearEquality(lhs,rhs);
      }
   }


private:
   const INDEX start_arc_;
   const INDEX no_arcs_;
};

} // end namespace LP_MP

#endif // LP_MP_MIN_COST_FLOW_FACTOR_CS2_HXX


