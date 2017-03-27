#ifndef LP_MP_MIN_COST_FLOW_FACTOR_CS2_HXX
#define LP_MP_MIN_COST_FLOW_FACTOR_CS2_HXX

//#include "mcmf.hxx"
#include "lib/MinCost/MinCost.h"
#include "config.hxx"

//do zrobienia: remove again
//#include <ilcplex/ilocplex.h>

namespace LP_MP {

enum class MessagePassingType {SRMP,HUNGARIAN}; // to distinguish between hungarian bp and our variant

// mcf holds reverse copies of edges. arcs are ordered lexicographically by (tail,node).
// this creates problems for some applications, where only one direction of arcs is needed, but it may be advantageous for other applications.
class MinCostFlowFactorCS2 {
public:
   using MinCostFlowSolverType = MinCost<SIGNED_INDEX,REAL>;
   struct Edge {
      INDEX start_node, end_node, lower_bound, upper_bound; // we assume that the minimum capacity is 0
      REAL cost;
   };
   
   MinCostFlowFactorCS2(const std::vector<Edge>& edges, const std::vector<SIGNED_INDEX>& supply_demand, const INDEX no_binary_edges) 
      : demand_(supply_demand),
      no_binary_edges_(no_binary_edges)
   {
      INDEX noNodes_ = supply_demand.size();
      INDEX noEdges_ = edges.size();

      minCostFlow_ = new MinCostFlowSolverType(noNodes_, noEdges_);
      repamUpdateFlow_ = new MinCostFlowSolverType(noNodes_, noEdges_);
      for(auto& e : edges) {
         minCostFlow_->AddEdge(e.start_node, e.end_node, e.upper_bound, e.lower_bound, e.cost);
         repamUpdateFlow_->AddEdge(e.start_node, e.end_node, std::numeric_limits<SIGNED_INDEX>::max(), 0, 0);
      }
      for(INDEX i=0; i<supply_demand.size(); ++i) {
         minCostFlow_->AddNodeExcess(i, supply_demand[i]);
      }
      //minCostFlow_->SortArcs();
      minCostFlow_->Solve(); // to initialize data structures
      primal_.resize(edges.size());
      //repamUpdateFlow_->SortArcs();
      //repamUpdateFlow_->run_cs2(); // to initialize data structures
   }

   ~MinCostFlowFactorCS2()
   {
      // delete minCostFlow_; deleted in repam storage: not very clean!
      // do zrobienia: memory leak
      //delete graph_; // one has to delete the graph at the end! However add copy constructor
   }

   REAL EvaluatePrimal() const
   {
      assert(primal_.size() == this->size());
      REAL cost = 0.0;
      for(INDEX i=0; i<primal_.size(); ++i) {
         cost += primal_[i]*minCostFlow_->GetCost(i);
      }
      return cost;
      // do zrobienia: we must check feasiblity of primal as regards flow constraints. For this, a feasible flow on auxiliary variables must be built.
      // check feasibility
      /*
      std::vector<SIGNED_INDEX> excess(minCostFlow_->GetNodeNum(),0);
      for(INDEX e=0; e<repam.size(); ++e) {
         const INDEX tail = minCostFlow_->GetTailNodeId(e);
         const INDEX head = minCostFlow_->GetHeadNodeId(e);
         excess[tail] += primal[e];
         excess[head] -= primal[e];
      }
      for(INDEX i=0; i<excess.size(); ++i) {
         if(excess[i] != demand_[i]) {
            return std::numeric_limits<REAL>::infinity();
         }
      }
      */
      // to do: this is a very hacky implementation and it need not be correct for anything except assignment problems!
      // first check whether primal belongs to a feasible flow
      /*
      for(INDEX i=0; i<minCostFlow_->GetNodeNum(); ++i) {
         SIGNED_INDEX excess = 0;
         INDEX a_idx = minCostFlow_->StartingArc(i);
         for(INDEX c=0; c<minCostFlow_->NoArcs(i); ++c, ++a_idx) {
            if(primal[a_idx] == unknownState) {
               excess = minCostFlow_->GetDemand(i);
               continue;
            }
            const SIGNED_INDEX sign = minCostFlow_->GetCap(a_idx) == 1 ? 1 : -1;
            excess += sign*SIGNED_INDEX(primal[a_idx]);
         }
         if(excess != minCostFlow_->GetDemand(i)) {
            std::cout << "excess = " << excess << " demand = " << minCostFlow_->GetDemand(i) << "\n";
            return std::numeric_limits<REAL>::infinity();
         }
      }
      REAL cost = 0.0;
      assert(repam.size() == minCostFlow_->GetArcNum());
      for(INDEX e=0; e<repam.size(); e++) {
         if(minCostFlow_->GetCap(e) != 0 && minCostFlow_->GetCap(e) != 1) {
            assert(primal[e] == false || repam[e] == 0.0);
         }
         if(minCostFlow_->GetCap(e) == 1 && repam[e] != 0.0) { // some edges are auxiliary and have cost zero. Primal value is not set for them.
            assert(minCostFlow_->GetCap( minCostFlow_->GetReverseArcId(e) ) == 0);
            assert(primal[e] == primal[ minCostFlow_->GetReverseArcId(e) ]);
            assert(primal[e] != unknownState);
            cost += repam[e]*primal[e]; 
         }
      }
      return 0.5*cost;
      */
   }

   void MaximizePotential()
   { 
      minCostFlow_->Solve();
      //for(INDEX e; e<repam.size(); ++e) {
      //   minCostFlow_->set_cost(e,repam[e]);
      //}
      //minCostFlow_->cs2_cost_restart();
      //for(INDEX e=0; e<arcs_.size(); ++e) {
      //   cost[arcs_[e]] = Scale(repam[e]);
      //}
      //auto ret = minCostFlow_->costMap(cost).run(); 
      //assert(ret == MinCostFlowSolverType::OPTIMAL);
   }

   void MaximizePotentialAndComputePrimal()
   { 
      std::cout << "round mcf\n";
      std::cout << "problem is with infinities!\n";
      // we assume here that repam is the current one as well
      //MaximizePotential(repam);
      //for(INDEX e; e<repam.size(); ++e) {
      //   minCostFlow_->set_cost(e,repam[e]);
      //}
      minCostFlow_->Solve();
      for(INDEX e=0; e<primal_.size(); ++e) {
         // here we assume that flow is 0/1
         // but there may be edges that are not 0/1, e.g. slack edges, which are auxiliary ones
         //assert(-1 <= minCostFlow_->GetFlow(e) && minCostFlow_->GetFlow(e) <= 1);
         primal_[e] = minCostFlow_->GetFlow(e);
         //std::cout << INDEX(primal_[e]) << ",";
      }
      //std::cout << std::endl;
   }

   // sometimes reverse edges have to be set as well
   void PropagatePrimal(PrimalSolutionStorage::Element primal)
   {
      return;
      assert(false); // do zrobiebia: remove
      assert(false);
      return; // gets called too often when unaries round primals!
      /*
      for(INDEX a=0; a<size(); ++a) {
         if(primal[a] != unknownState) {
            const INDEX a_rev = minCostFlow_->GetReverseArcId(a);
            assert(primal[a_rev] == primal[a] || primal[a_rev] == unknownState);
            primal[a_rev] = primal[a];
         }
      }
      */
   }

   REAL LowerBound() const
   {
      return minCostFlow_->Solve();
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
      return no_binary_edges_;
   }
   const REAL operator[](const INDEX i) const 
   {
      assert(i<no_binary_edges_);
      return minCostFlow_->GetCost(i);
   }

   // for the SendMessages update step, compute maximal cost change such that current labeling stays optimal
   // do zrobienia: is possibly an array with different values from which to compute the maximal perturbation
   // do zrobienia: rename active_arcs
   template<class REPAM_ARRAY, typename ACTIVE_EDGES_ARRAY>
   void MaximallyPerturbCosts(const REPAM_ARRAY& repam, const ACTIVE_EDGES_ARRAY& active_edges) const
   {
      assert(false);
      // first read reduced cost of mcf into repam update
      /*
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
      */
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
      )return repam_cost;
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

   INDEX GetNumberOfAuxVariables() const { return minCostFlow_->GetEdgeNum() - no_binary_edges_; }

  void CreateConstraints(LpInterfaceAdapter* lp) const {
     // the first no_binary_edges go into usual variables, the last ones into auxiliary
     for(INDEX e=no_binary_edges_; e<minCostFlow_->GetEdgeNum(); ++e) {
        auto var = lp->GetAuxVariable(e-no_binary_edges_);
        lp->SetVariableBound(var, -REAL(minCostFlow_->GetReverseCap(e)), REAL(minCostFlow_->GetCap(e)));
     }
     // flow conservation constraints
     std::vector<LinExpr> flow_conservation(minCostFlow_->GetNodeNum());
     for(auto& c : flow_conservation) {
        c = lp->CreateLinExpr();
     }
     for(INDEX e=0; e<no_binary_edges_; ++e) {
        auto var = lp->GetVariable(e);
        const INDEX i=minCostFlow_->GetTailNodeId(e);
        const INDEX j=minCostFlow_->GetHeadNodeId(e);
        flow_conservation[i] += var;
        flow_conservation[j] -= var;
     }
     for(INDEX e=no_binary_edges_; e<minCostFlow_->GetEdgeNum(); ++e) {
        auto var = lp->GetAuxVariable(e-no_binary_edges_);
        const INDEX i=minCostFlow_->GetTailNodeId(e);
        const INDEX j=minCostFlow_->GetHeadNodeId(e);
        flow_conservation[i] += var;
        flow_conservation[j] -= var;
     }
     for(INDEX i=0; i<minCostFlow_->GetNodeNum(); ++i) {
        auto rhs = lp->CreateLinExpr();
        rhs += demand_[i]; // do zrobienia: store this directly and not in min cost flow problem
        lp->addLinearEquality(flow_conservation[i],rhs);
     }
  }

  void init_primal() { std::fill(primal_.begin(), primal_.end(), 2); }
  template<class ARCHIVE> void serialize_primal(ARCHIVE& ar) { ar(primal_); }
  template<class ARCHIVE> void serialize_dual(ARCHIVE& ar) 
  { 
     assert(false);
  }

   std::vector<unsigned char> primal_;
private:
   // note: this is also given to the reparametrization storage and hence points to the reparametrized potential. Use shared_ptr?
   MinCostFlowSolverType* minCostFlow_;
   MinCostFlowSolverType* repamUpdateFlow_;
   const std::vector<SIGNED_INDEX> demand_;
   const INDEX no_binary_edges_;


};

/*
template<typename FACTOR_CONTAINER>
class MinCostFlowReparametrizationStorageCS2 {
public:
   using MinCostFlowSolverType = MinCostFlowFactorCS2::MinCostFlowSolverType;
   MinCostFlowReparametrizationStorageCS2(const MinCostFlowFactorCS2& mfc) 
      : minCostFlow_(mfc.GetMinCostFlowSolver()),
      repamUpdateFlow_(mfc.GetRepamFlowSolver()),
      no_binary_edges_(mfc.size())
   {
      assert(minCostFlow_ != nullptr);
   }

   ~MinCostFlowReparametrizationStorageCS2() {
      std::cout << "delete min cost flow factor\n";
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
      {
         assert(i < minCostFlow_->GetEdgeNum());
      }
      WriteBackProxy operator=(const REAL x) {
         assert(std::isfinite(x));
         minCostFlow_->SetCost(i_, x);
         return *this;
      }
      operator REAL() const { 
         return minCostFlow_->GetCost(i_);
      }

      // when we write to the mcf factor, the cost of the reverse arc is also modified
      // do zrobienia: introduce function for changing cost by constant
      WriteBackProxy operator+=(const REAL x) {
         minCostFlow_->UpdateCost(i_, x);
         return *this;
      }

      WriteBackProxy operator-=(const REAL x) {
         minCostFlow_->UpdateCost(i_, -x);
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
      return no_binary_edges_;
   }

private:
   // do zrobienia: what about shared_ptr?
   MinCostFlowSolverType* minCostFlow_; //possibly not needed here.
   MinCostFlowSolverType* repamUpdateFlow_;
   const INDEX no_binary_edges_;
};
*/


// used in graph matching via mcf. 
// do zrobienia: rename
template<INDEX COVERING_FACTOR, MessagePassingType MPT> // COVERING_FACTOR specifies how often an mcf edge is covered at most. Usually, either one or two
class UnaryToAssignmentMessageCS2 {
// right factor is MinCostFlowFactor describing an 1:1 assignment as constructed above
// left factors are the left and right simplex corresponding to it. 
// Those are ordered as the left and right nodes
// templatize edgeIndex_ to allow for more compact representation. Possibly hold reference to original edgeId structure instead of full vector
public:
   UnaryToAssignmentMessageCS2(const std::vector<INDEX>& edges) : edges_(edges) {
      assert(edges_.size() >= 2);
      for(auto e : edges) {
         std::cout << e << ",";
      } 
      std::cout << std::endl;
   }

   template<typename LEFT_POT, typename MSG_ARRAY>
   void MakeLeftFactorUniform(const LEFT_POT& leftPot, MSG_ARRAY& msg, const REAL omega = 1.0) 
   {
      assert(leftPot.size() == edges_.size());

      for(INDEX i=0; i<leftPot.size(); ++i) {
         msg[i] -= omega*(leftPot[i]);
      }
   }

   // for hungarian bp, we do not invoke this function
   template<typename LEFT_FACTOR, typename LEFT_REPAM, typename MSG, bool ENABLE = MPT==MessagePassingType::SRMP >
   typename std::enable_if<ENABLE,void>::type
   SendMessageToRight(LEFT_FACTOR* l, const LEFT_REPAM& leftRepam, MSG& msg, const REAL omega)
   {
      MakeLeftFactorUniform(leftRepam, msg, omega);
   }
   
   template<typename LEFT_FACTOR, typename G1, typename G2>
   void ReceiveMessageFromLeft(LEFT_FACTOR* l, const G1& leftPot, G2& msg)
   {
      //assert(leftPot.size() == no_arcs_);
      MakeLeftFactorUniform(leftPot, msg);
   }

   // for hungarian BP
   template<typename LEFT_FACTOR, typename G1, typename G2>
   void
   ReceiveRestrictedMessageFromLeft(LEFT_FACTOR* l, const G1& leftPot, G2& msg, typename PrimalSolutionStorage::Element leftPrimal)
   { 
      for(INDEX i=0; i<edges_.size(); ++i) {
         assert(leftPrimal[i] == unknownState);
      }
      MakeLeftFactorUniform(leftPot, msg);
   }

   // for rounding based on unaries
   //template<typename RIGHT_FACTOR, typename G1, typename G2>
   //void
   //ReceiveRestrictedMessageFromRight(RIGHT_FACTOR* r, const G1& rightPot, G2& msg, typename PrimalSolutionStorage::Element rightPrimal)
   //{ 
   //   auto* mcf = r->GetMinCostFlowSolver();
   //   // if an arc is true, enforce it. If an arc is false, forbid it. Otherwise, get reduced costs
   //   for(INDEX l=0; l<edges_.size(); ++l) {
   //      const INDEX e = edges_[l];
   //      if(rightPrimal[edges_[l]] == unknownState) {
   //         msg[l] -= (mcf->GetReducedCost(e) );  // do zrobienia: this only works for left side models, not for right hand side ones!
   //      } else if(rightPrimal[e] == true) {
   //         msg[l] -= -std::numeric_limits<REAL>::max();
   //      } else {
   //         assert(rightPrimal[e] == false);
   //         msg[l] -= std::numeric_limits<REAL>::max();
   //      }
   //   }
   //}

   template<typename RIGHT_FACTOR, typename MSG_ARRAY, typename RIGHT_REPAM, typename ITERATOR>
   static void
   SendMessagesToLeft(const RIGHT_FACTOR& rightFactor, const RIGHT_REPAM& rightRepam, MSG_ARRAY msg_begin, const MSG_ARRAY msg_end, ITERATOR omegaIt)
   {
      auto* mcf = rightFactor.GetMinCostFlowSolver();
      assert(rightRepam.size() == rightFactor.size());
      for(INDEX i=0; i<rightFactor.size(); ++i) {
         assert(std::abs(rightRepam[i] - rightFactor[i]) <= eps);
      }
      REAL omega_sum = 0.0; //std::accumulate(omegaIt, omegaIt+msgs.size(),0.0); 
      for(auto it= msg_begin; it!=msg_end; ++it, ++omegaIt) {
        omega_sum += *omegaIt;
      }
           
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
      //for(INDEX i=0; i<msgs.size(); ++i, ++omegaIt) {
      for(; msg_begin != msg_end; ++msg_begin, ++omegaIt) {

         /*
         INDEX start_arc = msgs[i].GetMessageOp().start_arc_;
         INDEX no_arcs = msgs[i].GetMessageOp().no_arcs_;
         INDEX node = mcf->GetTailNodeId(start_arc);
         INDEX l_tmp=0;
         //std::cout << "node = " << i << "\n";
         for(auto it = mcf->begin(node); it != mcf->end(node); ++it) {
                 const REAL delta = omega_sum*1.0/REAL(COVERING_FACTOR)*((*it)->GetRCost());
                 auto a = mcf->N_arc(*it);
                 //std::cout << "a = " << a << "\n";
                 msgs[i][l_tmp] -= delta;
                 //mcf->UpdateCost( a, -delta);
                 ++l_tmp;
         }
         */
         //for(INDEX l=0; l<msgs[i].size(); ++l) {
         for(INDEX l=0; l<(*msg_begin).size(); ++l) {
            //msgs[i][l] -= omega_sum*1.0/REAL(COVERING_FACTOR)*(-mcf->GetReducedCost(start_arc + l) + mcf->GetCost(start_arc + l));
            /*
            const INDEX e = start_arc + l;
            assert(mcf->GetTailNodeId(e) == node);
            //std::cout << "e = " << e << "\n";
            assert(rightFactor.Is01Arc(e) || rightFactor.Is01Arc(mcf->GetReverseArcId(e))); // only {0,1}-flow accepted here!
            const REAL sign = mcf->GetCap(e) == 1 ? 1.0 : -1.0;
            assert(sign == 1.0);
            msgs[i][l] -= omega_sum*1.0/REAL(COVERING_FACTOR)*(sign*mcf->GetReducedCost(e) );  
            */
            const INDEX e = (*msg_begin).GetMessageOp().edges_[l];
            (*msg_begin)[l] -= omega_sum*1.0/REAL(COVERING_FACTOR)*(mcf->GetReducedCost(e) );  
         }
      }
   }

   template<bool ENABLE = (MPT == MessagePassingType::HUNGARIAN), typename LEFT_FACTOR, typename RIGHT_FACTOR>
   typename std::enable_if<ENABLE,void>::type
   ComputeLeftFromRightPrimal(LEFT_FACTOR& l, const RIGHT_FACTOR& r)
   {
      for(INDEX i=0; i<edges_.size(); ++i) {
         if(r.primal_[ edges_[i] ] > 0) {
            l.primal() = i;
         }
      }
   }

   template<bool ENABLE = (MPT == MessagePassingType::SRMP), typename LEFT_FACTOR, typename RIGHT_FACTOR>
   typename std::enable_if<ENABLE,void>::type
   ComputeRightFromLeftPrimal(const LEFT_FACTOR& l, RIGHT_FACTOR& r)
   {
      for(INDEX e=0; e<edges_.size(); ++e) {
         if(l.primal() == e) {
            r.primal_[edges_[e]] = 1;
         } else {
            r.primal_[edges_[e]] = 0;
         } 
      }
   }

   INDEX size() const { return edges_.size(); }

   template<typename G>
   void RepamLeft(G& leftRepamPot, const REAL msg, const INDEX dim) 
   {
      assert(dim < edges_.size());
      assert(leftRepamPot.size() == edges_.size());
      leftRepamPot[dim] += msg;
   }

   template<typename G>
   void RepamRight(G& rightRepamPot, const REAL msg, const INDEX dim) 
   {
      assert(dim < edges_.size());
      rightRepamPot[edges_[dim]] += msg;
   }

   template<class LEFT_FACTOR_TYPE,class RIGHT_FACTOR_TYPE>
   void CreateConstraints(LpInterfaceAdapter* lp,LEFT_FACTOR_TYPE* l, RIGHT_FACTOR_TYPE* r) const
   { 
      auto* mcf = r->GetMinCostFlowSolver();
      for(INDEX i=0; i<edges_.size(); ++i) {
         LinExpr lhs = lp->CreateLinExpr();
         LinExpr rhs = lp->CreateLinExpr();
         lhs += lp->GetLeftVariable(i);
         rhs += lp->GetRightVariable(edges_[i]);
         lp->addLinearEquality(lhs,rhs);
      }
   }


private:
   std::vector<INDEX> edges_;
};

} // end namespace LP_MP

#endif // LP_MP_MIN_COST_FLOW_FACTOR_CS2_HXX


