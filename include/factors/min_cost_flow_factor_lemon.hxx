#ifndef LP_MP_MIN_COST_FLOW_FACTOR_LEMON_HXX
#define LP_MP_MIN_COST_FLOW_FACTOR_LEMON_HXX

#include "lemon/list_graph.h"
#include "lemon/network_simplex.h"
#include "lemon/cost_scaling.h"
#include "lemon/color.h"
#include "lemon/graph_to_eps.h"
#include <lemon/math.h>

#include <lemon/color.h>

// do zrobienia: this is needed for graph drawing
namespace lemon {

  const Color WHITE(1,1,1);

  const Color BLACK(0,0,0);
  const Color RED(1,0,0);
  const Color GREEN(0,1,0);
  const Color BLUE(0,0,1);
  const Color YELLOW(1,1,0);
  const Color MAGENTA(1,0,1);
  const Color CYAN(0,1,1);

  const Color GREY(0,0,0);
  const Color DARK_RED(.5,0,0);
  const Color DARK_GREEN(0,.5,0);
  const Color DARK_BLUE(0,0,.5);
  const Color DARK_YELLOW(.5,.5,0);
  const Color DARK_MAGENTA(.5,0,.5);
  const Color DARK_CYAN(0,.5,.5);

} //namespace lemon


namespace LP_MP {

// implement interface for MinCostFlowFactors

// build factor consisting of min cost network flow graph
// use as implementation the MinCost solver from the lemon project
// the reparametrization storage holds the min cost flow solver

// do zrobienia: a pointer to the min cost flow solver is held in the factor and in the ReparametrizationStorage, which is redundant. Possibly hold it just in one place via unique_ptr. Possibly do it via shared_ptr
class MinCostFlowFactorLemon {
public:
   using GraphType = lemon::ListDigraph;
   using MinCostFlowSolverType = lemon::NetworkSimplex<GraphType, SIGNED_INDEX, LONG_SIGNED_INDEX>;
   //using MinCostFlowSolverType = lemon::CostScaling<GraphType>;
   constexpr static REAL scalingFactor_ = 1000000000;
   // note: the scaling is somewhat brittle. Unless long unsigned int = 64 bit is used, monotone improvement often fails for eps in tolerance.hxx
   // do zrobienia: adaptive scaling?
   LONG_SIGNED_INDEX Scale(const REAL x) const { assert(std::abs(round(x)) < std::numeric_limits<LONG_SIGNED_INDEX>::max()/scalingFactor_ ); return std::round(scalingFactor_*x); }
   REAL Descale(const LONG_SIGNED_INDEX x) const { return REAL(x)/REAL(scalingFactor_); }
   struct Edge {
      const INDEX start_node, end_node, cap; // we assume that the reverse capacity is 0
      const REAL cost;
   };

   MinCostFlowFactorLemon(const std::vector<Edge>& edges, const std::vector<SIGNED_INDEX>& excess) 
   {
      noNodes_ = excess.size();
      noEdges_ = edges.size();

      assert(std::accumulate(excess.begin(), excess.end(),0) == 0);
      graph_ = new GraphType();
      graph_->clear();
      for(INDEX i=0; i<excess.size(); ++i) {
         nodes_.push_back(graph_->addNode());
      }
      for(INDEX e=0; e<edges.size(); ++e) {
         assert(edges[e].start_node < noNodes_ && edges[e].end_node < noNodes_ && edges[e].start_node != edges[e].end_node);
         arcs_.push_back(graph_->addArc(nodes_[edges[e].start_node], nodes_[edges[e].end_node]));
      }
      minCostFlow_ = new MinCostFlowSolverType(*graph_);
      // lower bounds will be zero, hence do not set them
      GraphType::ArcMap<SIGNED_INDEX> upper(*graph_);
      for(INDEX e=0; e<edges.size(); ++e) {
         upper[arcs_[e]] = edges[e].cap;
      }
      minCostFlow_->upperMap(upper);
      GraphType::ArcMap<LONG_SIGNED_INDEX> cost(*graph_);
      for(INDEX e=0; e<edges.size(); ++e) {
         cost[arcs_[e]] = Scale(edges[e].cost);
      }
      minCostFlow_->costMap(cost);
      GraphType::NodeMap<SIGNED_INDEX> supply(*graph_);
      for(INDEX i=0; i<excess.size(); ++i) {
         supply[nodes_[i]] = excess[i];
      }
      minCostFlow_->supplyMap(supply);
   }

   ~MinCostFlowFactorLemon()
   {
      // delete minCostFlow_; deleted in repam storage: not very clean!
      //delete graph_; // one has to delete the graph at the end! However add copy constructor
   }

   template<typename REPAM_ARRAY>
   REAL EvaluatePrimal(const REPAM_ARRAY& repam, const PrimalSolutionStorage::Element primal) const
   {
      REAL cost = 0.0;
      assert(repam.size() == arcs_.size());
      // assert that primal is indeed a flow
      for(INDEX e=0; e<repam.size(); e++) {
         cost += repam[e]*primal[e];
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
      GraphType::ArcMap<LONG_SIGNED_INDEX> cost(*graph_);
      for(INDEX e=0; e<arcs_.size(); ++e) {
         cost[arcs_[e]] = Scale(repam[e]);
      }
      auto ret = minCostFlow_->costMap(cost).run(); 
      assert(ret == MinCostFlowSolverType::OPTIMAL);
   }

   template<typename REPAM_ARRAY>
   void MaximizePotentialAndComputePrimal(const REPAM_ARRAY& repam, typename PrimalSolutionStorage::Element primal)
   { 
      MaximizePotential(repam);
      assert(repam.size() == arcs_.size());
      for(INDEX e=0; e<arcs_.size(); ++e) {
        primal[e] = minCostFlow_->flow(arcs_[e]);
      }
      //std::vector<INDEX> sol(arcs_.size());
      //for(INDEX e=0; e<arcs_.size(); ++e) {
      //   sol[e] = minCostFlow_->flow(arcs_[e]);
      //}
      //return std::make_pair(sol, Descale(minCostFlow_->totalCost()));
   }

   template<typename REPAM_ARRAY>
   REAL LowerBound(const REPAM_ARRAY& repamPot) const
   {
      assert(repamPot.size() == noEdges_);
      GraphType::ArcMap<LONG_SIGNED_INDEX> repam(*graph_);
      for(INDEX e=0; e<arcs_.size(); ++e) {
         repam[arcs_[e]] = Scale(repamPot[e]);
      }
      /*
      REAL min_pot = repamPot[0];
      REAL max_pot = repamPot[0];
      for(INDEX e=1;e<repamPot.size(); ++e) {
         min_pot = std::min(min_pot,repamPot[e]);
         max_pot = std::max(max_pot,repamPot[e]);
      }
      //std::cout << "minimum repam element " << min_pot << ", maximum repam element " << max_pot << "\n";
      */


      /*
      typedef lemon::dim2::Point<int> Point;
      GraphType::NodeMap<Point> coords(*graph_);
      for(INDEX i=0; i<nodes_.size()/2; ++i) {
         coords[nodes_[i]] = Point(3*i,1);
         coords[nodes_[nodes_.size() - i -1]] = Point(3*i,50);
      }
      GraphType::NodeMap<double> sizes(*graph_);
      for(INDEX i=0; i<nodes_.size(); ++i) {
         sizes[nodes_[i]] = 0.1;
      }
      lemon::Palette palette;
      GraphType::NodeMap<int> ncolors(*graph_);
      for(INDEX i=0; i<nodes_.size()/2; ++i) {
         ncolors[nodes_[i]] = 1;
         ncolors[nodes_[nodes_.size() - i -1]] = 2;
      }
      GraphType::ArcMap<int> acolors(*graph_);
      for(INDEX i=0; i<arcs_.size(); ++i) {
         if(minCostFlow_->flow(arcs_[i]) == 1) {
            acolors[arcs_[i]] = 1;
         } else {
            acolors[arcs_[i]] = 2;
         }
      }

      std::cout << "draw graph \n";
      graphToEps(*graph_,"test.eps").coords(coords).
         //absoluteNodeSizes().
         nodeSizes(sizes).
         arcWidthScale(.0001).
         nodeColors(lemon::composeMap(palette,ncolors)).
         arcColors(lemon::composeMap(palette,acolors)).
         run();
         */

      auto ret = minCostFlow_->costMap(repam).run();
      assert(ret == MinCostFlowSolverType::OPTIMAL);
      //std::cout << "kwaskwas: " << minCostFlow_->totalCost() << ", " << Descale(minCostFlow_->totalCost()) << "\n";
      
      return Descale(minCostFlow_->totalCost());
      // here we always assume that LowerBound is called with the current reparametrization
   }

   MinCostFlowSolverType* GetMinCostFlowSolver() const
   {
      return minCostFlow_;
   }

   const INDEX size() const 
   { 
      throw std::runtime_error("not supported"); return 0; 
      return noNodes_;
   }
   const REAL operator[](const INDEX i) const 
   {
      throw std::runtime_error("not supported"); return 0; 
      return 0.0;
   }

   /*
   template<typename REPAM_ARRAY, typename PERTURBATION_ARRAY>
   void MaximallyPerturb(const REPAM_ARRAY& repam, PERTURBATION_ARRAY& perturb)
   {
      // this function needs to be called after MaximizePotential, when messages need to be sent. Alternatively, when some message is received from this factor.
      // here we assume that all edges except the last one are under consideration.
      MinCostFlowSolverType minCostFlowRepamUpdate(*graph_);
      const SIGNED_INDEX max_cap = 100000;//minCostFlowRepamUpdate.INF; // note: INF is not good, as negative INF is handled incorrectly.
      GraphType::ArcMap<SIGNED_INDEX> lower(*graph_), upper(*graph_);
      GraphType::ArcMap<LONG_SIGNED_INDEX> cost(*graph_);
      for(INDEX e=0; e<repam.size(); ++e) { // do zrobienia: artifact. Factors should be designed such that repam is really the variables that are affected.
         const SIGNED_INDEX flow = minCostFlow_->flow(arcs_[e]);
         if(e < repam.size()-1) {
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
      //std::cout << "Solving reparametrization problem\n";
      auto ret = minCostFlowRepamUpdate.lowerMap(lower).upperMap(upper).costMap(cost).run();
      assert(ret == MinCostFlowSolverType::OPTIMAL);
      //std::cout << "done\n";
      // check now whether the flow reparametrization update does not violate the capacity constraints, i.e. whether max_cap is achieved somewhere
      for(INDEX e=0; e<arcs_.size(); ++e) {
         assert(minCostFlowRepamUpdate.flow(arcs_[e]) != max_cap && minCostFlowRepamUpdate.flow(arcs_[e]) != - max_cap);
      }

      for(INDEX e=0; e<repam.size(); ++e) {
         if(e < repam.size()-1) {
            const REAL reducedCost = Descale(Scale(repam[e]) + minCostFlowRepamUpdate.potential(graph_->source(arcs_[e])) - minCostFlowRepamUpdate.potential(graph_->target(arcs_[e])));
            if(minCostFlow_->flow(arcs_[e]) == 1) {
               perturb[e] = -reducedCost;
               assert(reducedCost < +eps);
               assert(reducedCost > -eps);
            } else if(minCostFlow_->flow(arcs_[e]) == 0) {
               perturb[e] = -reducedCost;
               assert(reducedCost > -eps);
            } else {
               assert(false); // flow is not 0/1
            }
         } else
            perturb[e] = 0.0;
      }
      // do zrobienia: construct test problem with perturbed cost and see whether solution is still optimal and just at edge
   }
   */

   // for the SendMessages update step, compute maximal cost change such that current labeling stays optimal
   // do zrobienia: is possibly an array with different values from which to compute the maximal perturbation
   template<class REPAM_ARRAY, typename ACTIVE_EDGES_ARRAY>
   std::vector<REAL> MaximallyPerturbCosts(const REPAM_ARRAY& repam, const ACTIVE_EDGES_ARRAY& active_edges) const
   {
      //std::cout << "maximum repam = " << *std::max_element(repam.begin(), repam.end()) << " minimum repam " << *std::min_element(repam.begin(), repam.end()) << "\n";
      assert(active_edges.size() == noEdges_);
      MinCostFlowSolverType minCostFlowRepamUpdate(*graph_);
      const SIGNED_INDEX max_cap = 100000;//minCostFlowRepamUpdate.INF; // note: INF is not good, as negative INF is handled incorrectly.
      //std::cout << "maximally perturb mcf\n";
      GraphType::ArcMap<SIGNED_INDEX> lower(*graph_), upper(*graph_);
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
               //lower[arcs_[e]] = -max_cap;
               upper[arcs_[e]] = +max_cap;
            } else {
               lower[arcs_[e]] = -max_cap;
               upper[arcs_[e]] = +max_cap;
            }
         }
         assert(lower[arcs_[e]] < upper[arcs_[e]]);
         //std::cout << repam[e] << ", ";
         cost[arcs_[e]] = Scale(repam[e]);
      }
      //std::cout << "\n";
      //std::cout << "Solving reparametrization problem\n";
      auto ret = minCostFlowRepamUpdate.lowerMap(lower).upperMap(upper).costMap(cost).run();
      assert(ret == MinCostFlowSolverType::OPTIMAL);
      if(ret != MinCostFlowSolverType::OPTIMAL) {
         throw std::runtime_error("Reparametrization problem could not be solved to optimality");
      }
      //std::cout << "done\n";
      // check now whether the flow reparametrization update does not violate the capacity constraints, i.e. whether max_cap is achieved somewhere
      for(INDEX e=0; e<arcs_.size(); ++e) {
         assert(minCostFlowRepamUpdate.flow(arcs_[e]) != max_cap && minCostFlowRepamUpdate.flow(arcs_[e]) != - max_cap);
      }

      std::vector<REAL> repam_cost;
      INDEX noActiveEdges = 0;
      for(INDEX e=0; e<active_edges.size(); ++e) {
         noActiveEdges += active_edges[e];
      }
      repam_cost.reserve(noActiveEdges);
      //std::cout << "max perturb\n";
      //repam_cost.reserve(std::count(active_edges.begin(),active_edges.end(),true));
      for(INDEX e=0; e<active_edges.size(); ++e) {
         if(active_edges[e] == true) {
            // do zrobienia: why should this be true?
            const REAL reducedCost = Descale(Scale(repam[e]) + minCostFlowRepamUpdate.potential(graph_->source(arcs_[e])) - minCostFlowRepamUpdate.potential(graph_->target(arcs_[e])));
            if(minCostFlow_->flow(arcs_[e]) == 1) {
               //std::cout << -reducedCost << ", ";
               repam_cost.push_back( -reducedCost );
               // these asserts are only valid if we make (-\infty,-1] and [1,\infty) the intervals for the flow in the reparametrization problem
               //assert(reducedCost < +eps);
               //assert(reducedCost > -eps);
            } else if(minCostFlow_->flow(arcs_[e]) == 0) {
               //std::cout << -reducedCost << ", ";
               repam_cost.push_back( -reducedCost );
               //assert(reducedCost > -eps);
            } else {
               assert(false); // flow is not 0/1
            }
         }
      }
      //std::cout << "\n";
      // do zrobienia: construct test problem with perturbed cost and see whether solution is still optimal and just at edge
      return repam_cost;
   }

private:
   // note: this is also given to the reparametrization storage and hence points to the reparametrized potential. Use shared_ptr?
   MinCostFlowSolverType* minCostFlow_;
   GraphType* graph_; // memory leak for this one!
   INDEX noNodes_;
   INDEX noEdges_;
   std::vector<GraphType::Node> nodes_;
   std::vector<GraphType::Arc> arcs_;

};

// Possibly can be replaced by RepamStorage for lemon mcf solver
template<typename FACTOR_CONTAINER>
class MinCostFlowReparametrizationStorageLemon {
public:
   using MinCostFlowSolverType = MinCostFlowFactorLemon::MinCostFlowSolverType;
   MinCostFlowReparametrizationStorageLemon(const MinCostFlowFactorLemon& mfc, const std::vector<REAL>& cost) 
      : minCostFlow_(mfc.GetMinCostFlowSolver()),
        repam_(cost)
   {
      assert(minCostFlow_ != nullptr);
   }

   ~MinCostFlowReparametrizationStorageLemon() {
      static_assert( std::is_same<typename FACTOR_CONTAINER::FactorType, MinCostFlowFactorLemon>::value, "");
      delete minCostFlow_;
   }
   
   //possibly not needed here.
   MinCostFlowSolverType* GetMinCostFlowSolver() const
   {
      return minCostFlow_;
   }

   const REAL operator[](const INDEX i) const 
   {
      return repam_[i];
   }

   REAL& operator[](const INDEX i) 
   {
      return repam_[i];
   }

   const INDEX size() const 
   {
      return repam_.size();
   }

private:
   // do zrobienia: what about shared_ptr?
   MinCostFlowSolverType* minCostFlow_; //possibly not needed here.
   std::vector<REAL> repam_;
   std::vector<REAL> maximumPerturbation_;
};

} // end namespace LP_MP

#endif // LP_MP_MIN_COST_FLOW_FACTOR_HXX

