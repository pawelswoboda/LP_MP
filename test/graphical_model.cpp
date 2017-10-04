#include "catch.hpp"
#include "solvers/graphical_model/graphical_model.h"
#include "LP_sat.hxx"
#include "visitors/standard_visitor.hxx"

#ifdef USE_GUROBI
#include "lp_interface/lp_gurobi.hxx"
#endif

// UAI test input. Note: not all unaries are present, hence zero unaries must be added.
std::string uai_test_input = 
R"(MARKOV
3
2 2 3
3
1 0
2 0 1
2 1 2

2
 0.436 0.564

4
 0.128 0.872
 0.920 0.080

6
 0.210 0.333 0.457
 0.811 0.000 0.189 
)";

std::vector<std::string> solver_options = {
   {"graphical model test"},
   {"--maxIter"}, {"100"},
   {"--timeout"}, {"60"}, // one minute
   {"--lowerBoundComputationInterval"}, {"1"},
   {"--primalComputationInterval"}, {"5"},
   {"--standardReparametrization"}, {"anisotropic"},
   {"--roundingReparametrization"}, {"anisotropic"},
   {"-v"}, {"0"},
   {"--inputFile"}, uai_test_input
};


TEST_CASE("gm", "[graphical model]") {
   auto sat_solver_options = solver_options;
   sat_solver_options[12] = std::string("damped_uniform");
   sat_solver_options.push_back("--satReductionMode");
   sat_solver_options.push_back("static");

   SECTION("MAP-estimation for a model given in uai format") {
      using FMC = FMC_SRMP;
      using VisitorType = StandardVisitor;
      using SolverType = MpRoundingSolver<Solver<FMC,LP,VisitorType>>;

      SolverType s(solver_options);
      UaiMrfInput::ParseString<SolverType,0>(uai_test_input, s); 
      s.Solve();

      REQUIRE(std::abs(s.lower_bound() - 0.564) < LP_MP::eps); // is this actually correct?
   }

   //std::vector<REAL> negPotts = {1, 0, 0 ,1};
   matrix<REAL> negPotts(2,2);
   negPotts(0,0) = 1.0; negPotts(1,0) = 0.0;
   negPotts(0,1) = 0.0; negPotts(1,1) = 1.0;
   //std::vector<REAL> posPotts = {0, 1, 1 ,0};
   matrix<REAL> posPotts(2,2);
   posPotts(0,0) = 0.0; posPotts(1,0) = 1.0;
   posPotts(0,1) = 1.0; posPotts(1,1) = 0.0;

   matrix<REAL> negPotts23(2,3);
   negPotts23(0,0) = 1.0; negPotts23(1,0) = 0.0; 
   negPotts23(0,1) = 0.0; negPotts23(1,1) = 1.0; 
   negPotts23(0,2) = 2.0; negPotts23(1,2) = 2.0; 
   //std::vector<REAL> negPotts23 = {1,0,2, 1,0,2};
   matrix<REAL> posPotts24(2,4);
   posPotts24(0,0) = 0.0; posPotts24(1,0) = 1.0;
   posPotts24(0,1) = 1.0; posPotts24(1,1) = 0.0;
   posPotts24(0,2) = 2.0; posPotts24(1,2) = 2.0;
   posPotts24(0,3) = 2.0; posPotts24(1,3) = 2.0;
   //std::vector<REAL> posPotts24 = {0,1,2,2, 1,0,2,2};
   matrix<REAL> posPotts34(3,4);
   posPotts34(0,0) = 0.0; posPotts34(1,0) = 1.0; posPotts34(2,0) = 2.0;
   posPotts34(0,1) = 1.0; posPotts34(1,1) = 0.0; posPotts34(2,1) = 2.0;
   posPotts34(0,2) = 2.0; posPotts34(1,2) = 2.0; posPotts34(2,2) = 2.0;
   posPotts34(0,3) = 2.0; posPotts34(1,3) = 2.0; posPotts34(2,3) = 2.0;
   //std::vector<REAL> posPotts34 = {0,1,2,2, 1,0,2,2, 2,2,2,2};

   SECTION("SAT rounding without duality gap") {
      using FMC = FMC_SRMP;
      using VisitorType = StandardVisitor;
      using SolverType = MpRoundingSolver<Solver<FMC,LP_sat<LP>,VisitorType>>;

      SolverType s(sat_solver_options);
      auto& mrf = s.template GetProblemConstructor<0>();

      mrf.AddUnaryFactor(std::vector<REAL>(2,0.0));
      mrf.AddUnaryFactor(std::vector<REAL>(2,0.0));
      mrf.AddUnaryFactor(std::vector<REAL>(2,0.0));
      mrf.AddUnaryFactor(std::vector<REAL>(2,0.0));
      mrf.AddUnaryFactor(std::vector<REAL>(2,0.0));

      // make a cycle of length 4 visiting each label once -> one negative Potts and three positive Potts
      mrf.AddPairwiseFactor(0,1,negPotts);
      mrf.AddPairwiseFactor(1,2,posPotts);
      mrf.AddPairwiseFactor(2,3,posPotts);
      mrf.AddPairwiseFactor(3,4,posPotts);

      s.Solve();
      REQUIRE(std::abs(s.lower_bound() - 0.0) <= eps);
      REQUIRE(std::abs(s.primal_cost() - 0.0) <= eps);
   }

   SECTION("SAT rounding with duality gap") {
      using FMC = FMC_SRMP;
      using VisitorType = StandardVisitor;
      using SolverType = MpRoundingSolver<Solver<FMC,LP_sat<LP>,VisitorType>>;

      SolverType s(sat_solver_options);
      auto& mrf = s.template GetProblemConstructor<0>();

      mrf.AddUnaryFactor(std::vector<REAL>(2,0.0));
      mrf.AddUnaryFactor(std::vector<REAL>(2,0.0));
      mrf.AddUnaryFactor(std::vector<REAL>(2,0.0));
      mrf.AddUnaryFactor(std::vector<REAL>(2,0.0));

      // make a cycle of length 4 visiting each label once -> one negative Potts and three positive Potts
      mrf.AddPairwiseFactor(0,1,negPotts);
      mrf.AddPairwiseFactor(1,2,posPotts);
      mrf.AddPairwiseFactor(2,3,posPotts);
      mrf.AddPairwiseFactor(0,3,posPotts);

      s.Solve();
      REQUIRE(std::abs(s.lower_bound() - 0.0) <= eps);
      REQUIRE(s.primal_cost() <= 3.0 + eps);
   }


   SECTION("Tightening") {
      using FMC = FMC_SRMP_T;
      using VisitorType = StandardTighteningVisitor;
      using SolverType = MpRoundingSolver<Solver<FMC,LP,VisitorType>>;

      auto tightening_solver_options = solver_options;
      tightening_solver_options.push_back("--tighten");
      tightening_solver_options.push_back("--tightenIteration");
      tightening_solver_options.push_back("5");
      tightening_solver_options.push_back("--tightenConstraintsMax");
      tightening_solver_options.push_back("1");
      tightening_solver_options.push_back("--tightenInterval");
      tightening_solver_options.push_back("10");
      tightening_solver_options.push_back("--tightenReparametrization");
      tightening_solver_options.push_back("damped_uniform");
      SolverType s(tightening_solver_options);
      auto& mrf = s.template GetProblemConstructor<0>();

      
      SECTION("binary triplet") {
         mrf.AddUnaryFactor(std::vector<REAL>(2,0.0));
         mrf.AddUnaryFactor(std::vector<REAL>(2,0.0));
         mrf.AddUnaryFactor(std::vector<REAL>(2,0.0));

         // make a cycle of length 3 visiting each label once -> one negative Potts and two positive Potts
         // make number of labels varying to check whether bounds work alright
         mrf.AddPairwiseFactor(0,1,negPotts);
         mrf.AddPairwiseFactor(0,2,posPotts);
         mrf.AddPairwiseFactor(1,2,posPotts);

         mrf.AddTighteningTriplet(0,1,2);
         s.Solve();
         REQUIRE(std::abs(s.lower_bound() - 1.0) <= eps);
      }

      SECTION("multi-label triplet") {
         mrf.AddUnaryFactor(std::vector<REAL>(2,0.0));
         mrf.AddUnaryFactor(std::vector<REAL>(3,0.0));
         mrf.AddUnaryFactor(std::vector<REAL>(4,0.0));

         // make a cycle of length 3 visiting each label once -> one negative Potts and two positive Potts
         // make number of labels varying to check whether bounds work alright
         mrf.AddPairwiseFactor(0,1,negPotts23);
         mrf.AddPairwiseFactor(0,2,posPotts24);
         mrf.AddPairwiseFactor(1,2,posPotts34);

         mrf.AddTighteningTriplet(0,1,2);
         s.Solve();
         REQUIRE(std::abs(s.lower_bound() - 1.0) <= eps);
      }

      SECTION("triplet search") {
         mrf.AddUnaryFactor(std::vector<REAL>(2,0.0));
         mrf.AddUnaryFactor(std::vector<REAL>(3,0.0));
         mrf.AddUnaryFactor(std::vector<REAL>(4,0.0));

         // make a cycle of length 3 visiting each label once -> one negative Potts and two positive Potts
         // make number of labels varying to check whether bounds work alright
         mrf.AddPairwiseFactor(0,1,negPotts23);
         mrf.AddPairwiseFactor(0,2,posPotts24);
         mrf.AddPairwiseFactor(1,2,posPotts34);

         s.Solve();
         REQUIRE(std::abs(s.lower_bound() - 1.0) <= eps);
      }

#ifdef USE_GUROBI
      SECTION("triplet with cplex") {
         auto cplex_options = solver_options;
         cplex_options.push_back("--onlyLp");
         cplex_options.push_back("--relax");
         VisitorLpSolver<SolverType,VisitorType,LpInterfaceGurobi> s(cplex_options);
         auto& mrf = s.template GetProblemConstructor<0>();

         mrf.AddUnaryFactor(std::vector<REAL>(2,0.0));
         mrf.AddUnaryFactor(std::vector<REAL>(3,0.0));
         mrf.AddUnaryFactor(std::vector<REAL>(4,0.0));

         // make a cycle of length 3 visiting each label once -> one negative Potts and two positive Potts
         // make number of labels varying to check whether bounds work alright
         mrf.AddPairwiseFactor(0,1,negPotts23);
         mrf.AddPairwiseFactor(0,2,posPotts24);
         mrf.AddPairwiseFactor(1,2,posPotts34);

         mrf.AddTighteningTriplet(0,1,2);

         s.Solve();
         REQUIRE(s.lower_bound() >= 1.0 - eps);
      }
#endif


      SECTION("cycle: k projection graph") { // k projection graph (individual labels) can find violated cycle
         mrf.AddUnaryFactor(std::vector<REAL>(2,0.0));
         mrf.AddUnaryFactor(std::vector<REAL>(2,0.0));
         mrf.AddUnaryFactor(std::vector<REAL>(2,0.0));
         mrf.AddUnaryFactor(std::vector<REAL>(2,0.0));

         // make a cycle of length 4 visiting each label once -> one negative Potts and three positive Potts
         mrf.AddPairwiseFactor(0,1,negPotts);
         mrf.AddPairwiseFactor(1,2,posPotts);
         mrf.AddPairwiseFactor(2,3,posPotts);
         mrf.AddPairwiseFactor(0,3,posPotts);

         s.Solve();
         REQUIRE(std::abs(s.lower_bound() - 1.0) <= eps);
      }

      // to do: add test problem where only create_expanded_projection_graph, but not create_k_projection_graph, can find a violated cycle.
      // to do: modified cycle inequalites do not work -> correct
      /*
      SECTION("cycle: expanded projection graph") {
         mrf.AddUnaryFactor(std::vector<REAL>(4,0.0));
         mrf.AddUnaryFactor(std::vector<REAL>(4,0.0));
         mrf.AddUnaryFactor(std::vector<REAL>(4,0.0));
         mrf.AddUnaryFactor(std::vector<REAL>(4,0.0));

         std::vector<REAL> Potts2 = {0,0,1,1, 0,0,1,1, 1,1,0,0, 1,1,0,0};
         std::vector<REAL> negPotts2;
         std::transform(Potts2.begin(), Potts2.end(), std::back_inserter(negPotts2), [](auto x) { return -x+1; });

         mrf.AddPairwiseFactor(0,1,negPotts2);
         mrf.AddPairwiseFactor(1,2,Potts2);
         mrf.AddPairwiseFactor(2,3,Potts2);
         mrf.AddPairwiseFactor(0,3,Potts2);

         mrf.AddTighteningTriplet(0,1,2);
         mrf.AddTighteningTriplet(0,1,3);
         mrf.AddTighteningTriplet(0,2,3);
         mrf.AddTighteningTriplet(1,2,3);

         s.Solve();
         assert(false);
         REQUIRE(s.lower_bound() == 1.0);
      }
      */

      SECTION("SAT rounding") {
         using SatSolverType = MpRoundingSolver<Solver<FMC,LP_sat<LP>,VisitorType>>;

         SatSolverType s_sat(solver_options);
         auto& mrf_sat = s_sat.template GetProblemConstructor<0>();

         mrf_sat.AddUnaryFactor(std::vector<REAL>(2,0.0));
         mrf_sat.AddUnaryFactor(std::vector<REAL>(2,0.0));
         mrf_sat.AddUnaryFactor(std::vector<REAL>(2,0.0));
         mrf_sat.AddUnaryFactor(std::vector<REAL>(2,0.0));

         // make a cycle of length 4 visiting each label once -> one negative Potts and three positive Potts
         mrf_sat.AddPairwiseFactor(0,1,negPotts);
         mrf_sat.AddPairwiseFactor(1,2,posPotts);
         mrf_sat.AddPairwiseFactor(2,3,posPotts);
         mrf_sat.AddPairwiseFactor(0,3,posPotts);

         mrf_sat.AddTighteningTriplet(0,1,2);
         mrf_sat.AddTighteningTriplet(0,1,3);
         mrf_sat.AddTighteningTriplet(0,2,3);
         mrf_sat.AddTighteningTriplet(1,2,3);

         s_sat.Solve();
         REQUIRE(std::abs(s_sat.lower_bound() - 1.0) <= eps);
         REQUIRE(std::abs(s_sat.primal_cost() - 1.0) <= eps);
      }
   }
}


