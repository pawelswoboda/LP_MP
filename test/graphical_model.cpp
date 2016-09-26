#include "catch.hpp"
#include "solvers/graphical_model/graphical_model.h"
#include "visitors/standard_visitor.hxx"

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
   {"--maxIter"}, {"10"},
   {"--timeout"}, {"60"}, // one minute
   {"--lowerBoundComputationInterval"}, {"1"},
   {"--standardReparametrization"}, {"anisotropic"},
   {"--roundingReparametrization"}, {"anisotropic"},
   {"--inputFile"}, uai_test_input
};


TEST_CASE("gm", "[graphical model]") {

   SECTION("MAP-estimation for a model given in uai format") {
      using FMC = FMC_SRMP;
      using VisitorType = StandardVisitor;
      using SolverType = Solver<FMC>;
      static auto Input = UaiMrfInput::ParseString<FMC,0>;

      VisitorSolver<SolverType,VisitorType> s(solver_options);
      s.ReadProblem(Input);
      s.Solve();

      REQUIRE(std::abs(s.lower_bound() - 0.564) < LP_MP::eps); // is this actually correct?
   }

   SECTION("Tightening") {
      using FMC = FMC_SRMP_T;
      using VisitorType = StandardTighteningVisitor;
      using SolverType = MpRoundingSolver<FMC>;

      auto tightening_solver_options = solver_options;
      tightening_solver_options.push_back("--tighten");
      tightening_solver_options.push_back("--tightenIteration");
      tightening_solver_options.push_back("5");
      tightening_solver_options.push_back("--tightenConstraintsMax");
      tightening_solver_options.push_back("1");
      tightening_solver_options.push_back("--tightenInterval");
      tightening_solver_options.push_back("10");
      tightening_solver_options.push_back("--tightenReparametrization");
      tightening_solver_options.push_back("uniform");
      VisitorSolver<SolverType,VisitorType> s(tightening_solver_options);
      auto& mrf = s.template GetProblemConstructor<0>();

      std::vector<REAL> negPotts = {1, 0, 0 ,1};
      std::vector<REAL> posPotts = {0, 1, 1 ,0};

      SECTION("triplet") {
         mrf.AddUnaryFactor(std::vector<REAL>(2,0.0));
         mrf.AddUnaryFactor(std::vector<REAL>(2,0.0));
         mrf.AddUnaryFactor(std::vector<REAL>(2,0.0));

         // make a cycle of length 3 visiting each label once -> one negative Potts and two positive Potts
         mrf.AddPairwiseFactor(0,1,negPotts);
         mrf.AddPairwiseFactor(0,2,posPotts);
         mrf.AddPairwiseFactor(1,2,posPotts);

         s.Solve();
         assert(s.lower_bound() == 1.0);
      }

      SECTION("cycle: k projection graph") { // k projection graph (individual labels) can find violated cycle
         mrf.AddUnaryFactor(std::vector<REAL>(2,0.0));
         mrf.AddUnaryFactor(std::vector<REAL>(2,0.0));
         mrf.AddUnaryFactor(std::vector<REAL>(2,0.0));
         mrf.AddUnaryFactor(std::vector<REAL>(2,0.0));

         // make a cycle of length 4 visiting each label once -> one negative Potts and two positive Potts
         mrf.AddPairwiseFactor(0,1,negPotts);
         mrf.AddPairwiseFactor(1,2,posPotts);
         mrf.AddPairwiseFactor(2,3,posPotts);
         mrf.AddPairwiseFactor(0,3,posPotts);

         s.Solve();
         assert(s.lower_bound() == 1.0);
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
         assert(s.lower_bound() == 1.0);
      }
      */
   }
}


