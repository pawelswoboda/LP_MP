#include "catch.hpp"
#include <vector>
#include "visitors/standard_visitor.hxx"
#include "lp_interface/lp_gurobi.hxx"
#include "solvers/discrete_tomography/discrete_tomography.h"

using namespace LP_MP;
using namespace std;

template<typename SOLVER>
void RunTestModel(SOLVER& s)
{
   const INDEX noLabels = 3;
   std::vector<REAL> PottsCost = {0.0,1.0,1.0,
                                  1.0,0.0,1.0,
                                  1.0,1.0,0.0};

   std::vector<INDEX> projectionVar {0,1,2,3};
   std::vector<REAL> projectionCost {std::numeric_limits<REAL>::infinity(),std::numeric_limits<REAL>::infinity(),0};

   auto& mrf = s.template GetProblemConstructor<0>();
   auto& dt = s.template GetProblemConstructor<1>();
   dt.SetNumberOfLabels(3);

   // first build chain with four nodes
   mrf.AddUnaryFactor(0,std::vector<REAL>{0.0,0.0,0.0});
   mrf.AddUnaryFactor(1,std::vector<REAL>{0.0,0.0,0.0});
   mrf.AddUnaryFactor(2,std::vector<REAL>{0.0,0.0,0.0});
   mrf.AddUnaryFactor(3,std::vector<REAL>{0.0,0.0,0.0});

   // now link with pairwise Potts factors
   mrf.AddPairwiseFactor(0,1,PottsCost);
   mrf.AddPairwiseFactor(1,2,PottsCost);
   mrf.AddPairwiseFactor(2,3,PottsCost);

   dt.AddProjection(projectionVar,projectionCost);
   s.Solve();

   std::vector<REAL> vars(4,0.0);
   REAL counting = 0.0;

   // do zrobienia: get solution and see whether it is feasible
}


TEST_CASE( "discrete tomography LP", "[dt lp]" ) {

   std::string tmp_file_name = std::tmpnam(nullptr);
   std::cout << "temporary file name: " << tmp_file_name << '\n';

   std::vector<std::string> options = {
      {"discrete tomography test"},
      //{"--lpFile"},{"model.lp"},
      {"-i"}, {""}, // note: we dot not have an input file, but argument is mandatory
      //{"--maxIter"}, {"10"},
      //{"--lowerBoundComputationInterval"}, {"1"}
   };
   
   SECTION( "LP interface tight relaxation" ) { 
      VisitorSolver<LpSolver<Solver<FMC_DT>,LpInterfaceGurobi>,StandardVisitor> s(options);
      RunTestModel(s);
   }

   SECTION( "LP interface naive relaxation" ) {
      VisitorSolver<LpSolver<Solver<FMC_DT_NAIVE>,LpInterfaceGurobi>,StandardVisitor> s_naive(options);
      RunTestModel(s_naive);
   }
}
