#include "catch.hpp"
#include <vector>
#include "visitors/standard_visitor.hxx"
#include "solvers/discrete_tomography/discrete_tomography.h"

using namespace LP_MP;
using namespace std;

TEST_CASE( "discrete tomography", "[dt]" ) {

   std::vector<std::string> options = {
      {"discrete tomography test"},
      {"--inputFile"}, {""}, // note: we dot not have an input file, but argument is mandatory
      {"--maxIter"}, {"10"},
      {"--lowerBoundComputationInterval"}, {"1"}
   };
   
   VisitorSolver<Solver<FMC_DT>,StandardVisitor> s(options);
   auto& mrf = s.GetProblemConstructor<0>();
   auto& dt = s.GetProblemConstructor<1>();
   dt.SetNumberOfLabels(3);


   SECTION( "tree" ) { // build a tree for a single projection and optimize. Check, whether global optimum could be achieved.
      // first build chain with four nodes
      mrf.AddUnaryFactor(0,std::vector<REAL>{0.0,1.0,0.5});
      mrf.AddUnaryFactor(1,std::vector<REAL>{0.0,1.0,0.5});
      mrf.AddUnaryFactor(2,std::vector<REAL>{0.0,1.0,0.5});
      mrf.AddUnaryFactor(3,std::vector<REAL>{0.0,1.0,0.5});

      // now link with pairwise Potts factors
      std::vector<REAL> PottsCost = {0.0,1.0,1.0,
                                     1.0,0.0,1.0,
                                     1.0,1.0,0.0};
      mrf.AddPairwiseFactor(0,1,PottsCost);
      mrf.AddPairwiseFactor(1,2,PottsCost);
      mrf.AddPairwiseFactor(2,3,PottsCost);

      std::vector<INDEX> projectionVar {0,1,2,3};
      std::vector<REAL> projectionCost {10,10,0,100};

      dt.AddProjection(projectionVar,projectionCost);
      s.Solve();
   }

}

