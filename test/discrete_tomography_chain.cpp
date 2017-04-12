#include "catch.hpp"
#include <vector>
#include "visitors/standard_visitor.hxx"
#include "solvers/discrete_tomography/discrete_tomography.h"

using namespace LP_MP;
using namespace std;

// to do: check correctness for pure sequential, pure recursive and mixed constructor on chains


template<typename SOLVER>
void add_projection_and_run(SOLVER& s, std::vector<INDEX> projection_var, const INDEX projection_sum) 
{
   auto& dt = s.template GetProblemConstructor<1>();

   std::vector<REAL> projectionCost(projection_sum+1);
   std::fill(projectionCost.begin(), projectionCost.end(), std::numeric_limits<REAL>::infinity());
   projectionCost.back() = 0.0;
   dt.AddProjection(projection_var, projectionCost); 

   s.Solve();
   if((projection_sum) % 10 != 0) {
      assert(std::abs(s.lower_bound() - 1.0) < eps);
      REQUIRE(std::abs(s.lower_bound() - 1.0) < eps);
      REQUIRE(std::abs(s.primal_cost() - 1.0) < eps);
   } else {
      assert(std::abs(s.lower_bound()) < eps);
      REQUIRE(std::abs(s.lower_bound()) < eps);
      REQUIRE(std::abs(s.primal_cost()) < eps);
   }
}



TEST_CASE("discrete tomography single chain", "[dt chain]") {

   char * options[5];
   options[0] = "";
   options[1] = "-i";
   options[2] = "";
   options[3] = "--maxIter";
   options[4] = "100";

   MpRoundingSolver<Solver<FMC_DT,LP_sat<LP>,StandardVisitor>> s(5,options);

   // add single Potts chain of length 10 with varying summation costs
   const INDEX noLabels = 3;
   matrix<REAL> PottsCost(3,3,0.0);
   for(INDEX x1=0; x1<noLabels; ++x1) {
      for(INDEX x2=0; x2<noLabels; ++x2) {
         if(x1 != x2) {
            PottsCost(x1,x2) = 1.0;
         }
      }
   }
   assert(PottsCost(0,0) == 0.0 && PottsCost(0,0) == 0.0 && PottsCost(0,0) == 0.0); 

   auto& mrf = s.template GetProblemConstructor<0>();
   for(INDEX i=0; i<10; ++i) {
      mrf.AddUnaryFactor(i,{0.0,0.0,0.0});
   }
   for(INDEX i=1; i<10; ++i) {
      mrf.AddPairwiseFactor(i-1,i,PottsCost);
   }
   auto& dt = s.template GetProblemConstructor<1>();
   dt.SetNumberOfLabels(3);

   std::vector<INDEX> projection_var {0,1,2,3,4,5,6,7,8,9};

   SECTION("sum = 0") { add_projection_and_run(s, projection_var, 0); }
   SECTION("sum = 1") { add_projection_and_run(s, projection_var, 1); }
   SECTION("sum = 2") { add_projection_and_run(s, projection_var, 2); }
   SECTION("sum = 3") { add_projection_and_run(s, projection_var, 3); }
   SECTION("sum = 4") { add_projection_and_run(s, projection_var, 4); }
   SECTION("sum = 5") { add_projection_and_run(s, projection_var, 5); }
   SECTION("sum = 6") { add_projection_and_run(s, projection_var, 6); }
   SECTION("sum = 7") { add_projection_and_run(s, projection_var, 7); }
   SECTION("sum = 8") { add_projection_and_run(s, projection_var, 8); }
   SECTION("sum = 9") { add_projection_and_run(s, projection_var, 9); }

   SECTION("sum = 10") { add_projection_and_run(s, projection_var, 10); }
   SECTION("sum = 11") { add_projection_and_run(s, projection_var, 11); }
   SECTION("sum = 12") { add_projection_and_run(s, projection_var, 12); }
   SECTION("sum = 13") { add_projection_and_run(s, projection_var, 13); }
   SECTION("sum = 14") { add_projection_and_run(s, projection_var, 14); }
   SECTION("sum = 15") { add_projection_and_run(s, projection_var, 15); }
   SECTION("sum = 16") { add_projection_and_run(s, projection_var, 16); }
   SECTION("sum = 17") { add_projection_and_run(s, projection_var, 17); }
   SECTION("sum = 18") { add_projection_and_run(s, projection_var, 18); }
   SECTION("sum = 19") { add_projection_and_run(s, projection_var, 19); }

   SECTION("sum = 20") { add_projection_and_run(s, projection_var, 20); }
}
