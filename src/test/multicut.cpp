#include "catch.hpp"

#include "solvers/multicut/multicut_unary_factor.hxx"
#include "solvers/multicut/multicut_triplet_factor.hxx"
#include "solvers/multicut/multicut_odd_wheel.hxx"
#include "solvers/multicut/lifted_multicut_factors_messages.hxx"

using namespace LP_MP;
using namespace std;

TEST_CASE( "multicut", "[multicut]" ) {

   MulticutUnaryFactor unary(1);
   SECTION("edge factor") {
      REQUIRE(unary.size() == 1);
      vector<double> cost(1,1.0);
      vector<unsigned char> primal(1,unknownState);
      REQUIRE(unary.LowerBound(cost) == 0);
      unary.MaximizePotentialAndComputePrimal(cost,primal.begin());
      REQUIRE(primal[0] == 0);
      cost[0] = -1;
      primal[0] = unknownState;
      REQUIRE(unary.LowerBound(cost) == -1);
      unary.MaximizePotentialAndComputePrimal(cost,primal.begin());
      REQUIRE(primal[0] == 1);

   }

   SECTION("triangle factor") {
      MulticutTripletFactor triangle;
      REQUIRE(triangle.size() == 4);
      REQUIRE(triangle.PrimalSize() == triangle.size() + 1);
      vector<double> cost({1.0,2,3.3,1.5});
      REQUIRE(triangle.LowerBound(cost) == 0.0);
      vector<unsigned char> primal({false,false,false,false,true});
      REQUIRE(triangle.EvaluatePrimal(cost,primal.begin()) == 0.0);
      cost[1] = -0.5;
      cost[2] = -0.3;
      REQUIRE(triangle.LowerBound(cost) == -0.5);
      REQUIRE(triangle.EvaluatePrimal(cost,primal.begin()) == 0.0);
      primal[4] = false;
      primal[2] = true;
      REQUIRE(triangle.EvaluatePrimal(cost,primal.begin()) == -0.3);
      primal[2] = false;
      REQUIRE(triangle.EvaluatePrimal(cost,primal.begin()) == std::numeric_limits<double>::infinity());

      primal[3] = unknownState;
      triangle.PropagatePrimal(primal.begin());
      REQUIRE(primal[3] == true);

      primal[1] = unknownState;
      primal[3] = unknownState;
      primal[2] = true;
      triangle.PropagatePrimal(primal.begin());
      REQUIRE(primal[3] == false);
      REQUIRE(primal[1] == false);

      primal = {unknownState, unknownState,false,false,false};
      triangle.PropagatePrimal(primal.begin());
      REQUIRE(primal[0] == unknownState);
      REQUIRE(primal[1] == unknownState);
      REQUIRE(primal[2] == false);
      REQUIRE(primal[3] == false);
      REQUIRE(primal[4] == false);
   }

   SECTION("spoke factor") {
      MulticutTripletPlusSpokeFactor spoke;
      vector<double> cost({1.0,2,3.3,1.5,
                          -1.0,-1.0,-0.9,0});
      REQUIRE(spoke.LowerBound(cost) == -1.0);
      for_each(cost.begin(), cost.end(), [](double& x) { x = std::abs(x); });
      REQUIRE(spoke.LowerBound(cost) == 0.0);
   }

}

TEST_CASE("lifted multicut", "[lifted multicut]") {
   SECTION("cut factor") {
      LiftedMulticutCutFactor cutFactor(3);
      cutFactor.IncreaseLifted();
      REQUIRE(cutFactor.size() == 4);

   }
}
