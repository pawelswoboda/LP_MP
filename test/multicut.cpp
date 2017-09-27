#include "catch.hpp"

#include "solvers/multicut/multicut_factors_messages.hxx"
#include "solvers/multicut/lifted_multicut_factors_messages.hxx"

using namespace LP_MP;

TEST_CASE( "multicut", "[multicut]" ) {

   multicut_edge_factor unary;
   SECTION("edge factor") {
      unary[0] = 1.0;
      REQUIRE(unary.size() == 1);
      REQUIRE(unary.LowerBound() == 0);
      unary[0] = -1.0;
      REQUIRE(unary.LowerBound() == -1);
   }

   SECTION("triangle factor") {
      multicut_triplet_factor triangle;
      REQUIRE(triangle.size() == 4);

      triangle[0] = 1.0;
      triangle[1] = 2.0;
      triangle[2] = 3.3;
      triangle[3] = 1.5;
      REQUIRE(triangle.LowerBound() == 0.0);

      triangle[1] = -0.5;
      triangle[2] = -0.3;
      REQUIRE(triangle.LowerBound() == -0.5);
   }
}

TEST_CASE("lifted multicut", "[lifted multicut]") {
   SECTION("cut factor") {
      //LiftedMulticutCutFactor cutFactor(3);
      //cutFactor.IncreaseLifted();
      //REQUIRE(cutFactor.size() == 4);

   }
}
