#include "catch.hpp"
#include <vector>
#include "factors/simplex_factor.hxx"

using namespace LP_MP;
using namespace std;

TEST_CASE( "simplex", "[simplex factor]" ) {
   std::vector<double> cost {0.1, 0.2, 0.05, 1};
   SimplexFactor<> simplex(cost);
   SECTION( "lower bound" ) {
      REQUIRE( simplex.LowerBound(cost) == 0.05 );
   }

   vector<unsigned char> primal {false, false, true, false};
   SECTION( "primal evaluation" ) {
      REQUIRE( simplex.EvaluatePrimal(cost, primal.begin()) == 0.05);
      primal[2] = false; primal[1] = true;
      REQUIRE( simplex.EvaluatePrimal(cost, primal.begin()) == 0.2);
   }

   SECTION( "primal computation" ) {
      simplex.MaximizePotentialAndComputePrimal(cost, primal.begin());
      REQUIRE( primal[0] == false );  
      REQUIRE( primal[1] == false );
      REQUIRE( primal[2] == true );
      REQUIRE( primal[3] == false );
   }

   SECTION( "primal propagation" ) {
      primal = {false,unknownState,false,false};
      simplex.PropagatePrimal(primal.begin());
      REQUIRE( primal[0] == false );  
      REQUIRE( primal[1] == true );
      REQUIRE( primal[2] == false );
      REQUIRE( primal[3] == false );
   }

   SECTION( "primal propagation" ) {
      primal = {unknownState,true,unknownState,unknownState};
      simplex.PropagatePrimal(primal.begin());
      REQUIRE( primal[0] == false );  
      REQUIRE( primal[1] == true );
      REQUIRE( primal[2] == false );
      REQUIRE( primal[3] == false );
   }
}
