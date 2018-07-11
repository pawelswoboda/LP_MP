#include "catch.hpp"
#include <vector>
#include "factors/simplex_factor.hxx"

using namespace LP_MP;

TEST_CASE( "unary simplex", "[unary simplex factor]" ) {
   std::vector<double> cost {0.1, 0.2, 0.05, 1};
   UnarySimplexFactor simplex(cost);
   SECTION( "lower bound" ) {
      REQUIRE( simplex.LowerBound() == 0.05 );
   }

   //std::vector<unsigned char> primal {false, false, true, false};
   SECTION( "primal evaluation" ) {
      //REQUIRE( simplex.EvaluatePrimal(primal.begin()) == 0.05);
      //primal[2] = false; primal[1] = true;
      //REQUIRE( simplex.EvaluatePrimal(primal.begin()) == 0.2);
   }

   SECTION( "primal computation" ) {
      //simplex.MaximizePotentialAndComputePrimal(primal.begin());
      //REQUIRE( primal[0] == false );  
      //REQUIRE( primal[1] == false );
      //REQUIRE( primal[2] == true );
      //REQUIRE( primal[3] == false );
   }

   /*
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
   */
}

TEST_CASE( "pairwise simplex", "[pairwise simplex factor]" ) {

   PairwiseSimplexFactor simplex(3,3);
   for(INDEX x1=0; x1<simplex.dim1(); ++x1) {
      for(INDEX x2=0; x2<simplex.dim2(); ++x2) {
         if(x1 != x2) {
            simplex.cost(x1,x2) = 0;
         } else {
            simplex.cost(x1,x2) = -REAL(x1)-1.0;
         }
      }
   }

   SECTION( "lower bound" ) {
      REQUIRE( simplex.LowerBound() == -3.0 );
   }
}
