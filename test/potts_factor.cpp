#include "catch.hpp"
#include <array>
#include "factors/simplex_factor.hxx"

using namespace LP_MP;
using namespace std;

void set_potts(PairwiseSimplexFactor& p, const REAL diff)
{
   for(INDEX x1=0; x1<p.dim1(); ++x1) { 
      for(INDEX x2=0; x2<p.dim2(); ++x2) { 
         if(x1 == x2) {
            p.cost(x1,x2) = 0.0;
         } else {
            p.cost(x1,x2) = diff;
         }
      }
   }
}

void test_factor_equal(pairwise_potts_factor& p1, PairwiseSimplexFactor p2)
{
   REQUIRE( p1.LowerBound() == p2.LowerBound() ); 

   std::array<REAL,3> msg;
   std::array<REAL,3> msg2;

   p1.min_marginal_1(msg); p2.min_marginal_1(msg2);
   REQUIRE(msg == msg2);

   p1.min_marginal_2(msg); p2.min_marginal_2(msg2);
   REQUIRE(msg == msg2);
}

TEST_CASE( "pairwise Potts", "[Potts factor]" ) {
   pairwise_potts_factor potts(3,1.0);
   PairwiseSimplexFactor potts2(3,3);
   set_potts(potts2,1.0);


   SECTION( "lower bound, positive coupling" ) {
      test_factor_equal(potts, potts2);
   }

   SECTION( "lower bound, negative coupling" ) {
      potts.diff_cost() = -1.0;
      set_potts(potts2,-1.0);
      test_factor_equal(potts, potts2);
   }

   potts.msg1(0) = -0.1; potts2.msg1(0) = -0.1;
   potts.msg1(1) = 0.5;  potts2.msg1(1) = 0.5;
   potts.msg1(2) = 0.8;  potts2.msg1(2) = 0.8;

   potts.msg2(0) = 1.5;  potts2.msg2(0) = 1.5;
   potts.msg2(1) = 1.0;  potts2.msg2(1) = 1.0;

   SECTION( "lower bound, positive coupling, with messages" ) {
      test_factor_equal(potts, potts2);
   }

   SECTION( "lower bound, negative coupling, with messages" ) {
      potts.diff_cost() = -1.0;
      set_potts(potts2,-1.0);
      test_factor_equal(potts, potts2);
   }

   SECTION( "primal evaluation" ) {
   }


}
