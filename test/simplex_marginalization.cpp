#include "catch.hpp"
#include <vector>
#include "factors/simplex_factor.hxx"
#include "messages/message_loops.hxx"
#include "messages/simplex_marginalization_message.hxx"

using namespace LP_MP;
using namespace std;

TEST_CASE( "unary/pairwise message between simplex factors", "[unary/pairwise message]" ) {
   vector<double> costUnary {0.1, 0.2, 0.05, 1};
   vector<double> costPairwise {  0.1,    0.2,  0.05,
                                  0.3,  0.001,  0.2 ,
                                 -0.3, -0.001, -0.2 ,
                                  0.3,  0.001,  0.2 };
   SimplexFactor<> simplexUnary(costUnary);
   SimplexFactor<> simplexPairwise(costPairwise);

   using UnaryLoopType = UnaryLoop<>;
   using LeftLoopType = PairwiseLoop<0>;
   using RightLoopType = PairwiseLoop<1>;
   UnaryLoopType unaryLoop(4);
   RightLoopType rightLoop({3,4});
   LeftLoopType leftLoop({3,4});
   SimplexMarginalizationMessage<UnaryLoopType,RightLoopType,true,false> rightMessage(unaryLoop,rightLoop);
   SimplexMarginalizationMessage<UnaryLoopType,LeftLoopType,true,false> leftMessage(unaryLoop,leftLoop);

   SECTION( "marginalize pairwise right" ) {
      vector<double> marg(4,0.0);
      rightMessage.ReceiveMessageFromRight(&simplexPairwise, costPairwise, marg);
      REQUIRE(marg[0] == -0.05);
      REQUIRE(marg[1] == -0.001);
      REQUIRE(marg[2] ==  0.3);
      REQUIRE(marg[3] == -0.001);
   }

   SECTION( "marginalize pairwise left" ) {
      vector<double> marg(3,0.0);
      leftMessage.ReceiveMessageFromRight(&simplexPairwise, costPairwise, marg);
      REQUIRE(marg[0] == 0.3);
      REQUIRE(marg[1] == 0.001);
      REQUIRE(marg[2] == 0.2);
   }

   // do zrobienia: test ComputeLeftFromRightPrimal and reverse, repamLeft,repamRight, Send{Left|Right}Message and also interchange unary and right/left loop

}


// do zrobienia: triplet to pairwise marginalization
TEST_CASE("pairwise/triplet message between simplex factors", "[pairwise/triplet message]") {

}
