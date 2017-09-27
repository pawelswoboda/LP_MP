#include "catch.hpp"
#include <vector>
#include "factors/simplex_factor.hxx"
#include "messages/simplex_marginalization_message.hxx"

using namespace LP_MP;

TEST_CASE( "unary/pairwise message", "[unary/pairwise message between simplex factors]" ) {
   std::vector<double> costUnary {0.1, 0.2, 0.05, 1};
   std::vector<double> costPairwise {  0.1,    0.2,  0.05,
                                  0.3,  0.001,  0.2 ,
                                 -0.3, -0.001, -0.2 ,
                                  0.3,  0.001,  0.2 };
   UnarySimplexFactor simplexUnary(costUnary);
   PairwiseSimplexFactor simplexPairwise(4,3);
   for(INDEX i=0; i<4; ++i) {
      for(INDEX j=0; j<3; ++j) {
         simplexPairwise.cost(i,j) = costPairwise[i*3 + j];
      }
   }

   UnaryPairwiseMessage<Chirality::left,true> leftMessage(4,3);
   UnaryPairwiseMessage<Chirality::right,true> rightMessage(4,3);

   // must add operators -= and += to vector to support the below things
   SECTION( "marginalize pairwise right" ) {
      vector<REAL> marg(4,0.0);
      leftMessage.send_message_to_left(simplexPairwise, marg);
      REQUIRE(marg[0] == -0.05);
      REQUIRE(marg[1] == -0.001);
      REQUIRE(marg[2] ==  0.3);
      REQUIRE(marg[3] == -0.001);
   }

   SECTION( "marginalize pairwise left" ) {
      vector<REAL> marg(3,0.0);
      rightMessage.send_message_to_left(simplexPairwise, marg);
      REQUIRE(marg[0] == 0.3);
      REQUIRE(marg[1] == 0.001);
      REQUIRE(marg[2] == 0.2);
   }

   // do zrobienia: test ComputeLeftFromRightPrimal and reverse, repamLeft,repamRight, Send{Left|Right}Message and also interchange unary and right/left loop

}

/*
TEST_CASE("pairwise/triplet message", "[pairwise/triplet message between simplex factors]") {
   PairwiseSimplexFactor simplexPairwise(4,3);
   for(INDEX i=0; i<3; ++i) {
      for(INDEX j=0; j<3; ++j) {
         simplexPairwise(i,j) = costPairwise[i*3 + j];
      }
   }
   SimpleTighteningTernarySimplexFactor simplexTriplet(3,3,3);

   PairwiseTripletMessage12<MessageSendingType::SRMP> message(3,3,3);

   SECTION( "marginalize triplet" ) {
      matrix marg(3,3,0.0);
      message.ReceiveMessageFromRight(simplexTriplet, marg);
   }

   SECTION( "marginalize pairwise" ) {
      matrix marg(3,3,0.0);
      message.ReceiveMessageFromLeft(simplexPairwise, marg);
   }
}
*/
