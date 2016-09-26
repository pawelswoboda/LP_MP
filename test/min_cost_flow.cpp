#include "catch.hpp"
#include "lib/MinCost/MinCost.h"
#include <random>

TEST_CASE( "ssp", "[min cost flow by successive shortest path]" ) {

   std::random_device rd;     // only used once to initialise (seed) engine
   std::mt19937 rng(rd());    // random-number engine used (Mersenne-Twister in this case)
   std::uniform_int_distribution<long> uni(0,10); // guaranteed unbiased

   SECTION("assignment problem") {
      const int num_nodes = 6;
      const int num_arcs = 9;
      MCF::SSP<int,int> mcf( num_nodes, num_arcs);

      for(int i=0; i<3; ++i) {
         for(int j=0; j<3; ++j) {
            mcf.AddEdge(i,3+j,0,1,uni(rng));
         }
      }
      for(int i=0; i<3; ++i) {
         mcf.AddNodeExcess(i,1);
         mcf.AddNodeExcess(3+i,-1);
      }

      mcf.SortArcs();

      REQUIRE(mcf.GetTailNodeId(0) == 0); REQUIRE(mcf.GetHeadNodeId(0) == 3);
      REQUIRE(mcf.GetTailNodeId(1) == 0); REQUIRE(mcf.GetHeadNodeId(1) == 4);
      REQUIRE(mcf.GetTailNodeId(2) == 0); REQUIRE(mcf.GetHeadNodeId(2) == 5);

      REQUIRE(mcf.GetTailNodeId(3) == 1); REQUIRE(mcf.GetHeadNodeId(3) == 3);
      REQUIRE(mcf.GetTailNodeId(4) == 1); REQUIRE(mcf.GetHeadNodeId(4) == 4);
      REQUIRE(mcf.GetTailNodeId(5) == 1); REQUIRE(mcf.GetHeadNodeId(5) == 5);

      REQUIRE(mcf.GetTailNodeId(6) == 2); REQUIRE(mcf.GetHeadNodeId(6) == 3);
      REQUIRE(mcf.GetTailNodeId(7) == 2); REQUIRE(mcf.GetHeadNodeId(7) == 4);
      REQUIRE(mcf.GetTailNodeId(8) == 2); REQUIRE(mcf.GetHeadNodeId(8) == 5);

      REQUIRE(mcf.GetCap(0) == 1); REQUIRE(mcf.GetCap(9) == 0);
      REQUIRE(mcf.GetCap(1) == 1); REQUIRE(mcf.GetCap(10) == 0);
      REQUIRE(mcf.GetCap(2) == 1); REQUIRE(mcf.GetCap(11) == 0);
      REQUIRE(mcf.GetCap(3) == 1); REQUIRE(mcf.GetCap(12) == 0);
      REQUIRE(mcf.GetCap(4) == 1); REQUIRE(mcf.GetCap(13) == 0);
      REQUIRE(mcf.GetCap(5) == 1); REQUIRE(mcf.GetCap(14) == 0);
      REQUIRE(mcf.GetCap(6) == 1); REQUIRE(mcf.GetCap(15) == 0);
      REQUIRE(mcf.GetCap(7) == 1); REQUIRE(mcf.GetCap(16) == 0);
      REQUIRE(mcf.GetCap(8) == 1); REQUIRE(mcf.GetCap(17) == 0);

      mcf.Solve();

      // exactly one unit was sent from each left node
      assert(mcf.GetFlow(0) + mcf.GetFlow(1) + mcf.GetFlow(2) == 1);
      assert(mcf.GetFlow(3) + mcf.GetFlow(4) + mcf.GetFlow(5) == 1);
      assert(mcf.GetFlow(6) + mcf.GetFlow(7) + mcf.GetFlow(8) == 1);
      REQUIRE(mcf.GetRCap(0) + mcf.GetRCap(1) + mcf.GetRCap(2) == 2);
      REQUIRE(mcf.GetRCap(3) + mcf.GetRCap(4) + mcf.GetRCap(5) == 2);
      REQUIRE(mcf.GetRCap(6) + mcf.GetRCap(7) + mcf.GetRCap(8) == 2);

      // exactly one unit was received in each right node
      REQUIRE(mcf.GetRCap(9+0) + mcf.GetRCap(9+1) + mcf.GetRCap(9+2) == 1);
      REQUIRE(mcf.GetRCap(9+3) + mcf.GetRCap(9+4) + mcf.GetRCap(9+5) == 1);
      REQUIRE(mcf.GetRCap(9+6) + mcf.GetRCap(9+7) + mcf.GetRCap(9+8) == 1);

      REQUIRE(mcf.NoArcs(0) == 3);
      REQUIRE(mcf.NoArcs(1) == 3);
      REQUIRE(mcf.NoArcs(2) == 3);
      REQUIRE(mcf.NoArcs(3) == 3);
      REQUIRE(mcf.NoArcs(4) == 3);
      REQUIRE(mcf.NoArcs(5) == 3);
   }
}
