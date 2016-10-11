#include "catch.hpp"
#include "lib/MinCost/MinCost.h"
#include <random>

TEST_CASE( "ssp", "[min cost flow by successive shortest path]" ) {

        SECTION( "test problem" ) {
                const int num_nodes = 6;
                const int num_arcs = 8;
                MCF::SSP<int,int> mcf( num_nodes, num_arcs);

                mcf.AddEdge( 0, 1, 0, 4, 1);
                mcf.AddEdge( 0, 2, 0, 8, 5);
                mcf.AddEdge( 1, 2, 0, 5, 0);
                mcf.AddEdge( 2, 4, 0, 10, 1);
                mcf.AddEdge( 3, 1, 0, 8, 1);
                mcf.AddEdge( 3, 5, 0, 8, 1);
                mcf.AddEdge( 4, 3, 0, 8, 0);
                mcf.AddEdge( 4, 5, 0, 8, 9);
                mcf.AddNodeExcess( 0, 10);
                mcf.AddNodeExcess( 5, -10);

                mcf.SortArcs();
                long obj = mcf.Solve();
                REQUIRE(obj == 70);
                // after solving, arcs are reordered lexicographically. Note that reverse arc is implicitly added by AddEdge as well.
                for(int e=0; e<16; ++e) {
                        std::cout << mcf.GetTailNodeId(e) << " -> " << mcf.GetHeadNodeId(e) << "; cost = " << mcf.GetCost(e) << "\n";
                }
                //REQUIRE(mcf.compute_objective_cost() == 70);

                // correct arc ordering
                REQUIRE(mcf.GetTailNodeId(0) == 0); REQUIRE(mcf.GetHeadNodeId(0) == 1);
                REQUIRE(mcf.GetTailNodeId(1) == 0); REQUIRE(mcf.GetHeadNodeId(1) == 2);
                REQUIRE(mcf.GetTailNodeId(2) == 1); REQUIRE(mcf.GetHeadNodeId(2) == 0);
                REQUIRE(mcf.GetTailNodeId(3) == 1); REQUIRE(mcf.GetHeadNodeId(3) == 2);
                REQUIRE(mcf.GetTailNodeId(4) == 1); REQUIRE(mcf.GetHeadNodeId(4) == 3);
                REQUIRE(mcf.GetTailNodeId(5) == 2); REQUIRE(mcf.GetHeadNodeId(5) == 0);
                REQUIRE(mcf.GetTailNodeId(6) == 2); REQUIRE(mcf.GetHeadNodeId(6) == 1);
                REQUIRE(mcf.GetTailNodeId(7) == 2); REQUIRE(mcf.GetHeadNodeId(7) == 4);
                REQUIRE(mcf.GetTailNodeId(8) == 3); REQUIRE(mcf.GetHeadNodeId(8) == 1);
                REQUIRE(mcf.GetTailNodeId(9) == 3); REQUIRE(mcf.GetHeadNodeId(9) == 4);
                REQUIRE(mcf.GetTailNodeId(10) == 3); REQUIRE(mcf.GetHeadNodeId(10) == 5);
                REQUIRE(mcf.GetTailNodeId(11) == 4); REQUIRE(mcf.GetHeadNodeId(11) == 2);
                REQUIRE(mcf.GetTailNodeId(12) == 4); REQUIRE(mcf.GetHeadNodeId(12) == 3);
                REQUIRE(mcf.GetTailNodeId(13) == 4); REQUIRE(mcf.GetHeadNodeId(13) == 5);
                REQUIRE(mcf.GetTailNodeId(14) == 5); REQUIRE(mcf.GetHeadNodeId(14) == 3);
                REQUIRE(mcf.GetTailNodeId(15) == 5); REQUIRE(mcf.GetHeadNodeId(15) == 4);

                REQUIRE(mcf.StartingArc(0) == 0); REQUIRE(mcf.NoArcs(0) == 2);
                REQUIRE(mcf.StartingArc(1) == 2); REQUIRE(mcf.NoArcs(1) == 3);
                REQUIRE(mcf.StartingArc(2) == 5); REQUIRE(mcf.NoArcs(2) == 3);
                REQUIRE(mcf.StartingArc(3) == 8); REQUIRE(mcf.NoArcs(3) == 3);
                REQUIRE(mcf.StartingArc(4) == 11); REQUIRE(mcf.NoArcs(4) == 3);
                REQUIRE(mcf.StartingArc(5) == 14); REQUIRE(mcf.NoArcs(5) == 2);

                // correct flow values
                // here I must know the original capacities
                //REQUIRE(mcf.GetFlow(0) == 4);
                //REQUIRE(mcf.GetFlow(1) == 6);
                //REQUIRE(mcf.GetFlow(2) == 4);
                //REQUIRE(mcf.GetFlow(3) == 10);
                //REQUIRE(mcf.GetFlow(4) == 8);
                //REQUIRE(mcf.GetFlow(5) == 0);
                //REQUIRE(mcf.GetFlow(6) == 8);
                //REQUIRE(mcf.GetFlow(7) == 2);

                REQUIRE(4-mcf.GetRCap(0) == 4);
                REQUIRE(8-mcf.GetRCap(1) == 6);
                REQUIRE(5-mcf.GetRCap(3) == 4);
                REQUIRE(10-mcf.GetRCap(7) == 10);
                REQUIRE(8-mcf.GetRCap(8) == 0);
                REQUIRE(8-mcf.GetRCap(10) == 8);
                REQUIRE(8-mcf.GetRCap(12) == 8);
                REQUIRE(8-mcf.GetRCap(13) == 2);

                REQUIRE(mcf.GetFlow(0) == 4);
                REQUIRE(mcf.GetFlow(1) == 6);
                REQUIRE(mcf.GetFlow(3) == 4);
                REQUIRE(mcf.GetFlow(7) == 10);
                REQUIRE(mcf.GetFlow(8) == 0);
                REQUIRE(mcf.GetFlow(10) == 8);
                REQUIRE(mcf.GetFlow(12) == 8);
                REQUIRE(mcf.GetFlow(13) == 2);

                // complementary slackness
                REQUIRE(mcf.GetReducedCost(0) <= 0);
                REQUIRE(mcf.GetReducedCost(1) == 0);
                REQUIRE(mcf.GetReducedCost(3) == 0);
                REQUIRE(mcf.GetReducedCost(7) <= 0);
                REQUIRE(mcf.GetReducedCost(8) >= 0);
                REQUIRE(mcf.GetReducedCost(10) <= 0);
                REQUIRE(mcf.GetReducedCost(12) <= 0);
                REQUIRE(mcf.GetReducedCost(13) == 0);
        }





        SECTION("assignment problem") {
                std::random_device rd;     // only used once to initialise (seed) engine
                std::mt19937 rng(rd());    // random-number engine used (Mersenne-Twister in this case)
                std::uniform_int_distribution<long> uni(0,10); // guaranteed unbiased

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


        SECTION( "dynamic problem" ) {
                std::random_device rd;     // only used once to initialise (seed) engine
                std::mt19937 rng(rd());    // random-number engine used (Mersenne-Twister in this case)
                std::uniform_int_distribution<long> uni(0,10000); // guaranteed unbiased
                // build assignment problems and change cost after having found a solution

                SECTION( "cost update" ) {
                        for(int run =0; run<100; ++run) {
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
                                long orig_cost = mcf.Solve();

                                // read out solutions
                                std::vector<int> flow(num_arcs);
                                for(int i=0; i<num_arcs; ++i) {
                                        flow[i] = mcf.GetFlow(i);
                                }

                                REQUIRE(flow[0] + flow[1] + flow[2] == 1);
                                for(int i=0; i<3; ++i) {
                                        if(flow[i] == 1) { // increase cost by one and check whether cost of optimal solution is also increased by one
                                                mcf.UpdateCost(i, 1);
                                        }
                                }
                                long new_cost = mcf.Solve();
                                REQUIRE(orig_cost <= new_cost);
                                REQUIRE(new_cost <= orig_cost + 1);
                        }
                }

                SECTION( "capacity update" ) {
                        for(int run =0; run<100; ++run) {
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
                                long orig_cost = mcf.Solve();
                                std::vector<int> flow(num_arcs);
                                for(int i=0; i<num_arcs; ++i) {
                                        flow[i] = mcf.GetFlow(i);
                                }


                                // double capacities and demands
                                {
                                        for(int i=0; i<3; ++i) {
                                                for(int j=0; j<3; ++j) {
                                                        mcf.SetCap(3*i + j, 2);
                                                }
                                        }
                                        for(int i=0; i<3; ++i) {
                                                mcf.AddNodeExcess(i,1);
                                                mcf.AddNodeExcess(3+i,-1);
                                        }

                                        long new_cost = mcf.Solve();
                                        for(int i=0; i<num_arcs; ++i) {
                                                REQUIRE(2*flow[i] == mcf.GetFlow(i));
                                        }
                                        auto nc = mcf.Objective();
                                        assert(nc == new_cost);
                                        assert(2*orig_cost == new_cost);
                                        REQUIRE(2*orig_cost == new_cost);
                                }

                                // halve capacities and demands
                                {
                                        for(int i=0; i<3; ++i) {
                                                for(int j=0; j<3; ++j) { 
                                                        mcf.SetCap(3*i + j, 1);
                                                }
                                        }
                                        for(int i=0; i<3; ++i) {
                                                mcf.AddNodeExcess(i,-1);
                                                mcf.AddNodeExcess(3+i,1);
                                        }

                                        long new_cost = mcf.Solve();
                                        for(int i=0; i<num_arcs; ++i) {
                                                assert(flow[i] == mcf.GetFlow(i));
                                                REQUIRE(flow[i] == mcf.GetFlow(i));
                                        }
                                        auto nc = mcf.Objective();
                                        assert(nc == new_cost);
                                        assert(orig_cost == new_cost);
                                        REQUIRE(orig_cost == new_cost);
                                }
                        }
                }
        }
}
