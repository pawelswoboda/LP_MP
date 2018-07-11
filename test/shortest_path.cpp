#include "catch.hpp"
#include <vector>
#include "solvers/multicut/multicut.h"
#include "solvers/multicut/multicut_constructor.hxx"

using namespace LP_MP;

TEST_CASE( "shortest path search in multicut factor", "[shortest path]" ) {

   using mc_constructor = MulticutConstructor<FMC_MULTICUT<MessageSendingType::SRMP>, 0, 1, 0, 1, 2, 2>;
   //using graph = mc_constructor::Graph2;
   //using bfs_data = mc_constructor::BfsData2; 
   
   Graph g(4,8, std::vector<INDEX>{2,2,2,2});

   g.add_arc(0,1,1.0);
   g.add_arc(1,0,1.0);

   g.add_arc(1,2,2.0);
   g.add_arc(2,1,2.0);

   g.add_arc(2,3,2.0);
   g.add_arc(3,2,2.0);

   g.add_arc(3,0,2.0);
   g.add_arc(0,3,2.0);

   g.sort();
   BfsData sp(g);

   auto c1 = sp.FindPath(0,2,g);
   assert(std::get<1>(c1).size() == 3);
   REQUIRE(std::get<1>(c1).size() == 3);

   auto c2 = sp.FindPath(0,1,g, 0.5);
   assert(std::get<1>(c2).size() == 2);
   assert(std::get<0>(c2) == 1.0);
   REQUIRE(std::get<1>(c2).size() == 2);

   auto c3 = sp.FindPath(0,1,g, 1.5);
   assert(std::get<1>(c3).size() == 4);
   assert(std::get<0>(c3) == 2.0);
   REQUIRE(std::get<1>(c3).size() == 4);

}
