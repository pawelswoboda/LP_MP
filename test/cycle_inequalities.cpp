#include "catch.hpp"
#include <vector>
#include "solvers/graphical_model/graphical_model.h"
#include "solver.hxx"
#include "visitors/standard_visitor.hxx"
#include "problem_constructors/cycle_inequalities.hxx"

using namespace LP_MP;

TEST_CASE( "cycle inequalities tightening for MAP-MRF", "[MAP-MRF tightening]" ) {

   SECTION("conditional minima") {
      matrix<REAL> m(4,4);
      INDEX c=0;
      for(INDEX x1=0; x1<4; ++x1) {
         for(INDEX x2=0; x2<4; ++x2) {
            m(x1,x2) = c;
            ++c;
         }
      }
      const auto row_minima = k_ary_cycle_inequalities_search<int>::row_minima(m);
      const auto column_minima = k_ary_cycle_inequalities_search<int>::column_minima(m);
      const auto principal_minima = k_ary_cycle_inequalities_search<int>::principal_minima(m, column_minima);

      REQUIRE(row_minima(0,0) == 0);
      REQUIRE(row_minima(0,1) == 1);
      REQUIRE(row_minima(1,0) == 4);
      REQUIRE(row_minima(1,1) == 5);
      REQUIRE(row_minima(2,0) == 8);
      REQUIRE(row_minima(2,1) == 9);
      REQUIRE(row_minima(3,0) == 12);
      REQUIRE(row_minima(3,1) == 13);

      REQUIRE(column_minima(0,0) == 0);
      REQUIRE(column_minima(0,1) == 4);
      REQUIRE(column_minima(1,0) == 1);
      REQUIRE(column_minima(1,1) == 5);
      REQUIRE(column_minima(2,0) == 2);
      REQUIRE(column_minima(2,1) == 6);
      REQUIRE(column_minima(3,0) == 3);
      REQUIRE(column_minima(3,1) == 7);

      REQUIRE(principal_minima(0,0) == 5);
      REQUIRE(principal_minima(0,1) == 4);
      REQUIRE(principal_minima(0,2) == 4);
      REQUIRE(principal_minima(0,3) == 4);

      REQUIRE(principal_minima(1,0) == 1);
      REQUIRE(principal_minima(1,1) == 0);
      REQUIRE(principal_minima(1,2) == 0);
      REQUIRE(principal_minima(1,3) == 0);

      REQUIRE(principal_minima(2,0) == 1);
      REQUIRE(principal_minima(2,1) == 0);
      REQUIRE(principal_minima(2,2) == 0);
      REQUIRE(principal_minima(2,3) == 0);
   
      REQUIRE(principal_minima(3,0) == 1);
      REQUIRE(principal_minima(3,1) == 0);
      REQUIRE(principal_minima(3,2) == 0);
      REQUIRE(principal_minima(3,3) == 0);
   }

   char * i[5];
   i[0] = "";
   i[1] = "-i";
   i[2] = "";
   i[3] = "-v";
   i[4] = "2";
   Solver<FMC_SRMP_T,LP,StandardVisitor> s(5,i);
   auto& mrf = s.template GetProblemConstructor<0>();
   s.GetLP().set_reparametrization(LPReparametrizationMode::DampedUniform); // setting reparametrization mode to anisotropic leads to suboptimal fixed point

   matrix<REAL> negPotts2(2,2);
   negPotts2(0,0) = 1.0;
   negPotts2(1,1) = 1.0;
   negPotts2(0,1) = 0.0;
   negPotts2(1,0) = 0.0;

   matrix<REAL> posPotts2(2,2); 
   posPotts2(0,0) = 0.0;
   posPotts2(1,1) = 0.0;
   posPotts2(0,1) = 1.0;
   posPotts2(1,0) = 1.0;

   SECTION("binary violated 4 cycle") {
      mrf.AddUnaryFactor({0,0});
      mrf.AddUnaryFactor({0,0});
      mrf.AddUnaryFactor({0,0});
      mrf.AddUnaryFactor({0,0});

      
      mrf.AddPairwiseFactor(0,1,negPotts2);
      mrf.AddPairwiseFactor(1,2,posPotts2);
      mrf.AddPairwiseFactor(2,3,posPotts2);
      mrf.AddPairwiseFactor(0,3,posPotts2);

      k_ary_cycle_inequalities_search<typename std::remove_reference<decltype(mrf)>::type> cycle_search(mrf);

      auto triplets = cycle_search.search();
      REQUIRE(triplets.size() >= 2);
      mrf.add_triplets(triplets);

      for(INDEX i=0; i<100; ++i) {
         s.GetLP().ComputePass(i);
      }
      REQUIRE(s.GetLP().LowerBound() > 1.0-eps);
   }

   SECTION("binary violated 5 cycle") {
      mrf.AddUnaryFactor({0,0});
      mrf.AddUnaryFactor({0,0});
      mrf.AddUnaryFactor({0,0});
      mrf.AddUnaryFactor({0,0});
      mrf.AddUnaryFactor({0,0});

      
      mrf.AddPairwiseFactor(0,1,negPotts2);
      mrf.AddPairwiseFactor(1,2,posPotts2);
      mrf.AddPairwiseFactor(2,3,posPotts2);
      mrf.AddPairwiseFactor(3,4,posPotts2);
      mrf.AddPairwiseFactor(0,4,posPotts2);

      k_ary_cycle_inequalities_search<typename std::remove_reference<decltype(mrf)>::type> cycle_search(mrf);

      auto triplets = cycle_search.search();

      REQUIRE(triplets.size() >= 3);
      mrf.add_triplets(triplets);

      for(INDEX i=0; i<100; ++i) {
         s.GetLP().ComputePass(i);
      }
      REQUIRE(s.GetLP().LowerBound() > 1.0-eps);
      
   }

   matrix<REAL> negPotts4(4,4);
   negPotts4(0,0) = 1.0; negPotts4(0,1) = 1.0; negPotts4(0,2) = 0.0; negPotts4(0,3) = 0.0;
   negPotts4(1,0) = 1.0; negPotts4(1,1) = 1.0; negPotts4(1,2) = 0.0; negPotts4(1,3) = 0.0;
   negPotts4(2,0) = 0.0; negPotts4(2,1) = 0.0; negPotts4(2,2) = 1.0; negPotts4(2,3) = 1.0;
   negPotts4(3,0) = 0.0; negPotts4(3,1) = 0.0; negPotts4(3,2) = 1.0; negPotts4(3,3) = 1.0;

   matrix<REAL> posPotts4(4,4);
   posPotts4(0,0) = 0.0; posPotts4(0,1) = 0.0; posPotts4(0,2) = 1.0; posPotts4(0,3) = 1.0;
   posPotts4(1,0) = 0.0; posPotts4(1,1) = 0.0; posPotts4(1,2) = 1.0; posPotts4(1,3) = 1.0;
   posPotts4(2,0) = 1.0; posPotts4(2,1) = 1.0; posPotts4(2,2) = 0.0; posPotts4(2,3) = 0.0;
   posPotts4(3,0) = 1.0; posPotts4(3,1) = 1.0; posPotts4(3,2) = 0.0; posPotts4(3,3) = 0.0;

   SECTION("expanded k-ary cycle search 4-cycle") {
      mrf.AddUnaryFactor({0,0,0,0});
      mrf.AddUnaryFactor({0,0,0,0});
      mrf.AddUnaryFactor({0,0,0,0});
      mrf.AddUnaryFactor({0,0,0,0});

      
      mrf.AddPairwiseFactor(0,1,negPotts4);
      mrf.AddPairwiseFactor(1,2,posPotts4);
      mrf.AddPairwiseFactor(2,3,posPotts4);
      mrf.AddPairwiseFactor(0,3,posPotts4);

      k_ary_cycle_inequalities_search<typename std::remove_reference<decltype(mrf)>::type,false> cycle_search(mrf);
      k_ary_cycle_inequalities_search<typename std::remove_reference<decltype(mrf)>::type,true> cycle_search2(mrf);

      auto triplets = cycle_search.search();
      REQUIRE(triplets.size() == 0);
      triplets = cycle_search2.search();
      REQUIRE(triplets.size() >= 2);

      mrf.add_triplets(triplets);

      for(INDEX i=0; i<100; ++i) {
         s.GetLP().ComputePass(i);
      }
      REQUIRE(s.GetLP().LowerBound() > 1.0-eps);
   }

   SECTION("expanded k-ary cycle search 5-cycle") {
      mrf.AddUnaryFactor({0,0,0,0});
      mrf.AddUnaryFactor({0,0,0,0});
      mrf.AddUnaryFactor({0,0,0,0});
      mrf.AddUnaryFactor({0,0,0,0});
      mrf.AddUnaryFactor({0,0,0,0});

      
      mrf.AddPairwiseFactor(0,1,negPotts4);
      mrf.AddPairwiseFactor(1,2,posPotts4);
      mrf.AddPairwiseFactor(2,3,posPotts4);
      mrf.AddPairwiseFactor(3,4,posPotts4);
      mrf.AddPairwiseFactor(0,4,posPotts4);

      k_ary_cycle_inequalities_search<typename std::remove_reference<decltype(mrf)>::type,false> cycle_search(mrf);
      k_ary_cycle_inequalities_search<typename std::remove_reference<decltype(mrf)>::type,true> cycle_search2(mrf);

      auto triplets = cycle_search.search();
      REQUIRE(triplets.size() == 0);
      triplets = cycle_search2.search();
      REQUIRE(triplets.size() >= 3);

      mrf.add_triplets(triplets);

      for(INDEX i=0; i<100; ++i) {
         s.GetLP().ComputePass(i);
      }
      REQUIRE(s.GetLP().LowerBound() > 1.0-eps);
   }
}

