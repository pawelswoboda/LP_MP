#include "catch.hpp"
#include <vector>
#include <algorithm>
#include "solvers/discrete_tomography/minConv.hxx"

using namespace LP_MP;
using namespace std;

TEST_CASE( "min conv", "[min sum convolution]" ) {
   std::vector<double> a {0.1, 0.2, 0.05, 1};
   std::vector<double> b {0.1, 0.2, 0.05, 1};
   std::reverse(b.begin(), b.end());

   auto left = [&](int i){ return a[i]; };
   auto right = [&](int i){ return b[i]; };
   auto op = [=](int i, int j){ return std::min(i+j, 6); };

   typename LP_MP::discrete_tomo::MinConv<double, size_t> mc(left,right,a.size(),b.size(),6);
   mc.CalcConv(op,left,right);

   for(auto sum=0;sum<6;sum++){
      REQUIRE(sum == op(mc.getIdxA(sum),mc.getIdxB(sum)));
      REQUIRE(sum == (mc.getIdxA(sum) + mc.getIdxB(sum)));

      double val_expl = std::numeric_limits<double>::infinity();
      for(auto i=0; i<=sum; ++i) { 
         if(i<a.size() && sum-i < b.size()) {
            val_expl = std::min(val_expl, a[i] + b[sum-i]); 
         }
      }
      assert(mc.getConv(sum) == val_expl);
      REQUIRE(mc.getConv(sum) == val_expl);
   }

   SECTION("difference") {
      auto op = [=](size_t i, size_t j){ return i-j; };

      typename LP_MP::discrete_tomo::MinConv<double, size_t> mc(left,right,a.size(),b.size(),4);
      mc.CalcConv(op,left,right);

      for(auto sum=0;sum<4;sum++){
         REQUIRE(sum == op(mc.getIdxA(sum),mc.getIdxB(sum)));

         double val_expl = std::numeric_limits<double>::infinity();
         for(auto i=0; i<4-sum; ++i) { 
            val_expl = std::min(val_expl, a[sum+i] + b[i]); 
         }
         assert(mc.getConv(sum) == val_expl);
         REQUIRE(mc.getConv(sum) == val_expl);
      }


   }
}

