#include "catch.hpp"
#include <vector>
#include <algorithm>
#include <random>
#include "vector"
#include "include/min_convolution/min_convolution.hxx"

using namespace LP_MP;
//using namespace std;

template<typename VEC1, typename VEC2>
void test_naive_Bussieck(const VEC1 a, const VEC2 b)
{
  auto result_naive = min_convolution::min_conv_naive(a.begin(), a.end(), b.begin(), b.end());
  auto result_Bussieck = min_convolution::min_conv_Bussieck_et_al(a.begin(), a.end(), b.begin(), b.end());
  REQUIRE(result_naive.size() == result_Bussieck.size());
  for(INDEX i=0; i<result_naive.size(); ++i) {
    assert(result_naive[i] == result_Bussieck[i]);
    REQUIRE(result_naive[i] == result_Bussieck[i]);
  }
  //return result_naive == result_Bussieck;
}

TEST_CASE( "min conv", "[min sum convolution]" ) {

  SECTION("artificial input") {
   std::vector<double> a {0.1, 0.2, 0.05, 1};
   std::vector<double> b {0.1, 0.2, 0.05, 1};
   std::reverse(b.begin(), b.end());

   test_naive_Bussieck(a,b);

   //auto left = [&](int i){ return a[i]; };
   //auto right = [&](int i){ return b[i]; };
   //auto op = [=](int i, int j){ return std::min(i+j, 6); };

   auto result = min_convolution::arg_min_conv_Bussieck_et_al(a.begin(), a.end(), b.begin(), b.end(), 7);
   auto result_val = std::get<0>(result);
   REQUIRE(result_val.size() == 7);


   for(auto sum=0;sum<7;sum++){
      //REQUIRE(sum == op(mc.getIdxA(sum),mc.getIdxB(sum)));
      //REQUIRE(sum == (mc.getIdxA(sum) + mc.getIdxB(sum)));

      double val_expl = std::numeric_limits<double>::infinity();
      for(auto i=0; i<=sum; ++i) { 
         if(i<a.size() && sum-i < b.size()) {
            val_expl = std::min(val_expl, a[i] + b[sum-i]); 
         }
      }
      //assert(result_val[sum] == val_expl);
      //REQUIRE(result_val[sum] == val_expl);
   }
  }

  SECTION("random input") {
    // initialize seed explicitly to make unit test reproducible

    // random numbers for size of underlying vectors
    std::mt19937 gen(1); //Standard mersenne_twister_engine seeded with 1
    std::uniform_int_distribution<> dis_int(2, 100);

    // random numbers for vector values
    std::uniform_real_distribution<> dis_real(1.0, 2.0);
    for(INDEX run=0; run<1000; ++run) {
      const INDEX a_size = dis_int(gen); 
      const INDEX b_size = dis_int(gen);

      std::vector<REAL> a(a_size);
      for(INDEX i=0; i<a.size(); ++i) {
        a[i] = dis_real(gen);
      }

      std::vector<REAL> b(b_size);
      for(INDEX i=0; i<b.size(); ++i) {
        b[i] = dis_real(gen);
      }

      test_naive_Bussieck(a,b);

    } 
  }
   /*
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
   */

  // generate random vectors and check whether naive implementation and Bussieck et al algorithms return same
  SECTION("random min convolution") {

  }
}

