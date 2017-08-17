#include "catch.hpp"
#include "serialization.hxx"

using namespace LP_MP;

TEST_CASE( "serialization", "[serialization]" ) {
   INDEX i=10;
   std::array<INDEX,2> a{20,30};
   std::vector<INDEX> v{30,40,50};
   INDEX p [4] = {60,70,80,90};

   REAL i_r=10;
   std::array<REAL,2> a_r{20,30};
   std::vector<REAL> v_r{30,40,50};
   REAL p_r [4] = {60,70,80,90};

   allocate_archive a_ar;
   a_ar(i);
   a_ar(a);
   a_ar(v);
   a_ar( binary_data<INDEX>(p,4) );

   SECTION("allocate archive") {
      {
         allocate_archive ar;
         ar( i );
         REQUIRE(ar.size() == sizeof(INDEX));
      }

      {
         allocate_archive ar;
         ar( a );
         REQUIRE(ar.size() == 2*sizeof(INDEX));
      }

      {
         allocate_archive ar;
         ar( v );
         REQUIRE(ar.size() == 3*sizeof(INDEX));
      }

      {
         allocate_archive ar;
         ar( binary_data<INDEX>(p,4) );
         REQUIRE(ar.size() == 4*sizeof(INDEX));
      }

      {
         allocate_archive ar;
         ar( i_r );
         REQUIRE(ar.size() == sizeof(REAL));
      }

      {
         allocate_archive ar;
         ar( a_r );
         REQUIRE(ar.size() == 2*sizeof(REAL));
      }

      {
         allocate_archive ar;
         ar( v_r );
         REQUIRE(ar.size() == 3*sizeof(REAL));
      }

      {
         allocate_archive ar;
         ar( binary_data<REAL>(p_r,4) );
         REQUIRE(ar.size() == 4*sizeof(REAL));
      }

   }


   SECTION("individual saving") {
      serialization_archive ar(a_ar);

      save_archive s_ar(ar);
      s_ar(i);
      s_ar(a);
      s_ar(v);
      s_ar( binary_data<INDEX>(p,4) );

      load_archive l_ar(ar);

      decltype(i) i_test;
      l_ar(i_test);
      REQUIRE(i_test == i);

      decltype(a) a_test;
      l_ar(a_test);
      REQUIRE(a_test == a);

      decltype(v) v_test(v.size());
      l_ar(v_test);
      assert(v_test == v);
      REQUIRE(v_test == v);

      decltype(p) p_test = {0,0,0,0};
      l_ar( binary_data<INDEX>(p_test, 4) );
      REQUIRE(p_test[0] == p[0]);
      REQUIRE(p_test[1] == p[1]);
      REQUIRE(p_test[2] == p[2]);
      REQUIRE(p_test[3] == p[3]);
   }

   SECTION("collective saving") {
      serialization_archive ar(a_ar);

      save_archive s_ar(ar);
      s_ar(i, a, v, binary_data<INDEX>(p,4));

      load_archive l_ar(ar);

      decltype(i) i_test;
      decltype(a) a_test;
      decltype(v) v_test(v.size());
      decltype(p) p_test = {0,0,0,0};

      l_ar(i_test, a_test, v_test, binary_data<INDEX>(p_test,4));

      REQUIRE(i_test == i);
      REQUIRE(a_test == a);
      REQUIRE(v_test == v);
      REQUIRE(p_test[0] == p[0]);
      REQUIRE(p_test[1] == p[1]);
      REQUIRE(p_test[2] == p[2]);
      REQUIRE(p_test[3] == p[3]);
   }
   
}

