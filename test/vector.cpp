#include "catch.hpp"
#include "vector.hxx"

using namespace LP_MP;

TEST_CASE("vector", "[vector]") {

  SECTION("vector minimum") {
    vector<REAL> v(5);
    v[0] = -1.0;
    v[1] = 0.0;
    v[2] = 1.0;
    v[3] = 2.0;
    v[4] = 3.0;

    REQUIRE(v.min() == -1.0); 
  }

  {
    matrix<REAL> m(5,6);
    m(0,0) = -2.0; m(0,1) = +0.0; m(0,2) = +2.0; m(0,3) = -0.5; m(0,4) = +0.0; m(0,5) = +0.5;
    m(1,0) = -1.0; m(1,1) = +0.0; m(1,2) = +1.0; m(1,3) = -0.5; m(1,4) = +0.0; m(1,5) = +0.5;
    m(2,0) = -0.0; m(2,1) = -4.0; m(2,2) = +0.5; m(2,3) = -0.5; m(2,4) = +0.0; m(2,5) = +0.5;
    m(3,0) = +1.0; m(3,1) = +0.0; m(3,2) = -1.0; m(3,3) = -0.5; m(3,4) = +0.0; m(3,5) = +0.5;
    m(4,0) = +2.0; m(4,1) = +0.0; m(4,2) = -2.0; m(4,3) = -0.5; m(4,4) = +0.0; m(4,5) = +0.5;

    SECTION("matrix columnwise minimum") {
      auto min_col = m.min2(); 
      REQUIRE(min_col.size() == 6);

      REQUIRE(min_col[0] == -2.0);
      REQUIRE(min_col[1] == -4.0);
      REQUIRE(min_col[2] == -2.0); 
      REQUIRE(min_col[3] == -0.5); 
      REQUIRE(min_col[4] == +0.0); 
      REQUIRE(min_col[5] == +0.5); 
    }

    SECTION("matrix rowwise minimum") {
      auto min_row = m.min1();
      REQUIRE(min_row.size() == 5);

      REQUIRE(min_row[0] == -2.0);
      REQUIRE(min_row[1] == -1.0);
      REQUIRE(min_row[2] == -4.0);
      REQUIRE(min_row[3] == -1.0);
      REQUIRE(min_row[4] == -2.0); 
    }
  } 
}

