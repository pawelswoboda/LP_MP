#include "catch.hpp"
#include "sat_solver.hxx"

using namespace LP_MP;

TEST_CASE( "sat solver", "[lingeling sat solver encapsulation]" ) {

   LP_MP::sat_solver sat;

   SECTION("add literal") {
      auto l = sat.add_literal();
      REQUIRE(sat.size() == 1);
      REQUIRE(l == 1);
   }

   SECTION("add literal vector") {
      auto l = sat.add_literal_vector(9);
      REQUIRE(sat.size() == 9);
      for(INDEX i=0; i<l.size(); ++i) {
         REQUIRE(l[i] == i+1);
      } 
   }

   SECTION("add literal matrix") {
      auto l = sat.add_literal_matrix(4,5);
      REQUIRE(sat.size() == 4*5);
      for(INDEX i=0; i<l.size(); ++i) {
         REQUIRE(l[i] == i+1);
      } 
   }

   SECTION("add literal tensor") {
      auto l = sat.add_literal_tensor(4,5,9);
      REQUIRE(sat.size() == 4*5*9);
      for(INDEX i=0; i<l.size(); ++i) {
         REQUIRE(l[i] == i+1);
      } 
   }

   SECTION("count added literals")
   {
      auto lit = sat.add_literal();
      auto vec = sat.add_literal_vector(11);
      auto matrix = sat.add_literal_matrix(7,2);
      auto tensor = sat.add_literal_tensor(4,9,6);
      REQUIRE(sat.size() == 1 + vec.size() + matrix.size() + tensor.size()); 
   }

   SECTION("feasible CNF") {
      // find out factorization of 7*5 = 35
      sat.add_literal_vector(29);
      sat.add_clause(2, 3 );
      sat.add_clause(5, 6 );
      sat.add_clause(-7 );
      sat.add_clause(8 );
      sat.add_clause(-18, 10, 12 );
      sat.add_clause(-18, -10, -12 );
      sat.add_clause(18, 10, -12 );
      sat.add_clause(18, -10, 12 );
      sat.add_clause(-10, -12, 19 );
      sat.add_clause(10, 12, -19 );
      sat.add_clause(10, -19 );
      sat.add_clause(12, -19 );
      sat.add_clause(-20, 11, 13, 19 );
      sat.add_clause(-20, -11, -13, 19 );
      sat.add_clause(-20, -11, 13, -19 );
      sat.add_clause(-20, 11, -13, -19 );
      sat.add_clause(20, -11, -13, -19 );
      sat.add_clause(20, 11, 13, -19 );
      sat.add_clause(20, 11, -13, 19 );
      sat.add_clause(20, -11, 13, 19 );
      sat.add_clause(-11, -13, 21 );
      sat.add_clause(-11, -19, 21 );
      sat.add_clause(-13, -19, 21 );
      sat.add_clause(11, 13, -21 );
      sat.add_clause(11, 19, -21 );
      sat.add_clause(13, 19, -21 );
      sat.add_clause(-22, 14, 21 );
      sat.add_clause(-22, -14, -21 );
      sat.add_clause(22, 14, -21 );
      sat.add_clause(22, -14, 21 );
      sat.add_clause(-14, -21, 23 );
      sat.add_clause(14, 21, -23 );
      sat.add_clause(14, -23 );
      sat.add_clause(21, -23 );
      sat.add_clause(-24, 20, 15 );
      sat.add_clause(-24, -20, -15 );
      sat.add_clause(24, 20, -15 );
      sat.add_clause(24, -20, 15 );
      sat.add_clause(-20, -15, 25 );
      sat.add_clause(20, 15, -25 );
      sat.add_clause(20, -25 );
      sat.add_clause(15, -25 );
      sat.add_clause(-26, 22, 16, 25 );
      sat.add_clause(-26, -22, -16, 25 );
      sat.add_clause(-26, -22, 16, -25 );
      sat.add_clause(-26, 22, -16, -25 );
      sat.add_clause(26, -22, -16, -25 );
      sat.add_clause(26, 22, 16, -25 );
      sat.add_clause(26, 22, -16, 25 );
      sat.add_clause(26, -22, 16, 25 );
      sat.add_clause(-22, -16, 27 );
      sat.add_clause(-22, -25, 27 );
      sat.add_clause(-16, -25, 27 );
      sat.add_clause(22, 16, -27 );
      sat.add_clause(22, 25, -27 );
      sat.add_clause(16, 25, -27 );
      sat.add_clause(-28, 23, 17, 27 );
      sat.add_clause(-28, -23, -17, 27 );
      sat.add_clause(-28, -23, 17, -27 );
      sat.add_clause(-28, 23, -17, -27 );
      sat.add_clause(28, -23, -17, -27 );
      sat.add_clause(28, 23, 17, -27 );
      sat.add_clause(28, 23, -17, 27 );
      sat.add_clause(28, -23, 17, 27 );
      sat.add_clause(-23, -17, 29 );
      sat.add_clause(-23, -27, 29 );
      sat.add_clause(-17, -27, 29 );
      sat.add_clause(23, 17, -29 );
      sat.add_clause(23, 27, -29 );
      sat.add_clause(17, 27, -29 );
      sat.add_clause(9, -1, -4 );
      sat.add_clause(-9, 1 );
      sat.add_clause(-9, 4 );
      sat.add_clause(10, -2, -4 );
      sat.add_clause(-10, 2 );
      sat.add_clause(-10, 4 );
      sat.add_clause(11, -3, -4 );
      sat.add_clause(-11, 3 );
      sat.add_clause(-11, 4 );
      sat.add_clause(12, -1, -5 );
      sat.add_clause(-12, 1 );
      sat.add_clause(-12, 5 );
      sat.add_clause(13, -2, -5 );
      sat.add_clause(-13, 2 );
      sat.add_clause(-13, 5 );
      sat.add_clause(14, -3, -5 );
      sat.add_clause(-14, 3 );
      sat.add_clause(-14, 5 );
      sat.add_clause(15, -1, -6 );
      sat.add_clause(-15, 1 );
      sat.add_clause(-15, 6 );
      sat.add_clause(16, -2, -6 );
      sat.add_clause(-16, 2 );
      sat.add_clause(-16, 6 );
      sat.add_clause(17, -3, -6 );
      sat.add_clause(-17, 3 );
      sat.add_clause(-17, 6 );
      sat.add_clause(9, -8 );
      sat.add_clause(-9, 8 );
      sat.add_clause(18, -8 );
      sat.add_clause(-18, 8 );
      sat.add_clause(24, -7 );
      sat.add_clause(-24, 7 );
      sat.add_clause(26, -7 );
      sat.add_clause(-26, 7 );
      sat.add_clause(28, -7 );
      sat.add_clause(-28, 7 );
      sat.add_clause(29, -8 );
      sat.add_clause(-29, 8 );

      REQUIRE(sat.solve() == true);
   }

   SECTION("infeasible CNF") {
      // pigeon hole principle: put 6 items in 5 holes
      sat.add_literal_vector(42);
      sat.add_clause(-1,-7);
      sat.add_clause(-1,-13);
      sat.add_clause(-1,-19);
      sat.add_clause(-1,-25);
      sat.add_clause(-1,-31);
      sat.add_clause(-1,-37);
      sat.add_clause(-7,-13);
      sat.add_clause(-7,-19);
      sat.add_clause(-7,-25);
      sat.add_clause(-7,-31);
      sat.add_clause(-7,-37);
      sat.add_clause(-13,-19);
      sat.add_clause(-13,-25);
      sat.add_clause(-13,-31);
      sat.add_clause(-13,-37);
      sat.add_clause(-19,-25);
      sat.add_clause(-19,-31);
      sat.add_clause(-19,-37);
      sat.add_clause(-25,-31);
      sat.add_clause(-25,-37);
      sat.add_clause(-31,-37);
      sat.add_clause(-2,-8);
      sat.add_clause(-2,-14);
      sat.add_clause(-2,-20);
      sat.add_clause(-2,-26);
      sat.add_clause(-2,-32);
      sat.add_clause(-2,-38);
      sat.add_clause(-8,-14);
      sat.add_clause(-8,-20);
      sat.add_clause(-8,-26);
      sat.add_clause(-8,-32);
      sat.add_clause(-8,-38);
      sat.add_clause(-14,-20);
      sat.add_clause(-14,-26);
      sat.add_clause(-14,-32);
      sat.add_clause(-14,-38);
      sat.add_clause(-20,-26);
      sat.add_clause(-20,-32);
      sat.add_clause(-20,-38);
      sat.add_clause(-26,-32);
      sat.add_clause(-26,-38);
      sat.add_clause(-32,-38);
      sat.add_clause(-3,-9);
      sat.add_clause(-3,-15);
      sat.add_clause(-3,-21);
      sat.add_clause(-3,-27);
      sat.add_clause(-3,-33);
      sat.add_clause(-3,-39);
      sat.add_clause(-9,-15);
      sat.add_clause(-9,-21);
      sat.add_clause(-9,-27);
      sat.add_clause(-9,-33);
      sat.add_clause(-9,-39);
      sat.add_clause(-15,-21);
      sat.add_clause(-15,-27);
      sat.add_clause(-15,-33);
      sat.add_clause(-15,-39);
      sat.add_clause(-21,-27);
      sat.add_clause(-21,-33);
      sat.add_clause(-21,-39);
      sat.add_clause(-27,-33);
      sat.add_clause(-27,-39);
      sat.add_clause(-33,-39);
      sat.add_clause(-4,-10);
      sat.add_clause(-4,-16);
      sat.add_clause(-4,-22);
      sat.add_clause(-4,-28);
      sat.add_clause(-4,-34);
      sat.add_clause(-4,-40);
      sat.add_clause(-10,-16);
      sat.add_clause(-10,-22);
      sat.add_clause(-10,-28);
      sat.add_clause(-10,-34);
      sat.add_clause(-10,-40);
      sat.add_clause(-16,-22);
      sat.add_clause(-16,-28);
      sat.add_clause(-16,-34);
      sat.add_clause(-16,-40);
      sat.add_clause(-22,-28);
      sat.add_clause(-22,-34);
      sat.add_clause(-22,-40);
      sat.add_clause(-28,-34);
      sat.add_clause(-28,-40);
      sat.add_clause(-34,-40);
      sat.add_clause(-5,-11);
      sat.add_clause(-5,-17);
      sat.add_clause(-5,-23);
      sat.add_clause(-5,-29);
      sat.add_clause(-5,-35);
      sat.add_clause(-5,-41);
      sat.add_clause(-11,-17);
      sat.add_clause(-11,-23);
      sat.add_clause(-11,-29);
      sat.add_clause(-11,-35);
      sat.add_clause(-11,-41);
      sat.add_clause(-17,-23);
      sat.add_clause(-17,-29);
      sat.add_clause(-17,-35);
      sat.add_clause(-17,-41);
      sat.add_clause(-23,-29);
      sat.add_clause(-23,-35);
      sat.add_clause(-23,-41);
      sat.add_clause(-29,-35);
      sat.add_clause(-29,-41);
      sat.add_clause(-35,-41);
      sat.add_clause(-6,-12);
      sat.add_clause(-6,-18);
      sat.add_clause(-6,-24);
      sat.add_clause(-6,-30);
      sat.add_clause(-6,-36);
      sat.add_clause(-6,-42);
      sat.add_clause(-12,-18);
      sat.add_clause(-12,-24);
      sat.add_clause(-12,-30);
      sat.add_clause(-12,-36);
      sat.add_clause(-12,-42);
      sat.add_clause(-18,-24);
      sat.add_clause(-18,-30);
      sat.add_clause(-18,-36);
      sat.add_clause(-18,-42);
      sat.add_clause(-24,-30);
      sat.add_clause(-24,-36);
      sat.add_clause(-24,-42);
      sat.add_clause(-30,-36);
      sat.add_clause(-30,-42);
      sat.add_clause(-36,-42);
      sat.add_clause( 6,5,4,3,2,1);
      sat.add_clause( 12,11,10,9,8,7);
      sat.add_clause( 18,17,16,15,14,13);
      sat.add_clause( 24,23,22,21,20,19);
      sat.add_clause( 30,29,28,27,26,25);
      sat.add_clause( 36,35,34,33,32,31);
      sat.add_clause( 42,41,40,39,38,37);

      REQUIRE(sat.solve() == false);
   }

   SECTION("simplex constraint") {
      auto literals = sat.add_literal_vector(10);
      sat.add_simplex_constraint(literals.begin(),literals.end());
      REQUIRE(sat.solve() == true);
      // forbid previous assignment and check whether new assignment can be found
      sat_vec assumptions;
      for(INDEX t=0; t<10; ++t) {
        if(!assumptions.empty()) {
          REQUIRE(sat.solve(assumptions.begin(), assumptions.end()) == true);
        } else {
          REQUIRE(sat.solve() == true);
        }
        auto sol = sat.solution();
        REQUIRE(std::count(sol.begin(), sol.begin()+10, true) == 1); 
        const auto sol_lit = literals[std::find(sol.begin(), sol.begin()+10, true) - sol.begin()];
        assumptions.push_back(-sol_lit);
      }
      assert(assumptions.size() == 10);
      REQUIRE(sat.solve(assumptions.begin(), assumptions.end()) == false);
   }

   SECTION("at most one constraint") {
      auto literals = sat.add_literal_vector(10);
      sat.add_at_most_one_constraint(literals.begin(),literals.end());
      REQUIRE(sat.solve() == true);
      auto sol = sat.solution();
      REQUIRE(std::count(sol.begin(), sol.begin()+10, true) <= 1); 
      // forbid previous assignment and check whether new assignment can be found
      sat_vec assumptions;
      for(INDEX t=0; t<10; ++t) {
        assumptions.push_back(-literals[t]);
      }
      assert(assumptions.size() == 10);
      REQUIRE(sat.solve(assumptions.begin(), assumptions.end()) == true);
   }

   SECTION("load literals") {
      auto vec_literal = sat.add_literal_vector(9);
      auto matrix_literal = sat.add_literal_matrix(9,3);
      auto tensor_literal = sat.add_literal_tensor(3,4,5);

      decltype(vec_literal) vec_literal_l(vec_literal.size());
      decltype(matrix_literal) matrix_literal_l(matrix_literal.dim(0), matrix_literal.dim(1));
      decltype(tensor_literal) tensor_literal_l(tensor_literal.dim(0), tensor_literal.dim(1), tensor_literal.dim(2));

      load_sat_literals(*(vec_literal.begin()), vec_literal_l, matrix_literal_l, tensor_literal_l);
      REQUIRE(vec_literal == vec_literal_l);
      REQUIRE(matrix_literal == matrix_literal_l);
      REQUIRE(tensor_literal == tensor_literal_l);
   } 

   SECTION("matrix slices") {
      auto matrix_literal = sat.add_literal_matrix(5,3);

      auto slice_left = matrix_literal.slice_left(1);
      REQUIRE(slice_left.size() == 3);
      REQUIRE(slice_left[0] == matrix_literal(1,0));
      REQUIRE(slice_left[1] == matrix_literal(1,1));
      REQUIRE(slice_left[2] == matrix_literal(1,2));

      auto slice_right = matrix_literal.slice_right(2);
      REQUIRE(slice_right.size() == 5);
      REQUIRE(slice_right[0] == matrix_literal(0,2));
      REQUIRE(slice_right[1] == matrix_literal(1,2));
      REQUIRE(slice_right[2] == matrix_literal(2,2));
      REQUIRE(slice_right[3] == matrix_literal(3,2));
      REQUIRE(slice_right[4] == matrix_literal(4,2));
   }

   SECTION("tensor slices") {
      auto tensor_literal = sat.add_literal_tensor(3,4,5);

      auto slice12 = tensor_literal.slice12(0,1);
      REQUIRE(slice12.size() == tensor_literal.dim(2));
      REQUIRE(slice12[0] == tensor_literal(0,1,0));
      REQUIRE(slice12[1] == tensor_literal(0,1,1));
      REQUIRE(slice12[2] == tensor_literal(0,1,2));
      REQUIRE(slice12[3] == tensor_literal(0,1,3));
      REQUIRE(slice12[4] == tensor_literal(0,1,4));

      auto slice13 = tensor_literal.slice13(1,2);
      REQUIRE(slice13.size() == tensor_literal.dim(1));
      REQUIRE(slice13[0] == tensor_literal(1,0,2));
      REQUIRE(slice13[1] == tensor_literal(1,1,2));
      REQUIRE(slice13[2] == tensor_literal(1,2,2));
      REQUIRE(slice13[3] == tensor_literal(1,3,2));

      auto slice23 = tensor_literal.slice23(3,4);
      REQUIRE(slice23.size() == tensor_literal.dim(0));
      REQUIRE(slice23[0] == tensor_literal(0,3,4));
      REQUIRE(slice23[1] == tensor_literal(1,3,4));
      REQUIRE(slice23[2] == tensor_literal(2,3,4));
   }
}

