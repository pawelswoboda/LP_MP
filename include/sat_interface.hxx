#ifndef LP_MP_SAT_INTETRFACE
#define LP_MP_SAT_INTETRFACE

#include "config.hxx"
//#include "simp/SimpSolver.h"
//#include "cryptominisat5/cryptominisat.h"
extern "C" {
#include "lglib.h"
}

namespace LP_MP {

  using sat_var = int;
  //using sat_var = Glucose::Var;
  using sat_literal = int;
  //using sat_literal = CMSat::Lit;
  template<typename T>
  using sat_vec = std::vector<T>;

  sat_literal to_literal(const sat_var v) 
  {
    return v+1;
    //return CMSat::Lit(v,false);
    //return CMSat::mkLit(v,false);
  }

  template<typename SAT_SOLVER>
  void make_sat_var_equal(SAT_SOLVER& s, const sat_literal i, const sat_literal j)
  {
    lgladd(s, i); lgladd(s, -j); lgladd(s, 0); 
    lgladd(s, -i); lgladd(s, j); lgladd(s, 0); 
      //s.add_clause({i,-j});
      //s.add_clause({-i,j}); 
  }

  template<typename SAT_SOLVER>
  std::vector<sat_var> create_sat_variables(SAT_SOLVER& s, const INDEX n)
  {
    std::vector<sat_var> v(n);
    for(INDEX i=0; i<n; ++i) {
      //v.push_back(s.nVars());
      //s.new_var(); 
      v.push_back( lglincvar(s)-1 );
    }
    return std::move(v);
  }
  template<typename SAT_SOLVER, typename ITERATOR>
  sat_var add_at_most_one_constraint_naive_sat(SAT_SOLVER& s, ITERATOR var_begin, ITERATOR var_end)
  {
    //assert(1 < std::distance(var_begin, var_end));
    const INDEX n = std::distance(var_begin, var_end);
    if(n == 1) {
      return *var_begin;
    }
    for(INDEX i=0; i<n; ++i) {
      for(INDEX j=i+1; j<n; ++j) {
        //s.add_clause({-to_literal(*(var_begin+i)), -to_literal(*(var_begin+j))});
        lgladd(s, -to_literal(*(var_begin+i)));
        lgladd(s, -to_literal(*(var_begin+j)));
        lgladd(s, 0);
      }
    }
    auto c = lglmaxvar(s);
    //auto c = s.nVars();
    //s.new_var();
    for(INDEX i=0; i<n; ++i) {
      //s.add_clause({to_literal(c),-to_literal(*(var_begin+i))});
      lgladd(s, to_literal(c));
      lgladd(s, -to_literal(*(var_begin+i)));
      lgladd(s, 0);
    }
    

    lgladd(s, -to_literal(c));
    for(INDEX i=0; i<n; ++i) {
      lgladd(s, to_literal(*(var_begin+i)));
    }
    lgladd(s, 0);
    //sat_vec<sat_literal> v({-to_literal(c)});
    //for(INDEX i=0; i<n; ++i) {
    //  v.push_back(to_literal(*(var_begin+i)));
    //}
    //s.add_clause(v); 

    return c; 
  }

  template<typename SAT_SOLVER, typename ITERATOR>
  sat_var add_at_most_one_constraint_sat(SAT_SOLVER& s, ITERATOR var_begin, ITERATOR var_end)
  {
    constexpr INDEX th = 3; // pursue a direct approach until 6 variables, and if more a recursive one
    const INDEX n = std::distance(var_begin, var_end);
    assert(n > 0);
    if(n <= th) {
      return add_at_most_one_constraint_naive_sat(s, var_begin, var_end);
    } else {
      if(n <= 3*th) {
        auto c1 = add_at_most_one_constraint_sat(s, var_begin, var_begin + n/2);
        auto c2 = add_at_most_one_constraint_sat(s, var_begin + n/2, var_end);
        std::array<sat_var,2> c_list {c1,c2};
        return add_at_most_one_constraint_naive_sat(s, c_list.begin(), c_list.end()); 
      } else {
        auto c1 = add_at_most_one_constraint_sat(s, var_begin, var_begin + n/3);
        auto c2 = add_at_most_one_constraint_sat(s, var_begin + n/3, var_begin + 2*n/3);
        auto c3 = add_at_most_one_constraint_sat(s, var_begin  + 2*n/3, var_end);
        std::array<sat_var,3> c_list {c1,c2,c3};
        return add_at_most_one_constraint_naive_sat(s, c_list.begin(), c_list.end());
      }
    }
  }



} // end namespace LP_MP

#endif
