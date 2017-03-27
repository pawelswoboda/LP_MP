//#include "graph_matching.h"
//#include "visitors/standard_visitor.hxx"
//using FMC_INST = FMC_MP_T<PairwiseConstruction::Left>;
//using BaseSolverType = Solver<FMC_INST,LP,StandardTighteningVisitor>;
//LP_MP_CONSTRUCT_SOLVER_WITH_INPUT_AND_VISITOR_MP_ROUNDING(FMC_INST, TorresaniEtAlInput::ParseProblemMP<BaseSolverType>, StandardTighteningVisitor);

#include <iostream>
#include <array>
#include <utility>

template<int M, int N>
struct compile_time_matrix
{
   constexpr int operator()(int m, int n) const { return p_[m*N + n]; }
   int p_[M*N];
};

constexpr compile_time_matrix<3,1> kwas {{0,0,0}};
// make int sequence out of array
constexpr auto make_integer_sequence(const compile_time_matrix<3,1> a)
{
   return std::integer_sequence<int, a(0,0), a(1,0), a(2,0)>{};
}

int main()
{
   std::cout << make_integer_sequence(kwas) << "\n";
   return 0;
}
