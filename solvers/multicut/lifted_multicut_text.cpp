#include "multicut.h"
#include "visitors/standard_visitor.hxx"
using namespace LP_MP;
using FMC = FMC_LIFTED_MULTICUT;
//LP_MP_CONSTRUCT_SOLVER_WITH_INPUT_AND_VISITOR(FMC, MulticutTextInput::ParseLiftedProblem<FMC>, StandardTighteningVisitor);
using SolverType = ProblemConstructorRoundingSolver<Solver<FMC,LP,StandardTighteningVisitor>>;
int main(int argc, char* argv[])
{
   SolverType solver(argc,argv);
   solver.ReadProblem(MulticutTextInput::ParseLiftedProblem<Solver<FMC,LP,StandardTighteningVisitor>>);
   return solver.Solve();
}



