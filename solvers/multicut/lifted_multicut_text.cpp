
#include "multicut.h"
#include "visitors/standard_visitor.hxx"
using namespace LP_MP;
int main(int argc, char* argv[])

{
CombinedMPProblemConstructorRoundingSolver<Solver<FMC_LIFTED_MULTICUT,LP,StandardTighteningVisitor>> solver(argc,argv);
solver.ReadProblem(MulticutTextInput::ParseLiftedProblem<Solver<FMC_LIFTED_MULTICUT,LP,StandardTighteningVisitor>>);
return solver.Solve();

}
