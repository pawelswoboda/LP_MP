
#include "multicut.h"
#include "visitors/standard_visitor.hxx"
using namespace LP_MP;
int main(int argc, char* argv[])

{
CombinedMPProblemConstructorRoundingSolver<Solver<FMC_LIFTED_MULTICUT,LP,StandardTighteningVisitor>> solver(argc,argv);
solver.ReadProblem(MulticutH5Input::ParseLiftedProblem<Solver<FMC_LIFTED_MULTICUT,LP,StandardTighteningVisitor>,false>);
return solver.Solve();

}
