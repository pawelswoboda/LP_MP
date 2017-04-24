
#include "multicut.h"
#include "visitors/standard_visitor.hxx"
using namespace LP_MP;
int main(int argc, char* argv[])

{
ProblemConstructorRoundingSolver<Solver<FMC_MULTIWAY_CUT,LP,StandardTighteningVisitor>> solver(argc,argv);
solver.ReadProblem(MulticutOpenGmInput::ParsePottsProblem<Solver<FMC_MULTIWAY_CUT,LP,StandardTighteningVisitor>>);
return solver.Solve();

}
