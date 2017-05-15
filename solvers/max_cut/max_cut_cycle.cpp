
#include "max_cut.h"
#include "visitors/standard_visitor.hxx"
using namespace LP_MP;
int main(int argc, char* argv[])

{
ProblemConstructorRoundingSolver<Solver<FMC_CYCLE_MAX_CUT,LP,StandardTighteningVisitor>> solver(argc,argv);
solver.ReadProblem(max_cut_simple_text_format::ParseProblemMaxCut<Solver<FMC_CYCLE_MAX_CUT,LP,StandardTighteningVisitor>>);
return solver.Solve();

}
