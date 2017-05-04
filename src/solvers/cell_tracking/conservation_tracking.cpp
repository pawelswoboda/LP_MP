
#include "solvers/cell_tracking/cell_tracking.h"
#include "visitors/standard_visitor.hxx"
using namespace LP_MP;
int main(int argc, char* argv[])

{
MpRoundingSolver<Solver<FMC_CONSERVATION_TRACKING,LP_sat<LP>,StandardVisitor>> solver(argc,argv);
solver.ReadProblem(conservation_tracking_parser::ParseProblem<Solver<FMC_CONSERVATION_TRACKING,LP_sat<LP>,StandardVisitor>>);
return solver.Solve();

}
