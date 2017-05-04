
#include "cell_tracking.h"
#include "visitors/standard_visitor.hxx"
using namespace LP_MP;
int main(int argc, char* argv[])

{
MpRoundingSolver<Solver<FMC_CELL_TRACKING_MOTHER_MACHINE,LP_sat<LP>,StandardVisitor>> solver(argc,argv);
solver.ReadProblem(cell_tracking_parser::ParseProblemMotherMachine<Solver<FMC_CELL_TRACKING_MOTHER_MACHINE,LP_sat<LP>,StandardVisitor>>);
return solver.Solve();

}
