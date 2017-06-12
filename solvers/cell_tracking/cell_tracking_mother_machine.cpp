
#include "cell_tracking.h"
#include "visitors/standard_visitor.hxx"
using namespace LP_MP;
int main(int argc, char* argv[])

{
Solver<FMC_CELL_TRACKING_MOTHER_MACHINE,LP_sat<LP>,StandardTighteningVisitor> solver(argc,argv);
solver.ReadProblem(cell_tracking_parser_mother_machine::ParseProblemMotherMachine<Solver<FMC_CELL_TRACKING_MOTHER_MACHINE,LP_sat<LP>,StandardTighteningVisitor>>);
return solver.Solve();

}
