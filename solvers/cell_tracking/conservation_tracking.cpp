
#include "cell_tracking.h"
#include "visitors/standard_visitor.hxx"
using namespace LP_MP;
int main(int argc, char* argv[])

{
Solver<FMC_CONSERVATION_TRACKING,LP_sat<LP>,StandardTighteningVisitor> solver(argc,argv);
solver.ReadProblem(conservation_tracking_parser::ParseProblem<Solver<FMC_CONSERVATION_TRACKING,LP_sat<LP>,StandardTighteningVisitor>>);
return solver.Solve();

}
