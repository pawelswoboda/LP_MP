
#include "cell_tracking.h"
#include "visitors/standard_visitor.hxx"
using namespace LP_MP;
int main(int argc, char* argv[])

{
MpRoundingSolver<Solver<FMC_CELL_TRACKING_WITH_DIVISION_DISTANCE,LP_sat<LP>,StandardTighteningVisitor>> solver(argc,argv);
solver.ReadProblem(cell_tracking_parser_2d::parse_problem_with_division_distance<Solver<FMC_CELL_TRACKING_WITH_DIVISION_DISTANCE,LP_sat<LP>,StandardTighteningVisitor>>);
return solver.Solve();

}
