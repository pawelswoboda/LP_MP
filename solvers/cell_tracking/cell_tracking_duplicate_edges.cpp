
#include "cell_tracking.h"
#include "visitors/standard_visitor.hxx"
using namespace LP_MP;
int main(int argc, char* argv[])

{
MpRoundingSolver<Solver<FMC_CELL_TRACKING_DUPLICATE_EDGES,LP,StandardTighteningVisitor>> solver(argc,argv);
solver.ReadProblem(cell_tracking_parser_2d::ParseProblem<Solver<FMC_CELL_TRACKING_DUPLICATE_EDGES,LP,StandardTighteningVisitor>>);
return solver.Solve();

}
