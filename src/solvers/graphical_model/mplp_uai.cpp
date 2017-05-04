
#include "solvers/graphical_model/graphical_model.h"
#include "visitors/standard_visitor.hxx"
int main(int argc, char* argv[])

{
MpRoundingSolver<Solver<FMC_MPLP,LP,StandardTighteningVisitor>> solver(argc,argv);
solver.ReadProblem(UaiMrfInput::ParseProblem<Solver<FMC_MPLP,LP,StandardTighteningVisitor>>);
return solver.Solve();

}
