
#include "graphical_model.h"
#include "visitors/standard_visitor.hxx"
int main(int argc, char* argv[])

{
MpRoundingSolver<Solver<FMC_SRMP,LP,StandardTighteningVisitor>> solver(argc,argv);
solver.ReadProblem(UaiMrfInput::ParseProblem<Solver<FMC_SRMP,LP,StandardTighteningVisitor>>);
return solver.Solve();

}
