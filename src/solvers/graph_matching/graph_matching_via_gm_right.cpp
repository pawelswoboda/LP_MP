
#include "solvers/graph_matching/graph_matching.h"
#include "visitors/standard_visitor.hxx"
int main(int argc, char* argv[])

{
MpRoundingSolver<Solver<FMC_GM<PairwiseConstruction::Right>,LP,StandardTighteningVisitor>> solver(argc,argv);
solver.ReadProblem(TorresaniEtAlInput::ParseProblemGM<Solver<FMC_GM<PairwiseConstruction::Right>,LP,StandardTighteningVisitor>>);
return solver.Solve();

}
