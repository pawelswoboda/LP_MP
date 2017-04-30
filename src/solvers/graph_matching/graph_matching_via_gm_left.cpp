
#include "graph_matching.h"
#include "visitors/standard_visitor.hxx"
int main(int argc, char* argv[])

{
MpRoundingSolver<Solver<FMC_GM<PairwiseConstruction::Left>,LP,StandardTighteningVisitor>> solver(argc,argv);
solver.ReadProblem(TorresaniEtAlInput::ParseProblemGM<Solver<FMC_GM<PairwiseConstruction::Left>,LP,StandardTighteningVisitor>>);
return solver.Solve();

}
