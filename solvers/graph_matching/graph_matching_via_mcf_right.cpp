
#include "graph_matching.h"
#include "visitors/standard_visitor.hxx"
int main(int argc, char* argv[])

{
MpRoundingSolver<Solver<FMC_MCF<PairwiseConstruction::Right>,LP,StandardTighteningVisitor>> solver(argc,argv);
solver.ReadProblem(TorresaniEtAlInput::ParseProblemMCF<Solver<FMC_MCF<PairwiseConstruction::Right>,LP,StandardTighteningVisitor>>);
return solver.Solve();

}
