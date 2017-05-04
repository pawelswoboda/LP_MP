
#include "solvers/graph_matching/graph_matching.h"
#include "visitors/standard_visitor.hxx"
int main(int argc, char* argv[])

{
MpRoundingSolver<Solver<FMC_MP<PairwiseConstruction::Left>,LP,StandardTighteningVisitor>> solver(argc,argv);
solver.ReadProblem(TorresaniEtAlInput::ParseProblemMP<Solver<FMC_MP<PairwiseConstruction::Left>,LP,StandardTighteningVisitor>>);
return solver.Solve();

}
