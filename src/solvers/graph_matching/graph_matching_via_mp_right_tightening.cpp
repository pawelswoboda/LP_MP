
#include "graph_matching.h"
#include "visitors/standard_visitor.hxx"
int main(int argc, char* argv[])

{
MpRoundingSolver<Solver<FMC_MP_T<PairwiseConstruction::Right>,LP,StandardTighteningVisitor>> solver(argc,argv);
solver.ReadProblem(TorresaniEtAlInput::ParseProblemMP<Solver<FMC_MP_T<PairwiseConstruction::Right>,LP,StandardTighteningVisitor>>);
return solver.Solve();

}
