#include "graph_matching.h"
#include "visitors/standard_visitor.hxx"
using FMC_INST = FMC_MP<PairwiseConstruction::BothSides>;
LP_MP_CONSTRUCT_SOLVER_WITH_INPUT_AND_VISITOR(FMC_INST, TorresaniEtAlInput::ParseProblemMP<FMC_INST>, StandardVisitor);
