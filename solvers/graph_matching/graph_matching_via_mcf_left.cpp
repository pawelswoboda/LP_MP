#include "graph_matching.h"
#include "visitors/standard_visitor.hxx"
using FMC_INST = FMC_MCF<PairwiseConstruction::Left>;
LP_MP_CONSTRUCT_SOLVER_WITH_INPUT_AND_VISITOR(FMC_INST, TorresaniEtAlInput::ParseProblemMCF<FMC_INST>, StandardVisitor);
