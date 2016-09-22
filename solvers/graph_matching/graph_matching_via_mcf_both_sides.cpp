#include "graph_matching.h"
#include "visitors/standard_visitor.hxx"
using namespace LP_MP;
using FMC_INST = FMC_MCF<PairwiseConstruction::BothSides>;
LP_MP_CONSTRUCT_SOLVER_WITH_INPUT_AND_VISITOR(FMC_INST, TorresaniEtAlInput::ParseProblemMCF<FMC_INST>, StandardVisitor);

