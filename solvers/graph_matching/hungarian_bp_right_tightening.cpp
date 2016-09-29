#include "graph_matching.h"
#include "visitors/standard_visitor.hxx"
using FMC_INST = FMC_HUNGARIAN_BP_T<PairwiseConstruction::Right>;
LP_MP_CONSTRUCT_SOLVER_WITH_INPUT_AND_VISITOR_MP_ROUNDING(FMC_INST, TorresaniEtAlInput::ParseProblemMCF<FMC_INST>, StandardTighteningVisitor);

