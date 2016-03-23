#include "graph_matching.h"
#include "visitors/standard_visitor.hxx"
// do zrobienia: the above does not look nice. Refactor into one solver class
using FMC_INST = FMC_MP<PairwiseConstruction::Left>;
LP_MP_CONSTRUCT_SOLVER_WITH_INPUT_AND_VISITOR(FMC_INST, UAIInput::ParseProblem<FMC_INST>, StandardVisitor<ProblemDecomposition<FMC_INST>>);


