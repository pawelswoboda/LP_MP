#include "graph_matching.h"
#include "visitors/ascii_plot_visitor.hxx"
using FMC_INST = FMC_MCF<PairwiseConstruction::Left>;
LP_MP_CONSTRUCT_SOLVER_WITH_INPUT_AND_VISITOR(FMC_INST, UAIInput::ParseProblem<FMC_INST>, AsciiPlotVisitor<ProblemDecomposition<FMC_INST>>);


