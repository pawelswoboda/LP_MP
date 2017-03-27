#include "cell_tracking.h"
#include "solver.hxx"
#include "visitors/standard_visitor.hxx"
using FMC_INST = LP_MP::FMC_CONSERVATION_TRACKING;
LP_MP_CONSTRUCT_SOLVER_WITH_INPUT_AND_VISITOR(FMC_INST, LP_MP::conservation_tracking_parser::ParseProblem<FMC_INST>, StandardVisitor);

