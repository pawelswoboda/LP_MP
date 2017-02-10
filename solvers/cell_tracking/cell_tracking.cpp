#include "cell_tracking.h"
#include "solver.hxx"
#include "visitors/standard_visitor.hxx"
using FMC_INST = LP_MP::FMC_CELL_TRACKING_MOTHER_MACHINE;
LP_MP_CONSTRUCT_SOLVER_WITH_INPUT_AND_VISITOR_MP_ROUNDING(FMC_INST, LP_MP::cell_tracking_parser::ParseProblemMotherMachine<FMC_INST>, StandardVisitor);

