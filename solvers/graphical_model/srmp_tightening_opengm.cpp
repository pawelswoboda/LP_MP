#include "graphical_model.h"
#include "visitors/standard_visitor.hxx"
using FMC_INST = FMC_SRMP_T;
LP_MP_CONSTRUCT_SOLVER_WITH_INPUT_AND_VISITOR(FMC_INST, HDF5Input::ParseProblem<FMC_INST>, StandardTighteningVisitor);
