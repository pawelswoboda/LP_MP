#include "multicut.h"
#include "visitors/standard_visitor.hxx"
using namespace LP_MP;
using FMC = FMC_LIFTED_MULTICUT;
LP_MP_CONSTRUCT_SOLVER_WITH_INPUT_AND_VISITOR(FMC, MulticutTextInput::ParseLiftedProblem<FMC>, StandardTighteningVisitor);


