#include "discrete_tomography.h"
#include "lp_interface/lp_localsolver.hxx"
#include "visitors/standard_visitor.hxx"
using namespace LP_MP;
LP_MP_CONSTRUCT_SOLVER_WITH_INPUT_AND_LPVISITOR(FMC_DT, DiscreteTomographyTextInput::ParseProblem, StandardTighteningVisitor, LpInterfaceLocalSolver);

