#include "discrete_tomography.h"
#include "visitors/standard_visitor.hxx"
using namespace LP_MP;
//using input_fct = DiscreteTomographyTextInput::ParseProblem;
LP_MP_CONSTRUCT_SOLVER_WITH_INPUT_AND_VISITOR_MP_ROUNDING(FMC_DT, DiscreteTomographyTextInput::ParseProblem, StandardTighteningVisitor);
