#include "discrete_tomography.h"
#include "visitors/standard_visitor.hxx"
using namespace LP_MP;
LP_MP_CONSTRUCT_SOLVER_WITH_INPUT_AND_VISITOR(FMC_DT, DiscreteTomographyTextInput::ParseProblem, StandardVisitor<Solver<FMC_DT> >);
