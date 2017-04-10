#include "discrete_tomography.h"
#include "visitors/standard_visitor.hxx"
using namespace LP_MP;
using SolverType = Solver<FMC_DT,LP,StandardTighteningVisitor>;
//LP_MP_CONSTRUCT_SOLVER_WITH_INPUT_AND_VISITOR_MP_ROUNDING(FMC_DT, DiscreteTomographyTextInput::ParseProblem<FMC_DT>, StandardTighteningVisitor);
LP_MP_CONSTRUCT_SOLVER_WITH_INPUT_AND_VISITOR(FMC_DT, DiscreteTomographyTextInput::ParseProblem<SolverType>, StandardTighteningVisitor);
