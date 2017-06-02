#include "solvers/discrete_tomography/discrete_tomography.h"
#include "lp_interface/lp_gurobi.hxx"
#include "visitors/standard_visitor.hxx"
using namespace LP_MP;
LP_MP_CONSTRUCT_SOLVER_WITH_INPUT_AND_VISITOR_LP(FMC_DT_COMBINED, DiscreteTomographyTextInput::ParseProblemCombined<FMC_DT_COMBINED>, StandardTighteningVisitor, LpInterfaceGurobi);
