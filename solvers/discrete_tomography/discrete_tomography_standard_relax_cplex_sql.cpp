#include "discrete_tomography.h"
#include "lp_interface/lp_cplex.hxx"
#include "visitors/sqlite_visitor.hxx"
using namespace LP_MP;
using VisitorType = SqliteVisitor<StandardTighteningVisitor>;
LP_MP_CONSTRUCT_SOLVER_WITH_INPUT_AND_VISITOR_LP(FMC_DT_NAIVE, DiscreteTomographyTextInput::ParseProblem<FMC_DT_NAIVE>, VisitorType, LpInterfaceCplex);

