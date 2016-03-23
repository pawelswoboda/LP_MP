#include "multicut.h"
#include "visitors/standard_visitor.hxx"
using namespace LP_MP;
using FMC = FMC_MULTICUT<MessageSending::SRMP>;
LP_MP_CONSTRUCT_SOLVER_WITH_INPUT_AND_VISITOR(FMC, MulticutTextInput::ParseProblem<FMC>, StandardTighteningVisitor<ProblemDecomposition<FMC>>);

