#include "multicut.h"
#include "visitors/standard_visitor.hxx"
using namespace LP_MP;
using FMC = FMC_MULTICUT<MessageSendingType::SRMP>;
//using FMC = FMC_MULTICUT<MessageSending::SRMP>;
//LP_MP_CONSTRUCT_SOLVER_WITH_INPUT_AND_VISITOR_AND_SOLVER(FMC, MulticutOpenGmInput::ParseProblem<Solver<FMC>>, StandardTighteningVisitor, ProblemConstructorRoundingSolver<FMC>);
using SolverType = ProblemConstructorRoundingSolver<FMC>;
int main(int argc, char* argv[])
{
   VisitorSolver<SolverType,StandardTighteningVisitor> solver(argc,argv);
   solver.ReadProblem(MulticutOpenGmInput::ParseProblem<FMC>);
   return solver.Solve();
}


