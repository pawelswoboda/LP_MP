#include "multicut.h"
#include "visitors/standard_visitor.hxx"
using namespace LP_MP;
using FMC = FMC_MULTICUT<MessageSendingType::SRMP>;
//using FMC = FMC_MULTICUT<MessageSending::SRMP>;
//LP_MP_CONSTRUCT_SOLVER_WITH_INPUT_AND_VISITOR_AND_SOLVER(FMC, MulticutOpenGmInput::ParseProblem<Solver<FMC>>, StandardTighteningVisitor, ProblemConstructorRoundingSolver<FMC>);
using SolverType = ProblemConstructorRoundingSolver<Solver<FMC,LP,StandardTighteningVisitor>>;
int main(int argc, char* argv[])
{
   SolverType solver(argc,argv);
   solver.ReadProblem(MulticutOpenGmInput::ParseProblem<Solver<FMC,LP,StandardTighteningVisitor>>);
   return solver.Solve();
}


