#include "multicut.h"
#include "visitors/standard_visitor.hxx"
using namespace LP_MP;
using FMC = FMC_MULTICUT<MessageSendingType::SRMP>;
//using FMC = FMC_MULTICUT<MessageSending::SRMP>;
//LP_MP_CONSTRUCT_SOLVER_WITH_INPUT_AND_VISITOR_AND_SOLVER(FMC, MulticutOpenGmInput::ParseProblem<Solver<FMC>>, StandardTighteningVisitor, ProblemConstructorRoundingSolver<FMC>);
//#ifdef LP_PM_PARALLEL
//using SolverType = ProblemConstructorRoundingSolver<Solver<FMC,LP_concurrent<LP>,StandardTighteningVisitor>>;
//#else
using SolverType = ProblemConstructorRoundingSolver<Solver<FMC,LP,StandardTighteningVisitor>>;
//#endif
int main(int argc, char* argv[])
{
   SolverType solver(argc,argv);
//#ifdef LP_PM_PARALLEL
//   solver.ReadProblem(MulticutOpenGmInput::ParseProblem<Solver<FMC,LP_concurrent<LP>,StandardTighteningVisitor>>);
//#else
   solver.ReadProblem(MulticutOpenGmInput::ParseProblem<Solver<FMC,LP,StandardTighteningVisitor>>);
//#endif
   return solver.Solve();
}


