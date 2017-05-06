#include "solvers/multicut/multicut.h"
#include "visitors/standard_visitor.hxx"

using namespace LP_MP;
int main(int argc, char* argv[])
{
    typedef KlRounder Rounder;
    typedef Solver<FMC_ODD_WHEEL_MULTICUT<MessageSendingType::SRMP,Rounder>,LP,StandardTighteningVisitor,Rounder> SolverBase;
    ProblemConstructorRoundingSolver<SolverBase> solver(argc,argv);
    solver.ReadProblem(MulticutOpenGmInput::ParseProblem<SolverBase>);
    return solver.Solve();
}
