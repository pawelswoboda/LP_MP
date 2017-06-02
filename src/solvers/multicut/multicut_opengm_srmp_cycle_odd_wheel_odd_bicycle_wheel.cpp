#include "solvers/multicut/multicut.h"
#include "visitors/standard_visitor.hxx"

using namespace LP_MP;
int main(int argc, char* argv[])
{
    typedef KlRounder Rounder;
    typedef Solver<FMC_ODD_BICYCLE_WHEEL_MULTICUT<Rounder>,LP,StandardTighteningVisitor,Rounder> SolverBase;
    ProblemConstructorRoundingSolver<SolverBase> solver(argc,argv);
    solver.ReadProblem(MulticutOpenGmInput::ParseProblem<SolverBase>);
    return solver.Solve();
}
