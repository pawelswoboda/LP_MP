#include "solvers/multicut/multicut.h"
#include "visitors/standard_visitor.hxx"

using namespace LP_MP;
int main(int argc, char* argv[])
{
    using SolverType = Solver<FMC_MULTIWAY_CUT<KlRounder>,LP,StandardTighteningVisitor>;
    ProblemConstructorRoundingSolver<SolverType> solver(argc,argv);
    solver.ReadProblem(MulticutOpenGmInput::ParsePottsProblem<SolverType>);
    return solver.Solve();
}
