#include "solvers/multicut/multicut.h"
#include "visitors/standard_visitor.hxx"
using namespace LP_MP;

int main(int argc, char* argv[])
{
    using SolverType = Solver<FMC_LIFTED_MULTICUT<KlRounder,LiftedKlRounder>,LP,StandardTighteningVisitor>;
    ProblemConstructorRoundingSolver<SolverType> solver(argc,argv);
    solver.ReadProblem(MulticutTextInput::ParseLiftedProblem<SolverType>);
    return solver.Solve();
}
