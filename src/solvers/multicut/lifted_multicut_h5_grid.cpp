
#include "solvers/multicut/multicut.h"
#include "visitors/standard_visitor.hxx"
using namespace LP_MP;
int main(int argc, char* argv[])

{
    using Rounder = LiftedKlRounder; 
    using SolverType = Solver<FMC_LIFTED_MULTICUT<Rounder>,LP,StandardTighteningVisitor,Rounder>;
    ProblemConstructorRoundingSolver<SolverType> solver(argc,argv);
    solver.ReadProblem(MulticutH5Input::ParseLiftedProblem<SolverType,true>);
    return solver.Solve();

}
