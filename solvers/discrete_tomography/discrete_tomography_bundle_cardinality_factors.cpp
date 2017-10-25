#include "discrete_tomography.h"
#include "LP_conic_bundle.hxx"
#include "visitors/standard_visitor.hxx"
using namespace LP_MP;
using SolverType = Solver<FMC_DT_NAIVE,LP_conic_bundle,StandardTighteningVisitor>;
//LP_MP_CONSTRUCT_SOLVER_WITH_INPUT_AND_VISITOR(FMC_DT, DiscreteTomographyTextInput::ParseProblemDD<SolverType>, StandardTighteningVisitor);

int main(int argc, char* argv[])

{
SolverType solver(argc,argv);
solver.ReadProblem(DiscreteTomographyTextInput::ParseProblemDD<SolverType>);
return solver.Solve();

}
