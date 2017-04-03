#include "solver.hxx"
#include "visitors/standard_visitor.hxx"
#include "moral_lineage_tracing.h"

using namespace LP_MP;

int main(int argc, char** argv)
{
	// specify FMC struct, LP solver type, visitor type
	using SolverType = Solver<FMC_MLT,LP,StandardTighteningVisitor>;
	SolverType s(argc,argv);
	s.ReadProblem(mlt_input::ParseProblem<SolverType>);
	s.Solve();
}