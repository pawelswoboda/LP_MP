
#include "multicut.h"
#include "visitors/standard_visitor.hxx"
using namespace LP_MP;
int main(int argc, char* argv[])

{
CombinedMPProblemConstructorRoundingSolver<Solver<FMC_ODD_WHEEL_MULTICUT<MessageSendingType::SRMP>,LP,StandardTighteningVisitor>> solver(argc,argv);
solver.ReadProblem(MulticutOpenGmInput::ParseProblem<Solver<FMC_ODD_WHEEL_MULTICUT<MessageSendingType::SRMP>,LP,StandardTighteningVisitor>>);
return solver.Solve();

}
