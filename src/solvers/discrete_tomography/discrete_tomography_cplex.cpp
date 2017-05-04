#include "solvers/discrete_tomography/discrete_tomography.h"
#include "lp_interface/lp_cplex.hxx"
#include "visitors/standard_visitor.hxx"
using namespace LP_MP;
LP_MP_CONSTRUCT_SOLVER_WITH_INPUT_AND_LPVISITOR(FMC_DT, DiscreteTomographyTextInput::ParseProblem<FMC_DT>, StandardTighteningVisitor, LpInterfaceCplex);
/*
int main(int argc, char* argv[])                                      
{                                                                     
  //VSolver<FMC_DT,StandardTighteningVisitor<Solver<FMC_DT> >> solver(argc,argv);                             
  
  LpRoundingSolver<FMC_DT,LpInterfaceGurobi> solver(argc,argv);  
  solver.ReadProblem(DiscreteTomographyTextInput::ParseProblem);
  return solver.Solve();                                              
}
*/
