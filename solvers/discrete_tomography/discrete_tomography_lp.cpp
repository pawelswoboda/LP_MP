#include "discrete_tomography.h"
#include "lp_interface/lp_gurobi.hxx"
#include "visitors/standard_visitor.hxx"

using namespace LP_MP;

int main(int argc, char* argv[])                                      
{                                                                     
  //VSolver<FMC_DT,StandardTighteningVisitor<Solver<FMC_DT> >> solver(argc,argv);                             
  
  SolverLP<FMC_DT,LpInterfaceGurobi> solver(argc,argv);  
  solver.ReadProblem(DiscreteTomographyTextInput::ParseProblem);
  return solver.Solve();                                              
}
