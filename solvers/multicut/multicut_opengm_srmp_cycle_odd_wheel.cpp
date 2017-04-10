#include "multicut.h"
#include "visitors/standard_visitor.hxx"
using namespace LP_MP;
using FMC = FMC_ODD_WHEEL_MULTICUT<MessageSendingType::SRMP>;
//using FMC = FMC_ODD_WHEEL_MULTICUT<MessageSending::SRMP>;
//LP_MP_CONSTRUCT_SOLVER_WITH_INPUT_AND_VISITOR_AND_SOLVER(FMC, MulticutOpenGmInput::ParseProblem<Solver<FMC>>, StandardTighteningVisitor, ProblemConstructorRoundingSolver<FMC>);
using SolverType = ProblemConstructorRoundingSolver<Solver<FMC,LP,StandardTighteningVisitor>>;
int main(int argc, char* argv[])
{
  
   /*
      multicut_triplet_odd_3_wheel_message_012::print_matching();
      std::cout << "\n";
      multicut_triplet_odd_3_wheel_message_013::print_matching(); 
      std::cout << "\n";
      multicut_triplet_odd_3_wheel_message_023::print_matching();
      std::cout << "\n";
      multicut_triplet_odd_3_wheel_message_123::print_matching();
      std::cout << "\n";
    */ 

   SolverType solver(argc,argv);
   solver.ReadProblem(MulticutOpenGmInput::ParseProblem<Solver<FMC,LP,StandardTighteningVisitor>>);
   return solver.Solve();
}


