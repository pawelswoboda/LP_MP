#include "visitors/standard_visitor.hxx"
#include "multicut.h"

using namespace LP_MP;
using FMC = FMC_ODD_WHEEL_MULTICUT<MessageSendingType::SRMP>;
using SolverType =
ProblemConstructorRoundingSolver<Solver<FMC,LP,StandardTighteningVisitor>>;

int main(int argc, char* argv[])
{
   std::vector<std::string> options = {
      "mc_test_simple",
      "-i", " ",
      "--primalComputationInterval", "100",
      "--standardReparametrization", "anisotropic",
      "--roundingReparametrization", "damped_uniform",
      "--tightenReparametrization", "damped_uniform",
      "--tighten",
      "--tightenInterval", "100",
      "--tightenIteration", "2",
      "--tightenSlope", "0.05",
      "--tightenConstraintsPercentage", "0.1"
   };

   SolverType solver(options);
   auto& multicut_constructor = solver.template GetProblemConstructor<0>();
   const double cost = 0.5;
   multicut_constructor.AddUnaryFactor(0,1, cost);
   multicut_constructor.AddUnaryFactor(1,2, cost);
   multicut_constructor.AddUnaryFactor(2,3, cost);
   multicut_constructor.AddUnaryFactor(3,4, cost);
   solver.Solve();

   const bool cut1 = multicut_constructor.get_edge_label(0,1);
   const bool cut2 = multicut_constructor.get_edge_label(1,2);
   assert(!cut1);
   assert(!cut2);
   std::cout << "Passed" << std::endl;
}
