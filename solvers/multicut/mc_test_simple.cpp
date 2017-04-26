#include "visitors/standard_visitor.hxx"
#include "multicut.h"

using namespace LP_MP;
using FMC = FMC_ODD_WHEEL_MULTICUT<MessageSendingType::SRMP>;
using SolverType =
ProblemConstructorRoundingSolver<Solver<FMC,LP,StandardTighteningVisitor>>;

int main(int argc, char* argv[])
{
   std::vector<std::string> options = {
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
   solver.Solve();

   const bool cut = multicut_constructor.get_edge_label(0,1);
   assert(!cut);
   std::cout << "Passed" << std::endl;
}
