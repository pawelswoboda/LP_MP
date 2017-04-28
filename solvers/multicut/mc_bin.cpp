#include "visitors/standard_visitor.hxx"
#include "multicut.h"

using namespace LP_MP;
using FMC = FMC_ODD_WHEEL_MULTICUT<MessageSendingType::SRMP>;
using SolverType =
ProblemConstructorRoundingSolver<Solver<FMC,LP,StandardTighteningVisitor>>;

//template<typename CONSTRUCTOR>
//void write_edge_result(const & CONSTRUCTOR constructor,
//        const & std::string outfile) {
//    // TODO iterate over edges and write to outfile
//    //const bool cut1 = multicut_constructor.get_edge_label(0,1);
//    //const bool cut2 = multicut_constructor.get_edge_label(1,2);
//}
   

int main()
{
    // FIXME hacky
    std::string infile;
    
    std::cout << "Path to input file" << std::endl;
    std::cin >> infile;
    
    std::vector<std::string> options = {
          "mc_bin",
          "-i", infile,
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
    solver.Solve();
    //auto & multicut_constructor = solver.template GetProblemConstructor<0>();
    //write_edge_result(multicut_constructor, "out.out")
}
