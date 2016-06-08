#ifndef LP_MP_SOLVER_HXX
#define LP_MP_SOLVER_HXX

#include "LP_MP.h"
#include "spdlog/spdlog.h"
#include "tclap/CmdLine.h"

namespace LP_MP {

// class containing the FactorList, MessageList, ProblemConstructorList, Visitor, InputFunction
// contains dual bound, primal solution
// takes care of logging
// binds together constructors, solver and organizes input/output
template<typename FMC, typename INPUT_FUNCTION, typename VISITOR>
class Solver {
public:
   using InputFunction = INPUT_FUNTION; // do zrobienia: this is not so easy to do as a predefined template
   using VisitorType = VISITOR;

   Solver(int argc, char** argv)
      : cmd_(std::string("Command line options for ") + FMC::name, ' ', "0.0.1"),
      lp_(LP()),
      visitor_(cmd),
      // do zrobienia: use perfect forwarding or std::piecewise_construct
      problem_constructor_(tupleMaker(FMC::ProblemDecompositionList{}, *this)),
      // first we build the standard command line arguments
      inputFileArg_("i","inputFile","file from which to read problem instance",true,"","file name",cmd_),
      outputFileArg_("o","outputFile","file to write solution",false,"","file name",cmd_)
      {
         // now initialize command line arguments
         cmd.parse(argc,argv);

         // initialize logger
         std::vector<spdlog::sink_ptr> sinks;
         if(protocolateConsole_) {
            sinks.push_back(std::make_shared<spdlog::sinks::stdout_sink_st>());
         } 
         if(protocolateFile_ != "") {
            sinks.push_back(std::make_shared<spdlog::sinks::simple_file_sink_st>(protocolateFile_.c_str(),true));
         }
         spdlog::logger logger_("", std::begin(sinks), std::end(sinks));
         logger_.set_pattern("%v");
         //spdlog::register_logger(logger);

      }

   ~ProblemDecomposition() 
   { 
      //if(lp_ != nullptr) { delete lp_; }
   }

   template<typename INPUT_FUNCTION, typename... ARGS>
   bool ReadProblem(INPUT_FUNCTION inputFct, ARGS... args)
   {

   }

   spdlog::logger& GetLogger() const { return logger_; }

private:
   TCLAP::CmdLine cmd_;
   LP lp_;
   spdlog::logger logger_;
   VISITOR visitor_;
   REAL dualBound_;
   PrimalSolutionStorage bestPrimal_, curPrimal_;
   //tuple_from_list<ProblemDecompositionList> problemConstructor_;
   ProblemConstructor problemConstructor_;

   // command line arguments
   TCLAP::ValueArg<std::string> inputFileArg_;
   TCLAP::ValueArg<std::string> outputFileArg_;
   std::string inputFile_;
   std::string outputFile_;

}


} // end namespace LP_MP

#endif // LP_MP_SOLVER_HXX

