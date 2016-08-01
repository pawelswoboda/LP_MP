#ifndef LP_MP_SOLVER_HXX
#define LP_MP_SOLVER_HXX

#include "LP_MP.h"
#include "meta/meta.hpp"
#include "function_existence.hxx"
#include "template_utilities.hxx"
#include "spdlog/spdlog.h"
#include "tclap/CmdLine.h"

namespace LP_MP {

// class containing the LP, problem constructor list, input function
// takes care of logging
// binds together problem constructors and solver and organizes input/output
template<typename FMC>
class Solver {
   using ProblemDecompositionList = typename FMC::ProblemDecompositionList;
   using FactorMessageConnection = FMC;

   // initialize a tuple uniformly
   template <class T, class... ARGS>
      std::tuple<ARGS...> tupleMaker(meta::list<ARGS...>, T& t) { return std::make_tuple(ARGS(t)...); }

public:
   Solver(int argc, char** argv)
      : cmd_(std::string("Command line options for ") + FMC::name, ' ', "0.0.1"),
      lp_(LP()),
      //logger_();
      // do zrobienia: use perfect forwarding or std::piecewise_construct
      problemConstructor_(tupleMaker(ProblemDecompositionList{}, *this)),
      // build the standard command line arguments
      inputFileArg_("i","inputFile","file from which to read problem instance",true,"","file name",cmd_),
      outputFileArg_("o","outputFile","file to write solution",false,"","file name",cmd_)
      {
         // now initialize command line arguments
         //cmd_.parse(argc,argv); // problem is: 

         // initialize logger
         //std::vector<spdlog::sink_ptr> sinks;
         //if(protocolateConsole_) {
         //   sinks.push_back(std::make_shared<spdlog::sinks::stdout_sink_st>());
         //} 
         //if(protocolateFile_ != "") {
         //   sinks.push_back(std::make_shared<spdlog::sinks::simple_file_sink_st>(protocolateFile_.c_str(),true));
         //}
         //spdlog::logger logger_("", std::begin(sinks), std::end(sinks));
         //logger_.set_pattern("%v");
         //spdlog::register_logger(logger);

      }

   ~Solver() {}

   virtual int Solve() = 0; // must be implemented by solver taking some visitor.

   // needed, as more arguments could be passed to cmd_, and then we need to parse again
   void Init(int argc, char** argv)
   {
      cmd_.parse(argc,argv);
      try {  
         inputFile_ = inputFileArg_.getValue();
         outputFile_ = outputFileArg_.getValue();
      } catch (TCLAP::ArgException &e) {
         std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl; 
         exit(1);
      }
   }

   template<typename INPUT_FUNCTION, typename... ARGS>
   bool ReadProblem(INPUT_FUNCTION inputFct, ARGS... args)
   {
      const bool success = inputFct(inputFile_,*this,args...);

      assert(success);
      if(!success) throw std::runtime_error("could not parse problem file");

      //spdlog::get("logger")->info("loading file " + inputFile_ + " succeeded");
      return success;
   }

   //spdlog::logger& GetLogger() const { return logger_; }


   // invoke the corresponding functions of problem constructors
   LP_MP_FUNCTION_EXISTENCE_CLASS(HasCheckPrimalConsistency,CheckPrimalConsistency);
   template<INDEX PROBLEM_CONSTRUCTOR_NO>
   constexpr static bool
   CanCheckPrimalConsistency()
   {
      // do zrobienia: this is not nice. CanTighten should only be called with valid PROBLEM_CONSTRUCTOR_NO
      constexpr INDEX n = PROBLEM_CONSTRUCTOR_NO >= ProblemDecompositionList::size() ? 0 : PROBLEM_CONSTRUCTOR_NO;
      if(n < PROBLEM_CONSTRUCTOR_NO) return false;
      else return HasCheckPrimalConsistency<meta::at_c<ProblemDecompositionList,n>, bool, PrimalSolutionStorage::Element>();
      //static_assert(PROBLEM_CONSTRUCTOR_NO<ProblemDecompositionList::size(),"");
   }
   template<INDEX PROBLEM_CONSTRUCTOR_NO>
   typename std::enable_if<PROBLEM_CONSTRUCTOR_NO >= ProblemDecompositionList::size(),INDEX>::type
   CheckPrimalConsistency(PrimalSolutionStorage::Element primal) 
   {
      return true; 
   }
   template<INDEX PROBLEM_CONSTRUCTOR_NO>
   typename std::enable_if<PROBLEM_CONSTRUCTOR_NO < ProblemDecompositionList::size() && !CanCheckPrimalConsistency<PROBLEM_CONSTRUCTOR_NO>(),INDEX>::type
   CheckPrimalConsistency(PrimalSolutionStorage::Element primal)
   {
      return CheckPrimalConsistency<PROBLEM_CONSTRUCTOR_NO+1>(primal);
   }
   template<INDEX PROBLEM_CONSTRUCTOR_NO>
   typename std::enable_if<PROBLEM_CONSTRUCTOR_NO < ProblemDecompositionList::size() && CanCheckPrimalConsistency<PROBLEM_CONSTRUCTOR_NO>(),INDEX>::type
   CheckPrimalConsistency(PrimalSolutionStorage::Element primal)
   {
      if(std::get<PROBLEM_CONSTRUCTOR_NO>(problemConstructor_).CheckPrimalConsistency(primal)) {
         return CheckPrimalConsistency<PROBLEM_CONSTRUCTOR_NO+1>(primal);
      } else { 
         return false;
      }
   }
   bool CheckPrimalConsistency(PrimalSolutionStorage::Element primal) 
   {
      return CheckPrimalConsistency<0>(primal);
   }


   LP_MP_FUNCTION_EXISTENCE_CLASS(HasTighten,Tighten);
   template<INDEX PROBLEM_CONSTRUCTOR_NO>
   constexpr static bool
   CanTighten()
   {
      // do zrobienia: this is not nice. CanTighten should only be called with valid PROBLEM_CONSTRUCTOR_NO
      constexpr INDEX n = PROBLEM_CONSTRUCTOR_NO >= ProblemDecompositionList::size() ? 0 : PROBLEM_CONSTRUCTOR_NO;
      if(n < PROBLEM_CONSTRUCTOR_NO) return false;
      else return HasTighten<meta::at_c<ProblemDecompositionList,n>, INDEX, REAL, INDEX>();
      //static_assert(PROBLEM_CONSTRUCTOR_NO<ProblemDecompositionList::size(),"");
   }
   template<INDEX PROBLEM_CONSTRUCTOR_NO>
   typename std::enable_if<PROBLEM_CONSTRUCTOR_NO >= ProblemDecompositionList::size(),INDEX>::type
   Tighten(const REAL minDualIncrease, const INDEX maxConstraints) { return 0; }
   template<INDEX PROBLEM_CONSTRUCTOR_NO>
   typename std::enable_if<PROBLEM_CONSTRUCTOR_NO < ProblemDecompositionList::size() && !CanTighten<PROBLEM_CONSTRUCTOR_NO>(),INDEX>::type
   Tighten(const REAL minDualIncrease, const INDEX maxConstraints)
   {
      return Tighten<PROBLEM_CONSTRUCTOR_NO+1>(minDualIncrease,maxConstraints);
   }
   template<INDEX PROBLEM_CONSTRUCTOR_NO>
   typename std::enable_if<PROBLEM_CONSTRUCTOR_NO < ProblemDecompositionList::size() && CanTighten<PROBLEM_CONSTRUCTOR_NO>(),INDEX>::type
   Tighten(const REAL minDualIncrease, const INDEX maxConstraints) 
   {
      spdlog::get("logger")->info() << "Tighten for pc no " << PROBLEM_CONSTRUCTOR_NO;
      const INDEX noCuttingPlaneAdded = std::get<PROBLEM_CONSTRUCTOR_NO>(problemConstructor_).Tighten(minDualIncrease,maxConstraints);
      return noCuttingPlaneAdded + Tighten<PROBLEM_CONSTRUCTOR_NO+1>(minDualIncrease,maxConstraints);
   }
   // minDualIncrease says how small minimally must be the increase guaranteed by added constraints, while maxConstraints gives maximum number of constraints to add
   INDEX Tighten(const REAL minDualIncrease, const INDEX maxConstraints) 
   {
      return Tighten<0>(minDualIncrease, maxConstraints);
   }



   LP_MP_FUNCTION_EXISTENCE_CLASS(HasComputePrimal,ComputePrimal);
   template<INDEX PROBLEM_CONSTRUCTOR_NO>
   constexpr static bool
   CanComputePrimal()
   {
      // do zrobienia: this is not nice. CanComputePrimal should only be called with valid PROBLEM_CONSTRUCTOR_NO
      constexpr INDEX n = PROBLEM_CONSTRUCTOR_NO >= ProblemDecompositionList::size() ? 0 : PROBLEM_CONSTRUCTOR_NO;
      if(n < PROBLEM_CONSTRUCTOR_NO) return false;
      else return HasComputePrimal<meta::at_c<ProblemDecompositionList,n>, void, PrimalSolutionStorage::Element>();
      //static_assert(PROBLEM_CONSTRUCTOR_NO<ProblemDecompositionList::size(),"");
   }
   template<INDEX PROBLEM_CONSTRUCTOR_NO>
   typename std::enable_if<PROBLEM_CONSTRUCTOR_NO >= ProblemDecompositionList::size()>::type
   ComputePrimal(PrimalSolutionStorage::Element primal) { return; }
   template<INDEX PROBLEM_CONSTRUCTOR_NO>
   typename std::enable_if<PROBLEM_CONSTRUCTOR_NO < ProblemDecompositionList::size() && !CanComputePrimal<PROBLEM_CONSTRUCTOR_NO>()>::type
   ComputePrimal(PrimalSolutionStorage::Element primal)
   {
      return ComputePrimal<PROBLEM_CONSTRUCTOR_NO+1>(primal);
   }
   template<INDEX PROBLEM_CONSTRUCTOR_NO>
   typename std::enable_if<PROBLEM_CONSTRUCTOR_NO < ProblemDecompositionList::size() && CanComputePrimal<PROBLEM_CONSTRUCTOR_NO>()>::type
   ComputePrimal(PrimalSolutionStorage::Element primal)
   {
      spdlog::get("logger")->info() << "ComputePrimal for pc no " << PROBLEM_CONSTRUCTOR_NO;
      std::get<PROBLEM_CONSTRUCTOR_NO>(problemConstructor_).ComputePrimal(primal);
      return ComputePrimal<PROBLEM_CONSTRUCTOR_NO+1>(primal);
   }
   void ComputePrimal(PrimalSolutionStorage::Element primal)
   {
      ComputePrimal<0>(primal);
   }


   template<INDEX PROBLEM_CONSTRUCTOR_NO>
   meta::at_c<ProblemDecompositionList, PROBLEM_CONSTRUCTOR_NO>& GetProblemConstructor() 
   {
      return std::get<PROBLEM_CONSTRUCTOR_NO>(problemConstructor_);
   }

   LP& GetLP() { return lp_; }

protected:
   TCLAP::CmdLine cmd_;
   LP lp_;
   std::shared_ptr<spdlog::logger> logger_;
   tuple_from_list<ProblemDecompositionList> problemConstructor_;

   // command line arguments
   TCLAP::ValueArg<std::string> inputFileArg_;
   TCLAP::ValueArg<std::string> outputFileArg_;
   std::string inputFile_;
   std::string outputFile_;
};

// add visitor to solver. Problem is: problem constructors do not depend on visitor, hence shall be given only base solver.
template<typename FMC, typename VISITOR>
class VSolver : public Solver<FMC> {
public:
   using VisitorType = VISITOR;

   VSolver(int argc, char** argv)
      : Solver<FMC>(argc,argv),
      visitor_(Solver<FMC>::cmd_,*this)
   {
      Solver<FMC>::Init(argc,argv); // do zrobienia: not nice! but needed to parse command line, after arguments were added by visitor
   }

   int Solve()
   {
      int rt = Solver<FMC>::lp_.Solve(visitor_, [&](PrimalSolutionStorage::Element primal) -> bool { return this->CheckPrimalConsistency(primal); });
      // discrete tomo inspection
      /*
      auto& mrf = std::get<0>(Solver<FMC>::problemConstructor_);
      std::cout << "unary potentials:\n";
      for(INDEX i=0; i<mrf.GetNumberOfVariables();++i) {
         std::cout << i << ": ";
         auto* f = mrf.GetUnaryFactor(i);
         for(INDEX x=0; x<mrf.GetNumberOfLabels(i);++x) {
            std::cout << f->operator[](x) << ", ";
         }
         std::cout << "\n";
      }
      std::cout << "pairwise potentials:\n";
      for(INDEX c=0; c<mrf.GetNumberOfPairwiseFactors(); ++c) {
         auto* f = mrf.GetPairwiseFactor(c);
         auto ij = mrf.GetPairwiseVariables(c);
         const INDEX i=std::get<0>(ij);
         const INDEX j=std::get<1>(ij);
         std::cout << i << "," << j << ":\n";
         for(INDEX xi=0; xi<mrf.GetNumberOfLabels(i);++xi) {
            for(INDEX xj=0; xj<mrf.GetNumberOfLabels(j);++xj) {
               std::cout << f->operator[](xi + xj*mrf.GetNumberOfLabels(i)) << ", ";
            }
            std::cout << "\n";
         }
      }
      */
      return rt;
   }

private:
   VISITOR visitor_;
};


// Macro for generating main function 
// do zrobienia: get version number automatically from CMake 
// do zrobienia: version number for individual solvers?
#define LP_MP_CONSTRUCT_SOLVER_WITH_INPUT_AND_VISITOR(FMC,PARSE_PROBLEM_FUNCTION,VISITOR) \
using namespace LP_MP; \
int main(int argc, char* argv[]) \
{ \
   VSolver<FMC,VISITOR> solver(argc,argv); \
   solver.ReadProblem(PARSE_PROBLEM_FUNCTION); \
   return solver.Solve(); \
}


} // end namespace LP_MP

#endif // LP_MP_SOLVER_HXX

