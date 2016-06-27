/*
#ifndef LP_MP_PROBLEM_DECOMPOSITION_HXX
#define LP_MP_PROBLEM_DECOMPOSITION_HXX

#include "LP_MP.h"
#include "tolerance.hxx"
#include "function_existence.hxx"
#include "meta/meta.hpp"
#include "tclap/CmdLine.h"

#include <fstream>
#include <sstream>
#include <string>
#include <regex>

#include "visitors/standard_visitor.hxx"

// construct problem from factor-message network with problem decomposition by invoking the correct problem decomposition routines of the respective classes.
namespace LP_MP {


// do zrobienia: better name for ProblemDecomposition?
// do zrobienia: templatize for input reading function
template<typename FMC>
class ProblemDecomposition {

   // initialize a tuple uniformly
   template <class T, class... ARGS>
      std::tuple<ARGS...> tupleMaker(meta::list<ARGS...>, T& t) { return std::make_tuple(ARGS(t)...); }
public:
   using ProblemDecompositionList = typename FMC::ProblemDecompositionList;
   using FactorMessageConnection = FMC;

   ProblemDecomposition(TCLAP::CmdLine& cmd)
   : lp_(new LP()),
      // first we build the standard command line arguments
      inputFileArg_("i","inputFile","file from which to read problem instance",true,"","file name",cmd),
      outputFileArg_("o","outputFile","file to write solution",false,"","file name",cmd),
     // do zrobienia: use perfect forwarding or std::piecewise_construct
     problem_constructor_(tupleMaker(ProblemDecompositionList{}, *this))
   {}
   ~ProblemDecomposition() 
   { 
      if(lp_ != nullptr) { delete lp_; }
   }

   template<typename VISITOR>
   SIGNED_INDEX Run(VISITOR& visitor)
   {
      int rc = GetLP()->Solve(visitor);
      if(rc == 0 && outputFile_ != "") {
         WritePrimal(outputFile_);
      }
      return rc;
   }

   template<typename PARSE_PROBLEM_FUNCTION>
   void ParseProblem(const PARSE_PROBLEM_FUNCTION& f)
   {
      ReadCommandLine();

      const bool success = f(inputFile_,*this);
      assert(success);
      if(!success) throw std::runtime_error("could not parse problem file");

      //spdlog::get("logger")->info("loading file " + inputFile_ + " succeeded");
      return;
   }


   // do zrobienia: does not work anymore with custom input grammars.
   template<INDEX PROBLEM_CONSTRUCTOR_NO, typename PRIMAL_SOLUTION_ITERATOR, typename PC_LIST = ProblemDecompositionList>
   typename std::enable_if<PROBLEM_CONSTRUCTOR_NO == PC_LIST::size()>::type
   WritePrimal(PRIMAL_SOLUTION_ITERATOR primalIt, std::ofstream& fs) {}
   template<INDEX PROBLEM_CONSTRUCTOR_NO, typename PRIMAL_SOLUTION_ITERATOR, typename PC_LIST = ProblemDecompositionList>
   typename std::enable_if<PROBLEM_CONSTRUCTOR_NO < PC_LIST::size()>::type
   WritePrimal(PRIMAL_SOLUTION_ITERATOR primalIt, std::ofstream& fs)
   {
      const INDEX factorIndexBegin = factorIndex_[PROBLEM_CONSTRUCTOR_NO];
      INDEX factorIndexEnd;
      if(PROBLEM_CONSTRUCTOR_NO != PC_LIST::size()-1) {
         factorIndexEnd = factorIndex_[PROBLEM_CONSTRUCTOR_NO+1];
      } else {
         factorIndexEnd = lp_->GetNumberOfFactors();
      }
      fs << "problem " << PROBLEM_CONSTRUCTOR_NO << "\n";
      lp_->WritePrimal(factorIndexBegin, factorIndexEnd, primalIt, fs);
      WritePrimal<PROBLEM_CONSTRUCTOR_NO+1>(primalIt+ factorIndexEnd-factorIndexBegin, fs);
   }
   void WritePrimal(const std::string& filename)
   {
      assert(false); // do zrobienia: this does not work, if ParseProblem was called
      std::ofstream fs{filename};
      if(!fs.is_open()) { 
         throw std::runtime_error("Could not open output file " + filename);
      }
      //auto primalSolution = lp_->GetBestPrimal();
      //WritePrimal<0>(primalSolution.begin(),fs);
      fs.close();
   }

   template<INDEX PROBLEM_CONSTRUCTOR_NO>
   meta::at_c<ProblemDecompositionList, PROBLEM_CONSTRUCTOR_NO>& GetProblemConstructor() 
   {
      return std::get<PROBLEM_CONSTRUCTOR_NO>(problem_constructor_);
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
      const INDEX noCuttingPlaneAdded = std::get<PROBLEM_CONSTRUCTOR_NO>(problem_constructor_).Tighten(minDualIncrease,maxConstraints);
      return noCuttingPlaneAdded + Tighten<PROBLEM_CONSTRUCTOR_NO+1>(minDualIncrease,maxConstraints);
   }
   INDEX Tighten(const REAL minDualIncrease, const INDEX maxConstraints) // minDualIncrease says how small minimally must be the increase guaranteed by added constraints, while maxConstraints gives maximum number of constraints to add
   {
      return Tighten<0>(minDualIncrease, maxConstraints);
   }

   LP_MP_FUNCTION_EXISTENCE_CLASS(HasCheckPrimalConsistency,CheckPrimalConsistency);
   template<INDEX PROBLEM_CONSTRUCTOR_NO>
   constexpr static bool
   CanCheckPrimalConsistency()
   {
      // do zrobienia: this is not nice. CanTighten should only be called with valid PROBLEM_CONSTRUCTOR_NO
      constexpr INDEX n = PROBLEM_CONSTRUCTOR_NO >= ProblemDecompositionList::size() ? 0 : PROBLEM_CONSTRUCTOR_NO;
      if(n < PROBLEM_CONSTRUCTOR_NO) return false;
      else return HasCheckPrimalConsistency<meta::at_c<ProblemDecompositionList,n>, bool, std::vector<bool>>();
      //static_assert(PROBLEM_CONSTRUCTOR_NO<ProblemDecompositionList::size(),"");
   }
   template<INDEX PROBLEM_CONSTRUCTOR_NO>
   typename std::enable_if<PROBLEM_CONSTRUCTOR_NO >= ProblemDecompositionList::size(),INDEX>::type
   CheckPrimalConsistency(const std::vector<bool>& primal) { return true; }
   template<INDEX PROBLEM_CONSTRUCTOR_NO>
   typename std::enable_if<PROBLEM_CONSTRUCTOR_NO < ProblemDecompositionList::size() && !CanCheckPrimalConsistency<PROBLEM_CONSTRUCTOR_NO>(),INDEX>::type
   CheckPrimalConsistency(const std::vector<bool>& primal)
   {
      return CheckPrimalConsistency<PROBLEM_CONSTRUCTOR_NO+1>(primal);
   }
   template<INDEX PROBLEM_CONSTRUCTOR_NO>
   typename std::enable_if<PROBLEM_CONSTRUCTOR_NO < ProblemDecompositionList::size() && CanCheckPrimalConsistency<PROBLEM_CONSTRUCTOR_NO>(),INDEX>::type
   CheckPrimalConsistency(std::vector<bool>& primal)
   {
      if(std::get<PROBLEM_CONSTRUCTOR_NO>(problem_constructor_).CheckPrimalConsistency(primal)) {
         return CheckPrimalConsistency<PROBLEM_CONSTRUCTOR_NO+1>(primal);
      } 
      return false;
   }
   bool CheckPrimalConsistency(const std::vector<bool>& primal) // minDualIncrease says how small minimally must be the increase guaranteed by added constraints, while maxConstraints gives maximum number of constraints to add
   {
      return CheckPrimalConsistency<0>(primal);
   }



   // do zrobienia: remove PC_LIST from templates and directly use ProblemDecompositionList as above
   template<INDEX PROBLEM_CONSTRUCTOR_NO, typename PC_LIST = ProblemDecompositionList>
   typename std::enable_if<PROBLEM_CONSTRUCTOR_NO == PC_LIST::size()>::type
   Construct() {}
   template<INDEX PROBLEM_CONSTRUCTOR_NO, typename PC_LIST = ProblemDecompositionList>
   typename std::enable_if<PROBLEM_CONSTRUCTOR_NO < PC_LIST::size()>::type
   Construct()
   {
      factorIndex_[PROBLEM_CONSTRUCTOR_NO] = lp_->GetNumberOfFactors();
      GetProblemConstructor<PROBLEM_CONSTRUCTOR_NO>().Construct(*this);
      return Construct<PROBLEM_CONSTRUCTOR_NO+1>();
   }

   // call Construct for each problem decomposition
   void Construct()
   {
      return Construct<0>();
   }

   LP* GetLP() const { return lp_; }
   
   void ReadCommandLine() 
   {
      try {  
         inputFile_ = inputFileArg_.getValue();
         outputFile_ = outputFileArg_.getValue();
      } catch (TCLAP::ArgException &e) {
         std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl; 
         exit(1);
      }
   }

   const std::string& GetInputFileName() const { return inputFile_; }

private:

   LP* lp_;
   tuple_from_list<ProblemDecompositionList> problem_constructor_;
   // do zrobienia: this does no work currently: for ParseProblem we call the Construct function explicitly, but then factorIndex is not called
   std::array<INDEX,meta::size<ProblemDecompositionList>::value> factorIndex_; // gives for each problem constructor the offset, from which it started to add factors to the LP class

   // command line arguments
   TCLAP::ValueArg<std::string> inputFileArg_;
   TCLAP::ValueArg<std::string> outputFileArg_;
   std::string inputFile_;
   std::string outputFile_;
};


// Macro for generating main function reading in file and then solving problem, do zrobinia: obsolete
#define LP_MP_CONSTRUCT_SOLVER(FMC) \
using namespace LP_MP; \
int main(int argc, char * argv[]) \
{ \
   ProblemDecomposition<FMC> pd(argc,argv); \
   pd.Construct(); \
   return pd.Run(); \
}

// Macro for generating main function reading in file and then solving problem, do zrobinia: Supply standard visitor
#define LP_MP_CONSTRUCT_SOLVER_WITH_INPUT(FMC,PARSE_PROBLEM_FUNCTION) \
using namespace LP_MP; \
int main(int argc, char * argv[]) \
{ \
   TCLAP::CmdLine cmd(std::string("Command line options for ") + FMC::name, ' ', "0.0.1"); \
   ProblemDecomposition<FMC> pd(cmd); \
   StandardVisitor visitor(cmd); \
   cmd.parse(argc,argv); \
   pd.ParseProblem(PARSE_PROBLEM_FUNCTION); \
   return pd.Run(visitor); \
}

// Macro which takes visitor as third argument
// do zrobienia: get version number automatically from CMake 
// do zrobienia: version number for individual solvers?
#define LP_MP_CONSTRUCT_SOLVER_WITH_INPUT_AND_VISITOR(FMC,PARSE_PROBLEM_FUNCTION,VISITOR) \
using namespace LP_MP; \
int main(int argc, char* argv[]) \
{ \
   TCLAP::CmdLine cmd(std::string("Command line options for ") + FMC::name, ' ', "0.0.1"); \
   ProblemDecomposition<FMC> pd(cmd); \
   VISITOR visitor(cmd,pd); \
   cmd.parse(argc,argv); \
   pd.ParseProblem(PARSE_PROBLEM_FUNCTION); \
   return pd.Run(visitor); \
}

} // end namespace LP_MP

#endif // LP_MP_PROBLEM_DECOMPOSITION_HXX
*/
