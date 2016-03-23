#ifndef LP_MP_SOLVER_HXX
#define LP_MP_SOLVER_HXX

namespace LP_MP {

// abstract class solver containing the FactorList, MessageList, ProblemConstructorList, Visitor, InputFunction
// coordinates between these three entities
template<typename FMC, typename INPUT_FUNCTION, typename VISITOR>
class Solver {
   std::array<INDEX,meta::size<ProblemDecompositionList>::value> factorIndex_; // gives for each problem constructor the offset, from which it started to add factors to the LP class

public:
   using InputFunction = INPUT_FUNTION; // do zrobienia: this is not so easy to do as a predefined template
   using VisitorType = VISITOR;

   Solver(TCLAP::CmdLine& cmd)
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



private:
   LP* lp_;
   tuple_from_list<ProblemDecompositionList> problem_constructor_;

   // command line arguments
   TCLAP::ValueArg<std::string> inputFileArg_;
   TCLAP::ValueArg<std::string> outputFileArg_;
   std::string inputFile_;
   std::string outputFile_;

}


} // end namespace LP_MP

#endif // LP_MP_SOLVER_HXX

