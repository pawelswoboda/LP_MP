#ifndef LP_MP_SOLVER_HXX
#define LP_MP_SOLVER_HXX

#include "LP_MP.h"
#include "meta/meta.hpp"
#include "function_existence.hxx"
#include "template_utilities.hxx"
#include "tclap/CmdLine.h"
#include "lp_interface/lp_interface.h"
#include "boost/hana.hpp"

namespace LP_MP {

namespace hana = boost::hana; // put this somewhere more central
// class containing the LP, problem constructor list, input function
// takes care of logging
// binds together problem constructors and solver and organizes input/output
// base class for solvers with primal rounding, e.g. LP-based rounding heuristics, message passing rounding and rounding provided by problem constructors.
template<typename FMC>
class Solver {
   // initialize a tuple uniformly
   //template <class T, class LIST>
   //   tuple_from_list<LIST> tupleMaker(LIST l, T& t) { return tuple_from_list<LIST>(t); }
   template <class T, class... ARGS>
      std::tuple<ARGS...> tupleMaker(meta::list<ARGS...>, T& t) { return std::make_tuple(ARGS(t)...); }

public:
   using ProblemDecompositionList = typename FMC::ProblemDecompositionList;
   using FactorMessageConnection = FMC;

   Solver()
      : cmd_(std::string("Command line options for ") + FMC::name, ' ', "0.0.1"),
      lp_(LP()),
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

   virtual ~Solver() {}

   // needed, as more arguments could be passed to cmd_, and then we need to parse again
   void Init(std::vector<std::string> arg)
   { 
      cmd_.parse(arg);
      Init_();
   }

   void Init(int argc, char** argv)
   {
      cmd_.parse(argc,argv);
      Init_();
   }
   void Init_()
   {
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
      else return HasTighten<meta::at_c<ProblemDecompositionList,n>, INDEX, INDEX>();
      //static_assert(PROBLEM_CONSTRUCTOR_NO<ProblemDecompositionList::size(),"");
   }
   template<INDEX PROBLEM_CONSTRUCTOR_NO>
   typename std::enable_if<PROBLEM_CONSTRUCTOR_NO >= ProblemDecompositionList::size(),INDEX>::type
   Tighten(const INDEX maxConstraints) { return 0; }
   template<INDEX PROBLEM_CONSTRUCTOR_NO>
   typename std::enable_if<PROBLEM_CONSTRUCTOR_NO < ProblemDecompositionList::size() && !CanTighten<PROBLEM_CONSTRUCTOR_NO>(),INDEX>::type
   Tighten(const INDEX maxConstraints)
   {
      return Tighten<PROBLEM_CONSTRUCTOR_NO+1>(maxConstraints);
   }
   template<INDEX PROBLEM_CONSTRUCTOR_NO>
   typename std::enable_if<PROBLEM_CONSTRUCTOR_NO < ProblemDecompositionList::size() && CanTighten<PROBLEM_CONSTRUCTOR_NO>(),INDEX>::type
   Tighten(const INDEX maxConstraints) 
   {
      std::cout << "Tighten for pc no " << PROBLEM_CONSTRUCTOR_NO << "\n";
      const INDEX noCuttingPlaneAdded = std::get<PROBLEM_CONSTRUCTOR_NO>(problemConstructor_).Tighten(maxConstraints);
      return noCuttingPlaneAdded + Tighten<PROBLEM_CONSTRUCTOR_NO+1>(maxConstraints);
   }
   // maxConstraints gives maximum number of constraints to add
   INDEX Tighten(const INDEX maxConstraints) 
   {
      INDEX noConstraintsAdded = Tighten<0>(maxConstraints);
      if(noConstraintsAdded > 0) { // tell lp to rebuild omegas etc
         lp_.Init();
      }
      return noConstraintsAdded;
   }
   
   template<INDEX PROBLEM_CONSTRUCTOR_NO>
   meta::at_c<ProblemDecompositionList, PROBLEM_CONSTRUCTOR_NO>& GetProblemConstructor() 
   {
      return std::get<PROBLEM_CONSTRUCTOR_NO>(problemConstructor_);
   }

   /* for hana
   template<INDEX PROBLEM_CONSTRUCTOR_NO>
   auto& GetProblemConstructor()
   {   
      std::get<PROBLEM_CONSTRUCTOR_NO>(problemConstructor_); 
   }
   */

   LP& GetLP() { return lp_; }
   
   // what to do before improving lower bound, e.g. setting reparametrization mode
   virtual void PreIterate(LpControl c) 
   {
      lp_.ComputeWeights(c.repam);
   } 

   // what to do for improving lower bound, typically ComputePass or ComputePassAndPrimal
   virtual void Iterate(LpControl c) {
      lp_.ComputePass();
   } 

   // what to do after one iteration of message passing, e.g. primal computation and/or tightening
   virtual void PostIterate(LpControl c) 
   {
      if(c.computeLowerBound) {
         lowerBound_ = lp_.LowerBound();
      }
      if(c.tighten) {
         Tighten(c.tightenConstraints);
      }
   } 

   void RegisterPrimal(PrimalSolutionStorage& p)
   {
      //std::cout << "size of primal = " << p.size() << "\n";
      const REAL cost = lp_.EvaluatePrimal(p.begin());
      if(cost <= bestPrimalCost_) {
         // check constraints
         if(CheckPrimalConsistency(p.begin())) {
            //std::cout << "primal cost = " << cost << ",, solution improved, primal solution feasible. in register primal\n";
            bestPrimalCost_ = cost;
            std::swap(bestPrimal_, p); // note: the best primal need not be admissible for the current lp, i.e. after tightening, the lp has changed, while best primal possibly has steyed the same.
         } else {
            //std::cout << "primal cost = " << cost << ", solution improved, primal solution infeasible. in register primal\n";
         }
      } else {
         //std::cout << "primal cost = " << cost << ", solution not improved. in register primal\n";
      }
   }
protected:
   TCLAP::CmdLine cmd_;
   LP lp_;
   //std::shared_ptr<spdlog::logger> logger_;
   //typename FMC::ProblemDecompositionListHana problemConstructor_;

   tuple_from_list<ProblemDecompositionList> problemConstructor_;

   // command line arguments
   TCLAP::ValueArg<std::string> inputFileArg_;
   TCLAP::ValueArg<std::string> outputFileArg_;
   std::string inputFile_;
   std::string outputFile_;

   REAL lowerBound_;
   // while Solver does not know how to compute primal, derived solvers do know. After computing a primal, they are expected to register their primals with the base solver
   REAL bestPrimalCost_ = std::numeric_limits<REAL>::infinity();
   PrimalSolutionStorage bestPrimal_; // these vectors are stored in the order of forwardOrdering_

};

// local rounding interleaved with message passing 
template<typename FMC>
class MpRoundingSolver : public Solver<FMC>
{
public:
   void Iterate(LpControl c)
   {
      if(c.computePrimal) {
         Solver<FMC>::lp_.ComputePassAndPrimal(forwardPrimal_, backwardPrimal_);
         RegisterPrimal(forwardPrimal_);
         RegisterPrimal(backwardPrimal_);
      } else {
         Solver<FMC>::Iterate(c);
      }
   }

private:
   PrimalSolutionStorage forwardPrimal_, backwardPrimal_;
};

// rounding based on primal heuristics provided by problem constructor
template<typename FMC>
class ProblemConstructorRoundingSolver : public Solver<FMC>
{
public:
   using Solver<FMC>::Solver;
   LP_MP_FUNCTION_EXISTENCE_CLASS(HasComputePrimal,ComputePrimal);
   template<INDEX PROBLEM_CONSTRUCTOR_NO>
   constexpr static bool
   CanComputePrimal()
   {
      // do zrobienia: this is not nice. CanComputePrimal should only be called with valid PROBLEM_CONSTRUCTOR_NO
      constexpr INDEX n = PROBLEM_CONSTRUCTOR_NO >= Solver<FMC>::ProblemDecompositionList::size() ? 0 : PROBLEM_CONSTRUCTOR_NO;
      if(n < PROBLEM_CONSTRUCTOR_NO) return false;
      else return HasComputePrimal<meta::at_c<typename Solver<FMC>::ProblemDecompositionList,n>, void, PrimalSolutionStorage::Element>();
      //static_assert(PROBLEM_CONSTRUCTOR_NO<ProblemDecompositionList::size(),"");
   }
   template<INDEX PROBLEM_CONSTRUCTOR_NO>
   typename std::enable_if<PROBLEM_CONSTRUCTOR_NO >= Solver<FMC>::ProblemDecompositionList::size()>::type
   ComputePrimal(PrimalSolutionStorage::Element primal) { return; }
   template<INDEX PROBLEM_CONSTRUCTOR_NO>
   typename std::enable_if<PROBLEM_CONSTRUCTOR_NO < Solver<FMC>::ProblemDecompositionList::size() && !CanComputePrimal<PROBLEM_CONSTRUCTOR_NO>()>::type
   ComputePrimal(PrimalSolutionStorage::Element primal)
   {
      return ComputePrimal<PROBLEM_CONSTRUCTOR_NO+1>(primal);
   }
   template<INDEX PROBLEM_CONSTRUCTOR_NO>
   typename std::enable_if<PROBLEM_CONSTRUCTOR_NO < Solver<FMC>::ProblemDecompositionList::size() && CanComputePrimal<PROBLEM_CONSTRUCTOR_NO>()>::type
   ComputePrimal(PrimalSolutionStorage::Element primal)
   {
      std::cout << "ComputePrimal for pc no " << PROBLEM_CONSTRUCTOR_NO << "\n";
      std::get<PROBLEM_CONSTRUCTOR_NO>(this->problemConstructor_).ComputePrimal(primal);
      return ComputePrimal<PROBLEM_CONSTRUCTOR_NO+1>(primal);
   }
   void ComputePrimal(PrimalSolutionStorage::Element primal)
   {
      ComputePrimal<0>(primal);
   }

   virtual void PostIterate(LpControl c)
   {
      if(c.computePrimal) {
         this->lp_.InitializePrimalVector(primal_); // do zrobienia: this is not nice: reallocation might occur
         ComputePrimal(primal_.begin());
         this->RegisterPrimal(primal_);
      }
      Solver<FMC>::PostIterate(c);
   }
   
private:
   PrimalSolutionStorage primal_;
};



// solver holding visitor. We do not want to have this in base class, as this would entail passing visitor information to constructors etc.
template<typename SOLVER, typename VISITOR>
class VisitorSolver : public SOLVER {
public:
   VisitorSolver(int argc, char** argv) :
      SOLVER(),
      visitor_(SOLVER::cmd_)
   {
      this->Init(argc,argv);
   }
   VisitorSolver(const std::vector<std::string>& arg) :
      SOLVER(),
      visitor_(SOLVER::cmd_)
   {
      this->Init(arg);
   }
   int Solve()
   {
      this->lp_.Init();
      LpControl c = visitor_.begin(this->lp_);
      while(!c.end) {
         this->PreIterate(c);
         this->Iterate(c);
         this->PostIterate(c);
         c = visitor_.visit(c, this->lowerBound_, this->bestPrimalCost_);
      }
      return c.error;
   }

private:
  VISITOR visitor_;
};

template<typename SOLVER, typename VISITOR, typename LpSolver>
class VisitorLpSolver : public SOLVER {
public:
   VisitorLpSolver(int argc, char** argv) :
      SOLVER(),
      visitor_(SOLVER::cmd_),
      LP_("","lp","use reparametrization for generating an lp",SOLVER::cmd_),
      LPOnly_("","onlyLp","using lp solver without reparametrization",SOLVER::cmd_,true),
      RELAX_("","relax","solve the mip relaxation",SOLVER::cmd_),
      timelimit_("","LpTimelimit","timelimit for the lp solver",false,3600.0,"positive real number",SOLVER::cmd_),
      threads_("","LpThreads","number of threads used by the lp solver",false,1,"integer",SOLVER::cmd_)
   {
      this->Init(argc,argv);
   }
   VisitorLpSolver(const std::vector<std::string>& arg) :
      SOLVER(),
      visitor_(SOLVER::cmd_),
      LP_("","lp","use reparametrization for generating an lp",SOLVER::cmd_),
      LPOnly_("","onlyLp","using lp solver without reparametrization",SOLVER::cmd_,true),
      RELAX_("","relax","solve the mip relaxation",SOLVER::cmd_),
      timelimit_("","LpTimelimit","timelimit for the lp solver",false,3600.0,"positive real number",SOLVER::cmd_),
      threads_("","LpThreads","number of threads used by the lp solver",false,1,"integer",SOLVER::cmd_)
   {
      this->Init(arg);
   }
   
   ~VisitorLpSolver(){
     delete lpSolver_;
   }
   
    template<class T,class E>
    class FactorMessageIterator {
    public:
      FactorMessageIterator(T w,INDEX i)
        : wrapper_(w),idx(i){ }
      FactorMessageIterator& operator++() {++idx;return *this;}
      FactorMessageIterator operator++(int) {FactorMessageIterator tmp(*this); operator++(); return tmp;}
      bool operator==(const FactorMessageIterator& rhs) {return idx==rhs.idx;}
      bool operator!=(const FactorMessageIterator& rhs) {return idx!=rhs.idx;}
      E* operator*() {return wrapper_(idx);}
      E* operator->() {return wrapper_(idx);}
      INDEX idx;
    private:
      T wrapper_;
    };
   
   int Solve() 
   {
      this->lp_.Init(); 
      LpControl c = visitor_.begin(this->lp_);
      if(LPOnly_.getValue()){        
        while(!c.end) {
           this->PreIterate(c);
           this->Iterate(c);
           this->PostIterate(c);
           c = visitor_.visit(c, this->lowerBound_, this->bestPrimalCost_);
        }
      }
      
      if(LP_.getValue() || !LPOnly_.getValue()){
        auto FactorWrapper = [&](INDEX i){ return SOLVER::lp_.GetFactor(i);};
        FactorMessageIterator<decltype(FactorWrapper),FactorTypeAdapter> FactorItBegin(FactorWrapper,0);
        FactorMessageIterator<decltype(FactorWrapper),FactorTypeAdapter> FactorItEnd(FactorWrapper,SOLVER::lp_.GetNumberOfFactors());

        auto MessageWrapper = [&](INDEX i){ return SOLVER::lp_.GetMessage(i);};
        FactorMessageIterator<decltype(MessageWrapper),MessageTypeAdapter> MessageItBegin(MessageWrapper,0);
        FactorMessageIterator<decltype(MessageWrapper),MessageTypeAdapter> MessageItEnd(MessageWrapper,SOLVER::lp_.GetNumberOfMessages());
        
        lpSolver_ = new LpSolver(FactorItBegin,FactorItEnd,MessageItBegin,MessageItEnd,!RELAX_.getValue());
        lpSolver_->SetTimeLimit(timelimit_.getValue());
        lpSolver_->SetNumberOfThreads(threads_.getValue());
        
        lpSolver_->solve();
      }
      
      return c.error;
   }

private:
  VISITOR visitor_;
  LpInterfaceAdapter* lpSolver_;
  
  TCLAP::ValueArg<REAL> timelimit_;
  TCLAP::ValueArg<INDEX> threads_;
  TCLAP::SwitchArg LP_;
  TCLAP::SwitchArg LPOnly_;
  TCLAP::SwitchArg RELAX_;
  
};

  // solver for rounding with standard (I)LP-solver
  template<typename FMC,class LpSolver>
  class LpRoundingSolver : public Solver<FMC> {

  public:
    LpRoundingSolver(int argc, char** argv)
      : Solver<FMC>(),
      LpOutputFile_("","lpFile","write LP model",false,"","file name",Solver<FMC>::cmd_),
      MIP_("","relax","solve the mip relaxation",Solver<FMC>::cmd_,true){
      Solver<FMC>::Init(argc,argv);
    }

    LpRoundingSolver(const std::vector<std::string>& arg)
      : Solver<FMC>(),
        LpOutputFile_("","lpFile","write LP model",false,"","file name",Solver<FMC>::cmd_),
        MIP_("","relax","solve the mip relaxation",Solver<FMC>::cmd_,true){
      Solver<FMC>::Init(arg);
    }
    
    ~LpRoundingSolver(){
      delete lpSolver_;
    }
    
    template<class T,class E>
    class FactorMessageIterator {
    public:
      FactorMessageIterator(T w,INDEX i)
        : wrapper_(w),idx(i){ }
      FactorMessageIterator& operator++() {++idx;return *this;}
      FactorMessageIterator operator++(int) {FactorMessageIterator tmp(*this); operator++(); return tmp;}
      bool operator==(const FactorMessageIterator& rhs) {return idx==rhs.idx;}
      bool operator!=(const FactorMessageIterator& rhs) {return idx!=rhs.idx;}
      E* operator*() {return wrapper_(idx);}
      E* operator->() {return wrapper_(idx);}
      INDEX idx;
    private:
      T wrapper_;
    };

    LpInterfaceAdapter* GetLpSolver(){ return lpSolver_; };
    
    int Solve(){
      
      Solver<FMC>::lp_.Init();
      int status = 0;
      
      auto FactorWrapper = [&](INDEX i){ return Solver<FMC>::lp_.GetFactor(i);};
      FactorMessageIterator<decltype(FactorWrapper),FactorTypeAdapter> FactorItBegin(FactorWrapper,0);
      FactorMessageIterator<decltype(FactorWrapper),FactorTypeAdapter> FactorItEnd(FactorWrapper,Solver<FMC>::lp_.GetNumberOfFactors());

      auto MessageWrapper = [&](INDEX i){ return Solver<FMC>::lp_.GetMessage(i);};
      FactorMessageIterator<decltype(MessageWrapper),MessageTypeAdapter> MessageItBegin(MessageWrapper,0);
      FactorMessageIterator<decltype(MessageWrapper),MessageTypeAdapter> MessageItEnd(MessageWrapper,Solver<FMC>::lp_.GetNumberOfMessages());
      
      lpSolver_ = new LpSolver(FactorItBegin,FactorItEnd,MessageItBegin,MessageItEnd,MIP_.getValue());
      
      if( LpOutputFile_.isSet() ){
        lpSolver_->WriteLpModel(LpOutputFile_.getValue());
      } else {
        status = lpSolver_->solve();
      }
      
      return status;
    }
  
  private:
    LpInterfaceAdapter* lpSolver_;
    PrimalSolutionStorage lpPrimal_;

    TCLAP::ValueArg<std::string> LpOutputFile_;
    TCLAP::SwitchArg MIP_;
    //TCLAP::SwitchArg intConstr ("","integer","solving a mip",cmd);
  
  };

  template<class FMC,class LPSolver>
  class SolverLP : public Solver<FMC> {
  public:
    SolverLP(int argc, char** argv)
      : Solver<FMC>(),LpOutputFile_("","lpFile","write LP model",false,"","file name",Solver<FMC>::cmd_){
      Solver<FMC>::Init(argc,argv);
    }

    ~SolverLP(){
      delete lpSolver_;
    }
    
    void ComputePrimal(PrimalSolutionStorage::Element primal)
    {
      lpSolver_->solve(lpPrimal_.begin());
    }

    template<class T,class E>
    class FactorMessageIterator {
    public:
      FactorMessageIterator(T w,INDEX i)
        : wrapper_(w),idx(i){ }
      FactorMessageIterator& operator++() {++idx;return *this;}
      FactorMessageIterator operator++(int) {FactorMessageIterator tmp(*this); operator++(); return tmp;}
      bool operator==(const FactorMessageIterator& rhs) {return idx==rhs.idx;}
      bool operator!=(const FactorMessageIterator& rhs) {return idx!=rhs.idx;}
      E* operator*() {return wrapper_(idx);}
      E* operator->() {return wrapper_(idx);}
      INDEX idx;
    private:
      T wrapper_;
    };
    
    int Solve(){
            
      auto FactorWrapper = [&](INDEX i){ return Solver<FMC>::lp_.GetFactor(i);};
      FactorMessageIterator<decltype(FactorWrapper),FactorTypeAdapter> FactorItBegin(FactorWrapper,0);
      FactorMessageIterator<decltype(FactorWrapper),FactorTypeAdapter> FactorItEnd(FactorWrapper,Solver<FMC>::lp_.GetNumberOfFactors());

      auto MessageWrapper = [&](INDEX i){ return Solver<FMC>::lp_.GetMessage(i);};
      FactorMessageIterator<decltype(MessageWrapper),MessageTypeAdapter> MessageItBegin(MessageWrapper,0);
      FactorMessageIterator<decltype(MessageWrapper),MessageTypeAdapter> MessageItEnd(MessageWrapper,Solver<FMC>::lp_.GetNumberOfMessages());
      
      lpSolver_ = new LPSolver(FactorItBegin,FactorItEnd,MessageItBegin,MessageItEnd);
      if( LpOutputFile_.isSet() ){
        lpSolver_->WriteLpModel(LpOutputFile_.getValue());
      }
      
      return 0;
    }
    
  private:
    LpInterfaceAdapter* lpSolver_;
    PrimalSolutionStorage lpPrimal_;

    TCLAP::ValueArg<std::string> LpOutputFile_;
    
  };


// Macro for generating main function 
// do zrobienia: get version number automatically from CMake 
// do zrobienia: version number for individual solvers?
#define LP_MP_CONSTRUCT_SOLVER_WITH_INPUT_AND_VISITOR(FMC,PARSE_PROBLEM_FUNCTION,VISITOR) \
using namespace LP_MP; \
int main(int argc, char* argv[]) \
{ \
   VisitorSolver<Solver<FMC>,VISITOR> solver(argc,argv); \
   solver.ReadProblem(PARSE_PROBLEM_FUNCTION); \
   return solver.Solve(); \
}

// Macro for generating main function with specific primal rounding solver
#define LP_MP_CONSTRUCT_SOLVER_WITH_INPUT_VISITOR_AND_SOLVER(FMC,PARSE_PROBLEM_FUNCTION,VISITOR,SOLVER) \
using namespace LP_MP; \
int main(int argc, char* argv[]) \
{ \
   VisitorSolver<SOLVER,VISITOR> solver(argc,argv); \
   solver.ReadProblem(PARSE_PROBLEM_FUNCTION); \
   return solver.Solve(); \
}

// Macro for generating main function with specific solver with additional LP option
#define LP_MP_CONSTRUCT_SOLVER_WITH_INPUT_AND_LPVISITOR(FMC,PARSE_PROBLEM_FUNCTION,VISITOR,LPSOLVER) \
using namespace LP_MP; \
int main(int argc, char* argv[]) \
{ \
   VisitorLpSolver<Solver<FMC>,VISITOR,LPSOLVER> solver(argc,argv); \
   solver.ReadProblem(PARSE_PROBLEM_FUNCTION); \
   return solver.Solve(); \
}


} // end namespace LP_MP

#endif // LP_MP_SOLVER_HXX

