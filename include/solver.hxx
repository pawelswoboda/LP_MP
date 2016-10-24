#ifndef LP_MP_SOLVER_HXX
#define LP_MP_SOLVER_HXX

#include <thread>
#include <mutex>
#include <condition_variable>
#include <fstream>

#include "LP_MP.h"
#include "meta/meta.hpp"
#include "function_existence.hxx"
#include "template_utilities.hxx"
#include "tclap/CmdLine.h"
#include "lp_interface/lp_interface.h"

namespace LP_MP {

// class containing the LP, problem constructor list, input function
// takes care of logging
// binds together problem constructors and solver and organizes input/output
// base class for solvers with primal rounding, e.g. LP-based rounding heuristics, message passing rounding and rounding provided by problem constructors.

// do zrobienia: currently we synchronize in all solvers. However, strictly speaking, this is only needed by  LpSolver and MpRoundingSolver. Should only they synchronize?

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

   ~Solver() {}

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
      std::unique_lock<std::mutex> guard(LpChangeMutex);
      const bool success = inputFct(inputFile_,*this,args...);

      assert(success);
      if(!success) throw std::runtime_error("could not parse problem file");

      //spdlog::get("logger")->info("loading file " + inputFile_ + " succeeded");
      return success;
   }

   LP_MP_FUNCTION_EXISTENCE_CLASS(HasWritePrimal,WritePrimal);
   template<INDEX PROBLEM_CONSTRUCTOR_NO>
   constexpr static bool
   CanWritePrimal()
   {
      // do zrobienia: this is not nice. CanTighten should only be called with valid PROBLEM_CONSTRUCTOR_NO
      constexpr INDEX n = PROBLEM_CONSTRUCTOR_NO >= ProblemDecompositionList::size() ? 0 : PROBLEM_CONSTRUCTOR_NO;
      if(n < PROBLEM_CONSTRUCTOR_NO) return false;
      else return HasWritePrimal<meta::at_c<ProblemDecompositionList,n>, void, std::ofstream, PrimalSolutionStorage>();
      //static_assert(PROBLEM_CONSTRUCTOR_NO<ProblemDecompositionList::size(),"");
   }
   template<INDEX PROBLEM_CONSTRUCTOR_NO>
   typename std::enable_if<PROBLEM_CONSTRUCTOR_NO >= ProblemDecompositionList::size()>::type
   WritePrimal(std::ofstream& s) {}
   template<INDEX PROBLEM_CONSTRUCTOR_NO>
   typename std::enable_if<PROBLEM_CONSTRUCTOR_NO < ProblemDecompositionList::size() && !CanWritePrimal<PROBLEM_CONSTRUCTOR_NO>()>::type
   WritePrimal(std::ofstream& s)
   {
      return WritePrimal<PROBLEM_CONSTRUCTOR_NO+1>(s);
   }
   template<INDEX PROBLEM_CONSTRUCTOR_NO>
   typename std::enable_if<PROBLEM_CONSTRUCTOR_NO < ProblemDecompositionList::size() && CanWritePrimal<PROBLEM_CONSTRUCTOR_NO>()>::type
   WritePrimal(std::ofstream& s) 
   {
      std::cout << "WritePrimal for pc no " << PROBLEM_CONSTRUCTOR_NO << "\n";
      std::get<PROBLEM_CONSTRUCTOR_NO>(problemConstructor_).WritePrimal(s,bestPrimal_);
   }

   void WritePrimal()
   {
      if(outputFileArg_.isSet()) {
         std::ofstream output_file;
         output_file.open(outputFile_);
         if(!output_file.is_open()) {
            throw std::runtime_error("could not open file " + outputFile_);
         }
         WritePrimal<0>(output_file);
      }
   }


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

   LP& GetLP() { return lp_; }
   
   // called before first iterations
   void Begin() {}

   // what to do before improving lower bound, e.g. setting reparametrization mode
   void PreIterate(LpControl c) 
   {
      std::unique_lock<std::mutex> guard(LpChangeMutex);
      lp_.ComputeWeights(c.repam);
   } 

   // what to do for improving lower bound, typically ComputePass or ComputePassAndPrimal
   void Iterate(LpControl c) {
      std::unique_lock<std::mutex> guard(LpChangeMutex);
      lp_.ComputePass();
   } 

   // what to do after one iteration of message passing, e.g. primal computation and/or tightening
   void PostIterate(LpControl c) 
   {
      if(c.computeLowerBound) {
         std::unique_lock<std::mutex> guard(LpChangeMutex);
         lowerBound_ = lp_.LowerBound();
      }
      if(c.tighten) {
         std::unique_lock<std::mutex> guard(LpChangeMutex);
         Tighten(c.tightenConstraints);
      }
   } 

   // called after last iteration
   void End() {}

   void RegisterPrimal(PrimalSolutionStorage& p)
   {
    std::unique_lock<std::mutex> guard(LpChangeMutex);
    //std::cout << "size of primal = " << p.size() << "\n";
    const REAL cost = lp_.EvaluatePrimal(p.begin());
    printf("cost: %.2f \n",cost);
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

   REAL lower_bound() const { return lowerBound_; }
protected:
   TCLAP::CmdLine cmd_;
   LP lp_;

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
   std::mutex LpChangeMutex;
};

// local rounding interleaved with message passing 
template<typename FMC>
class MpRoundingSolver : public Solver<FMC>
{
public:
   void Iterate(LpControl c)
   {
      if(c.computePrimal) {
         {
                 std::unique_lock<std::mutex> guard(this->LpChangeMutex);
                 Solver<FMC>::lp_.ComputePassAndPrimal(forwardPrimal_, backwardPrimal_);
         }
         this->RegisterPrimal(forwardPrimal_);
         this->RegisterPrimal(backwardPrimal_);
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

   void PostIterate(LpControl c)
   {
      if(c.computePrimal) {
         std::unique_lock<std::mutex> guard(this->LpChangeMutex);
         this->lp_.InitializePrimalVector(primal_); // do zrobienia: this is not nice: reallocation might occur
         // do zrobienia: possibly run this in own thread similar to lp solver
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
    this->Begin();
    std::unique_lock<std::mutex> VisitorGuard(VisitorMutex_);
    LpControl c = visitor_.begin(this->lp_);
    VisitorGuard.unlock();
      
    while(!c.end && !c.error) {
      this->PreIterate(c);
      this->Iterate(c);
      this->PostIterate(c);

      std::unique_lock<std::mutex> VisitorGuard(VisitorMutex_);
      c = visitor_.visit(c, this->lowerBound_, this->bestPrimalCost_);
    }
    if(!c.error) {
      this->End();
      // possibly primal has been computed in end. Call visitor again
      std::unique_lock<std::mutex> VisitorGuard(VisitorMutex_);
      visitor_.end(this->lowerBound_, this->bestPrimalCost_);
      this->WritePrimal();
    }
    return c.error;
  }

private:
  VISITOR visitor_;
  static std::mutex VisitorMutex_;
};
  template<typename SOLVER, typename VISITOR>
  std::mutex VisitorSolver<SOLVER,VISITOR>::VisitorMutex_;
  
template<typename SOLVER, typename LP_INTERFACE>
class LpSolver : public SOLVER {
public:
   LpSolver() :
      SOLVER(),
      LPOnly_("","onlyLp","using lowerbound from lp solver instead of message passing",SOLVER::cmd_),
      RELAX_("","relax","solve the mip relaxation",SOLVER::cmd_),
      timelimit_("","LpTimelimit","timelimit for the lp solver",false,3600.0,"positive real number",SOLVER::cmd_),
      roundBound_("","LpRoundValue","A small value removes many variables",false,std::numeric_limits<REAL>::infinity(),"positive real number",SOLVER::cmd_),
      threads_("","LpSolverThreads","number of threads to call Lp solver routine",false,1,"integer",SOLVER::cmd_),
      LpThreadsArg_("","LpThreads","number of threads used by the lp solver",false,1,"integer",SOLVER::cmd_),
      LpInterval_("","LpInterval","each n steps the lp solver will be executed if possible",false,1,"integer",SOLVER::cmd_)
    {}
      
    ~LpSolver(){
      runLp_ = false;
      std::unique_lock<std::mutex> WakeUpGuard(WakeLpSolverMutex_);
      WakeLpSolverCond.notify_all();
      WakeUpGuard.unlock();

      for(INDEX i=0;i<LpThreads_.size();i++){
        if(LpThreads_[i].joinable()){
          LpThreads_[i].join();
        }
      }
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
   
    void RunLpSolver(int displayLevel=0)
    {
            std::unique_lock<std::mutex> guard(this->LpChangeMutex);
            // Jan: lp can change its number of factors, e.g. after tightening. We must synchronize here as well!
            auto FactorWrapper = [&](INDEX i){ return SOLVER::lp_.GetFactor(i);};
            FactorMessageIterator<decltype(FactorWrapper),FactorTypeAdapter> FactorItBegin(FactorWrapper,0);
            FactorMessageIterator<decltype(FactorWrapper),FactorTypeAdapter> FactorItEnd(FactorWrapper,SOLVER::lp_.GetNumberOfFactors());

            auto MessageWrapper = [&](INDEX i){ return SOLVER::lp_.GetMessage(i);};
            FactorMessageIterator<decltype(MessageWrapper),MessageTypeAdapter> MessageItBegin(MessageWrapper,0);
            FactorMessageIterator<decltype(MessageWrapper),MessageTypeAdapter> MessageItEnd(MessageWrapper,SOLVER::lp_.GetNumberOfMessages());

      LP_INTERFACE solver(FactorItBegin,FactorItEnd,MessageItBegin,MessageItEnd,VariableThreshold_,!RELAX_.getValue());

      //solver.ReduceLp(FactorItBegin,FactorItEnd,MessageItBegin,MessageItEnd,VariableThreshold_);
      guard.unlock();

            solver.SetTimeLimit(timelimit_.getValue());
            solver.SetNumberOfThreads(LpThreadsArg_.getValue());
            solver.SetDisplayLevel(displayLevel);

            auto status = solver.solve();
            std::cout << "solved lp instance\n";
            //  TODO: we need a better way to decide, when the number of variables get increased/decreased
            //  e.g. What to do if the solver hit the timelimit? 
            if( status == 0 || status == 3 ){     
                    // Jan: this is only admissible when we have an integral solution, not in RELAXed mode!       
                    PrimalSolutionStorage x(FactorItBegin,FactorItEnd);
                    for(INDEX i=0;i<x.size();i++){
                      x[i] = std::round(solver.GetVariableValue(i));
                    }
                    this->RegisterPrimal(x);
                    //std::lock_guard<std::mutex> WriteGuard(BuildLpMutex);
                    if(LPOnly_.getValue()){
                      this->lowerBound_ = solver.GetBestBound();
                    }

                    if( status == 0){
                            std::cout << "Optimal solution found by Lp Solver with value: " << solver.GetObjectiveValue() << std::endl;
                    } else {
                            std::cout << "Suboptimal solution found by Lp Solver with value: " << solver.GetObjectiveValue() << std::endl;
                    }

                    std::lock_guard<std::mutex> lck(UpdateUbLpMutex);
                    // if solution found, decrease number of variables
                    std::cout << "feas decrease" << std::endl;
                    VariableThreshold_ *= 0.75;
                    ub_ = std::min(ub_,solver.GetObjectiveValue());
            }
            else if( status == 4){ // hit the timelimit
                    std::lock_guard<std::mutex> lck(UpdateUbLpMutex);
                    VariableThreshold_ *= 0.6;
                    std::cout << "time decrease" << std::endl;
            }
            else {
                    std::lock_guard<std::mutex> lck(UpdateUbLpMutex);
                    // if ilp infeasible, increase number of variables
                    std::cout << "infeas increase" << std::endl;
                    VariableThreshold_ *= 1.3;
            }
    }

    void IterateLpSolver()
    {
            while(runLp_) {

                    std::unique_lock<std::mutex> WakeLpGuard(WakeLpSolverMutex_);
                    WakeLpSolverCond.wait(WakeLpGuard);
                    WakeLpGuard.unlock();
                    if(!runLp_) { break; }
                    RunLpSolver();
            }
    }

    void Begin()
    {
            LpThreads_.resize(threads_.getValue());
            for(INDEX i=0;i<threads_.getValue();i++){
                    LpThreads_[i] = std::thread([this] () { this->IterateLpSolver(); });
            }

            VariableThreshold_ = roundBound_.getValue();
            SOLVER::Begin();
    }

    void Iterate(LpControl c)
    {
          // wait if some lp solver constructs lp. Reparametrization is not allowed to change then?
          SOLVER::Iterate(c);
    }
    void PostIterate(LpControl c)
    {
      // wait, as in PostIterate tightening can take place
      SOLVER::PostIterate(c);

      if( curIter_ % LpInterval_.getValue() == 0 ){
        std::unique_lock<std::mutex> WakeUpGuard(WakeLpSolverMutex_);
        WakeLpSolverCond.notify_one();
        std::this_thread::sleep_for(std::chrono::milliseconds(50)); // Without sleep, the next iteration will block the execution of the LpInterface
      }
      ++curIter_;
    }

    void End()
    {
      runLp_ = false;
      std::cout << "stop lp solvers\n";
      std::unique_lock<std::mutex> WakeUpGuard(WakeLpSolverMutex_);
      WakeLpSolverCond.notify_all();
      WakeUpGuard.unlock();

      std::cout << "construct final lp solver\n";
      std::thread th([this]() { this->RunLpSolver(1); });
      th.join();


      for(INDEX i=0;i<threads_.getValue();i++){
              if(LpThreads_[i].joinable()){
                      LpThreads_[i].join();
              }
      }

      
      
      SOLVER::End();
      std::cout << "\n";
      std::cout << "The best objective computed by ILP is " << ub_ << "\n";
    }

private:
  TCLAP::ValueArg<REAL> timelimit_;
  TCLAP::ValueArg<REAL> roundBound_;
  TCLAP::ValueArg<INDEX> threads_;
  TCLAP::ValueArg<INDEX> LpThreadsArg_;
  TCLAP::ValueArg<INDEX> LpInterval_;
  TCLAP::SwitchArg LPOnly_;
  TCLAP::SwitchArg RELAX_;
  
  std::mutex UpdateUbLpMutex;
  std::mutex WakeLpSolverMutex_;
  std::condition_variable WakeLpSolverCond;

  bool runLp_ = true;
  std::vector<std::thread> LpThreads_;
  REAL VariableThreshold_;
  REAL ub_ = std::numeric_limits<REAL>::infinity();
  INDEX curIter_ = 1;
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

#define LP_MP_CONSTRUCT_SOLVER_WITH_INPUT_AND_VISITOR_MP_ROUNDING(FMC,PARSE_PROBLEM_FUNCTION,VISITOR) \
using namespace LP_MP; \
int main(int argc, char* argv[]) \
{ \
   VisitorSolver<MpRoundingSolver<FMC>,VISITOR> solver(argc,argv); \
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
#define LP_MP_CONSTRUCT_SOLVER_WITH_INPUT_AND_VISITOR_LP(FMC,PARSE_PROBLEM_FUNCTION,VISITOR,LPSOLVER) \
using namespace LP_MP; \
int main(int argc, char* argv[]) \
{ \
   VisitorSolver<LpSolver<Solver<FMC>,LPSOLVER>,VISITOR> solver(argc,argv); \
   solver.ReadProblem(PARSE_PROBLEM_FUNCTION); \
   return solver.Solve(); \
}

} // end namespace LP_MP

#endif // LP_MP_SOLVER_HXX

