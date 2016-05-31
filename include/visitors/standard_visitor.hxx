#ifndef LP_MP_STANDARD_VISITOR_HXX
#define LP_MP_STANDARD_VISITOR_HXX

#include "LP_MP.h"
#include "tolerance.hxx"
#include "mem_use.c"
#include "tclap/CmdLine.h"
#include <chrono>
#include "spdlog/spdlog.h"

/*
 minimal visitor class:
class Visitor {
public:
   Visitor(TCLAP::CmdLine& cmd, PROBLEM_DECOMPOSITION& pd);
   LPVisitorReturnType begin(const LP* lp);
   template<LPVisitorReturnType LP_STATE>
   LPVisitorReturnType visit(LP* lp)

};
*/

namespace LP_MP {
   class PositiveRealConstraint : public TCLAP::Constraint<REAL>
   {
      public:
         std::string description() const { return "positive real constraint"; };
         std::string shortID() const { return "positive real number"; };
         bool check(const REAL& value) const { return value >= 0.0; };
   };
   class OpenUnitIntervalConstraint: public TCLAP::Constraint<REAL>
   {
      public:
         std::string description() const { return "0<x<1 real constraint"; };
         std::string shortID() const { return "positive real number smaller 1"; };
         bool check(const REAL& value) const { return value > 0.0 && value < 1.0; };
   };

   // standard visitor class for LP_MP solver, when no custom visitor is given
   // do zrobienia: add xor arguments primalBoundComputationInterval, dualBoundComputationInterval with boundComputationInterval
   template<class PROBLEM_DECOMPOSITION>
   class StandardVisitor {
      
      public:
      StandardVisitor(TCLAP::CmdLine& cmd, PROBLEM_DECOMPOSITION& pd)
         :
            maxIterArg_("","maxIter","maximum number of iterations of LP_MP, default = 1000",false,1000,"positive integer",cmd),
            maxMemoryArg_("","maxMemory","maximum amount of memory (MB) LP_MP is allowed to use",false,std::numeric_limits<INDEX>::max(),"positive integer",cmd),
            timeoutArg_("","timeout","time after which algorithm is stopped, in seconds, default = never, should this be type double?",false,std::numeric_limits<INDEX>::max(),"positive integer",cmd),
            // xor those //
            //boundComputationIntervalArg_("","boundComputationInterval","lower bound computation performed every x-th iteration, default = 5",false,5,"positive integer",cmd),
            primalComputationIntervalArg_("","primalComputationInterval","primal computation performed every x-th iteration, default = 5",false,5,"positive integer",cmd),
            lowerBoundComputationIntervalArg_("","lowerBoundComputationInterval","lower bound computation performed every x-th iteration, default = 1",false,1,"positive integer",cmd),
            ///////////////
            posConstraint_(),
            minDualImprovementArg_("","minDualImprovement","minimum dual improvement between iterations of LP_MP",false,0.0,&posConstraint_,cmd),
            standardReparametrizationArg_("","standardReparametrization","mode of reparametrization: {anisotropic,uniform}",false,"anisotropic","{anisotropic|uniform}",cmd),
            roundingReparametrizationArg_("","roundingReparametrization","mode of reparametrization for rounding primal solution: {anisotropic|uniform}",false,"uniform","{anisotropic|uniform}",cmd),
            protocolateConsoleArg_("","protocolateConsole","protocolate on console (stdout)",cmd,false),
            protocolateFileArg_("","protocolateFile","file into which to protocolate progress of algorithm",false,"","file name",cmd),
            pd_(pd)
      {}

      LPVisitorReturnType begin(const LP* lp) // called, after problem is constructed. 
      {
         try {
            maxIter_ = maxIterArg_.getValue();
            maxMemory_ = maxMemoryArg_.getValue();
            remainingIter_ = maxIter_;
            minDualImprovement_ = minDualImprovementArg_.getValue();
            timeout_ = timeoutArg_.getValue();
            //boundComputationInterval_ = boundComputationIntervalArg_.getValue();
            primalComputationInterval_ = primalComputationIntervalArg_.getValue();
            lowerBoundComputationInterval_ = lowerBoundComputationIntervalArg_.getValue();

            standardReparametrization_ = standardReparametrizationArg_.getValue();
            roundingReparametrization_ = roundingReparametrizationArg_.getValue();
            protocolateConsole_ = protocolateConsoleArg_.getValue();
            protocolateFile_ = protocolateFileArg_.getValue();
         } catch (TCLAP::ArgException &e) {
            std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl; 
            exit(1);
         }
         curDualBound_ = lp->LowerBound();

         try {
            std::vector<spdlog::sink_ptr> sinks;
            if(protocolateConsole_) {
               sinks.push_back(std::make_shared<spdlog::sinks::stdout_sink_st>());
            } 
            if(protocolateFile_ != "") {
               sinks.push_back(std::make_shared<spdlog::sinks::simple_file_sink_st>(protocolateFile_.c_str(),true));
            }
            auto logger = std::make_shared<spdlog::logger>("logger", std::begin(sinks), std::end(sinks));
            logger->set_pattern("%v");
            spdlog::register_logger(logger);
         } catch(const spdlog::spdlog_ex& ex) {
            std::cerr << "instantiating logger class failed: " << ex.what();
            throw std::runtime_error("could not instantiate logger");
         }


         if( ! (standardReparametrization_ == "anisotropic" || standardReparametrization_ == "uniform") ) {
            throw std::runtime_error("standard repararametrization mode must be {anisotropic|uniform}, is " + standardReparametrization_);
         }
         if( ! (roundingReparametrization_ == "anisotropic" || roundingReparametrization_ == "uniform") ) {
            throw std::runtime_error("rounding repararametrization mode must be {anisotropic|uniform}, is " + roundingReparametrization_);
         }

         spdlog::get("logger")->info() << "Initial number of factors = " << lp->GetNumberOfFactors();
         spdlog::get("logger")->info() << "Initial lower bound before optimizing = " << curDualBound_;
         beginTime_ = std::chrono::steady_clock::now();

         if(roundingReparametrization_ == "anisotropic") {
            return LPVisitorReturnType::ReparametrizeLowerBoundPrimalAnisotropic;
         } else if(roundingReparametrization_ == "uniform") {
            return LPVisitorReturnType::ReparametrizeLowerBoundPrimalUniform;
         } else {
            throw std::runtime_error("parametrization not recognized");
         }
      }


      // the default
      // do zrobienia: it makes little sense to have visit templatised. Just pass parameter as an argument
      template<LPVisitorReturnType LP_STATE>
      LPVisitorReturnType visit(LP* lp)
      {

         const INDEX timeElapsed = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - beginTime_).count();
         auto logger = spdlog::get("logger");

         if(LP_STATE == LPVisitorReturnType::ReparametrizeUniform || LP_STATE == LPVisitorReturnType::ReparametrizeAnisotropic) {
            // output nothing
         } else if(LP_STATE == LPVisitorReturnType::ReparametrizeLowerBoundUniform || LP_STATE == LPVisitorReturnType::ReparametrizeLowerBoundAnisotropic) {
            logger->info() << "iteration = " << curIter_ << ", lower bound = " << lp->BestLowerBound() << ", time elapsed = " << timeElapsed << " milliseconds";
         } else if(LP_STATE == LPVisitorReturnType::ReparametrizePrimalUniform || LP_STATE == LPVisitorReturnType::ReparametrizePrimalAnisotropic) {
            logger->info() << "iteration = " << curIter_ << ", upper bound = " << lp->BestPrimalBound() << ", time elapsed = " << timeElapsed << " milliseconds";
         } else if(LP_STATE == LPVisitorReturnType::ReparametrizeLowerBoundPrimalUniform || LP_STATE == LPVisitorReturnType::ReparametrizeLowerBoundPrimalAnisotropic) {
            logger->info() << "iteration = " << curIter_ << ", lower bound = " << lp->BestLowerBound() << ", upper bound = " << lp->BestPrimalBound() << ", time elapsed = " << timeElapsed << " milliseconds";
         } else if(LP_STATE == LPVisitorReturnType::Break) {
            auto endTime = std::chrono::steady_clock::now();
            logger->info() << "Optimization took " << std::chrono::duration_cast<std::chrono::milliseconds>(endTime - beginTime_).count() << " milliseconds and " << curIter_ << " iterations";
            return LPVisitorReturnType::Break;
         } else if(LP_STATE == LPVisitorReturnType::Error) {
            assert(false); // this case is currently not handled
         } else {
            assert(false); // unknown case
         }

         curIter_++;
         remainingIter_--;

         if(LP_STATE == LPVisitorReturnType::ReparametrizeLowerBoundAnisotropic
               || LP_STATE == LPVisitorReturnType::ReparametrizeLowerBoundPrimalAnisotropic
               || LP_STATE == LPVisitorReturnType::ReparametrizeLowerBoundUniform
               || LP_STATE == LPVisitorReturnType::ReparametrizeLowerBoundPrimalUniform) {
            prevDualBound_ = curDualBound_;
            curDualBound_ = lp->BestLowerBound();
         }
         // check if optimization has to be terminated
         if(remainingIter_ == 0) {
            logger->info("One iteration remaining");
            return LPVisitorReturnType::Break;
         } 
         if(lp->BestPrimalBound() <= lp->BestLowerBound() + eps) {
            assert(lp->BestPrimalBound() + eps >= lp->BestLowerBound());
            logger->info() << "Primal cost equals lower bound";
            return LPVisitorReturnType::Break;
         }
         if(timeout_ != std::numeric_limits<REAL>::max() && timeElapsed >= timeout_) {
            logger->info() << "Timeout reached after " << timeElapsed << " seconds";
            remainingIter_ = std::min(INDEX(1),remainingIter_);
         }
         if(maxMemory_ > 0) {
            const INDEX memoryUsed = memory_used()/(1024*1024);
            if(maxMemory_ < memoryUsed) {
               remainingIter_ = std::min(INDEX(1),remainingIter_);
               logger->info() << "Solver uses " << memoryUsed << " MB memory, aborting optimization";
            }
         }
         if(minDualImprovement_ > 0 && curDualBound_ - prevDualBound_ < minDualImprovement_) {
            logger->info() << "Dual improvement smaller than " << minDualImprovement_;
            remainingIter_ = std::min(INDEX(1),remainingIter_);
         }

         if(remainingIter_ == 1) {
            if(roundingReparametrization_ == "anisotropic") {
               return LPVisitorReturnType::ReparametrizeLowerBoundPrimalAnisotropic;
            } else if(roundingReparametrization_ == "uniform") {
               return LPVisitorReturnType::ReparametrizeLowerBoundPrimalUniform;
            }
         }


         // determine next state of solver
         if(curIter_ % primalComputationInterval_ == 0 && curIter_ % lowerBoundComputationInterval_ == 0) {
            if(roundingReparametrization_ == "anisotropic") {
               return LPVisitorReturnType::ReparametrizeLowerBoundPrimalAnisotropic;
            } else if(roundingReparametrization_ == "uniform") {
               return LPVisitorReturnType::ReparametrizeLowerBoundPrimalUniform;
            } else {
               throw std::runtime_error("unknown reparametrization mode");
            }
         } else if(curIter_ % primalComputationInterval_ == 0) {
            if(roundingReparametrization_ == "anisotropic") {
               return LPVisitorReturnType::ReparametrizePrimalAnisotropic;
            } else if(roundingReparametrization_ == "uniform") {
               return LPVisitorReturnType::ReparametrizePrimalUniform;
            } else {
               throw std::runtime_error("unknown reparametrization mode");
            }
         } else if(curIter_ % lowerBoundComputationInterval_ == 0) {
            if(standardReparametrization_ == "anisotropic") {
               return LPVisitorReturnType::ReparametrizeLowerBoundAnisotropic;
            } else if(roundingReparametrization_ == "uniform") {
               return LPVisitorReturnType::ReparametrizeLowerBoundUniform;
            } else {
               throw std::runtime_error("unknown reparametrization mode");
            }
         } else {
            if(standardReparametrization_ == "anisotropic") {
               return LPVisitorReturnType::ReparametrizeAnisotropic;
            } else if(roundingReparametrization_ == "uniform") {
               return LPVisitorReturnType::ReparametrizeUniform;
            } else {
               throw std::runtime_error("unknown reparametrization mode");
            }
         }
         
         assert(false);
      }

      
      using TimeType = decltype(std::chrono::steady_clock::now());
      TimeType GetBeginTime() const { return beginTime_; }
      REAL GetLowerBound() const { return curDualBound_; }
      INDEX GetIter() const { return curIter_; }

      protected:
      // command line arguments TCLAP
      TCLAP::ValueArg<INDEX> maxIterArg_;
      TCLAP::ValueArg<INDEX> maxMemoryArg_;
      TCLAP::ValueArg<INDEX> timeoutArg_;
      //TCLAP::ValueArg<INDEX> boundComputationIntervalArg_;
      TCLAP::ValueArg<INDEX> primalComputationIntervalArg_;
      TCLAP::ValueArg<INDEX> lowerBoundComputationIntervalArg_;
      PositiveRealConstraint posConstraint_;
      TCLAP::ValueArg<REAL> minDualImprovementArg_;
      TCLAP::ValueArg<std::string> standardReparametrizationArg_;
      TCLAP::ValueArg<std::string> roundingReparametrizationArg_;
      TCLAP::ValueArg<std::string> protocolateFileArg_;
      TCLAP::SwitchArg protocolateConsoleArg_;

      // command line arguments read out
      INDEX maxIter_;
      INDEX maxMemory_;
      INDEX timeout_;
      //INDEX boundComputationInterval_;
      INDEX primalComputationInterval_;
      INDEX lowerBoundComputationInterval_;
      REAL minDualImprovement_;
      std::string protocolateFile_;
      // do zrobienia: make enum for reparametrization mode
      std::string standardReparametrization_;
      std::string roundingReparametrization_;
      bool protocolateConsole_;

      // internal state of visitor
      INDEX remainingIter_;
      INDEX curIter_ = 0;
      REAL prevDualBound_ = -std::numeric_limits<REAL>::max();
      REAL curDualBound_ = -std::numeric_limits<REAL>::max();
      TimeType beginTime_;

      PROBLEM_DECOMPOSITION& pd_;
   };

   template<class PROBLEM_DECOMPOSITION>
   class StandardTighteningVisitor : public StandardVisitor<PROBLEM_DECOMPOSITION>
   {
      using BaseVisitorType = StandardVisitor<PROBLEM_DECOMPOSITION>;
      public:
      StandardTighteningVisitor(TCLAP::CmdLine& cmd, PROBLEM_DECOMPOSITION& pd)
         :
            BaseVisitorType(cmd,pd),
            tightenArg_("","tighten","enable tightening",cmd,false),
            tightenIterationArg_("","tightenIteration","number of iterations after which tightening is performed for the first time, default = never",false,std::numeric_limits<INDEX>::max(),"positive integer", cmd),
            tightenIntervalArg_("","tightenInterval","number of iterations between tightenings",false,std::numeric_limits<INDEX>::max(),"positive integer", cmd),
            tightenConstraintsMaxArg_("","tightenConstraintsMax","maximal number of constraints to be added during tightening",false,20,"positive integer"),
            tightenConstraintsPercentageArg_("","tightenConstraintsPercentage","maximal number of constraints to be added during tightening as percentage of number of initial factors",false,0.01,"positive real"),
            posConstraint_(),
            tightenMinDualIncreaseArg_("","tightenMinDualIncrease","minimum increase which additional constraint must guarantee",false,0.0,&posConstraint_, cmd),
            unitIntervalConstraint_(),
            tightenMinDualDecreaseFactorArg_("","tightenMinDualDecreaseFactor","factor by which to decrease minimum dual increase during tightening",false,0.5,&unitIntervalConstraint_, cmd)
      {
         cmd.xorAdd(tightenConstraintsMaxArg_,tightenConstraintsPercentageArg_);
      }

      LPVisitorReturnType begin(const LP* lp) // called, after problem is constructed. 
      {
         try {
            tighten_ = tightenArg_.getValue();
            tightenIteration_ = tightenIterationArg_.getValue();
            tightenInterval_ = tightenIntervalArg_.getValue();
            if(tightenConstraintsPercentageArg_.isSet()) {
               tightenConstraintsPercentage_ = tightenConstraintsPercentageArg_.getValue();
               tightenConstraintsMax_ = INDEX(tightenConstraintsPercentage_ * lp->GetNumberOfFactors());
            } else {
               tightenConstraintsMax_ = tightenConstraintsMaxArg_.getValue();
            }
            tightenMinDualIncrease_ = tightenMinDualIncreaseArg_.getValue();
            tightenMinDualDecreaseFactor_ = tightenMinDualDecreaseFactorArg_.getValue();
         } catch (TCLAP::ArgException &e) {
            std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl; 
            exit(1);
         }

         return BaseVisitorType::begin(lp);
      }

      // the default
      template<LPVisitorReturnType LP_STATE>
      LPVisitorReturnType visit(LP* lp)
      {
         const auto ret = BaseVisitorType::template visit<LP_STATE>(lp);
         auto logger = spdlog::get("logger");
         // do zrobienia: introduce tighten reparametrization
         if(tighten_) {
            if(LP_STATE != LPVisitorReturnType::Break) {
               // do zrobienia: if one specifies tightenIteration = 0, this code does not work. 
               if(this->GetIter() == tightenIteration_ || this->GetIter() >= lastTightenIteration_ + tightenInterval_) {
                  const INDEX noConstraintsAdded = Tighten();
                  if(noConstraintsAdded <= tightenConstraintsMax_/2) {
                     logger->info() << "Added only " << noConstraintsAdded << " constraints";
                     logger->info() << "Decrease minimal dual increase in tightening, current minimal increase = " << tightenMinDualIncrease_;
                     tightenMinDualIncrease_ *= tightenMinDualDecreaseFactor_;
                  }
                  if(noConstraintsAdded <= tightenConstraintsMax_/4) {
                     logger->info() << "Tighten again";
                     Tighten();
                  }
                  lastTightenIteration_ = this->GetIter();
                  lp->Init(); // reinitialize factors, as they might have changed
                  logger->info() << "New number of factors = " << lp->GetNumberOfFactors() << ", tightening took " << tightenTime_ << "ms";
               }
            } else if(LP_STATE == LPVisitorReturnType::Break) {
               logger->info() << "Tightening took " << tightenTime_ << " milliseconds";
            }
         }

         return ret;
      }
      INDEX Tighten()
      {
         auto tightenBeginTime = std::chrono::steady_clock::now();
         const INDEX constraintsAdded = BaseVisitorType::pd_.Tighten(tightenMinDualIncrease_, tightenConstraintsMax_);
         auto tightenEndTime = std::chrono::steady_clock::now();
         tightenTime_ += std::chrono::duration_cast<std::chrono::milliseconds>(tightenEndTime - tightenBeginTime).count();
         return constraintsAdded;
      }

      protected:
      TCLAP::SwitchArg tightenArg_;
      TCLAP::ValueArg<INDEX> tightenIterationArg_; // after how many iterations shall tightening be performed
      TCLAP::ValueArg<INDEX> tightenIntervalArg_; // interval between tightening operations.
      TCLAP::ValueArg<INDEX> tightenConstraintsMaxArg_; // How many constraints to add in tightening maximally
      TCLAP::ValueArg<REAL> tightenConstraintsPercentageArg_; // How many constraints to add in tightening maximally
      PositiveRealConstraint posConstraint_;
      TCLAP::ValueArg<REAL> tightenMinDualIncreaseArg_; // only include constraints which guarantee increase larger than specified value
      OpenUnitIntervalConstraint unitIntervalConstraint_;
      TCLAP::ValueArg<REAL> tightenMinDualDecreaseFactorArg_; 

      bool tighten_;
      bool tightenInNextIteration_ = false;
      bool resumeInNextIteration_ = false;
      LPVisitorReturnType retBeforeTightening_;
      LPReparametrizationMode repamModeBeforeTightening_;
      
      INDEX lastTightenIteration_ = std::numeric_limits<INDEX>::max()/2; // otherwise overflow occurs and tightening start immediately
      INDEX tightenIteration_;
      INDEX tightenInterval_;
      INDEX tightenConstraintsMax_;
      REAL tightenConstraintsPercentage_;
      REAL tightenMinDualIncrease_;
      REAL tightenMinDualDecreaseFactor_;
      
      INDEX tightenTime_ = 0;

   };
} // end namespace LP_MP

#endif // LP_MP_STANDARD_VISITOR_HXX
