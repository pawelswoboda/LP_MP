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
   Visitor(TCLAP::CmdLine& cmd, SOLVER& pd);
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
   // do zrobienia: shall visitor depend on solver?
   //template<class SOLVER>
   class StandardVisitor {
      
      public:
      StandardVisitor(TCLAP::CmdLine& cmd)
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
            primalTime_(0)
      {}

      LpControl begin() // called, after problem is constructed. 
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

            if(standardReparametrizationArg_.getValue() == "anisotropic") {
               standardReparametrization_ = LPReparametrizationMode::Anisotropic;
            } else if(standardReparametrizationArg_.getValue() == "uniform") {
               standardReparametrization_ = LPReparametrizationMode::Uniform;
            } else {
               throw std::runtime_error("reparametrization " + standardReparametrizationArg_.getValue() + " unknown");
            }
            if(roundingReparametrizationArg_.getValue() == "anisotropic") {
               roundingReparametrization_ = LPReparametrizationMode::Anisotropic;
            } else if(roundingReparametrizationArg_.getValue() == "uniform") {
               roundingReparametrization_ = LPReparametrizationMode::Uniform;
            } else {
               throw std::runtime_error("reparametrization " + roundingReparametrizationArg_.getValue() + " unknown");
            }
            protocolateConsole_ = protocolateConsoleArg_.getValue();
            protocolateFile_ = protocolateFileArg_.getValue();
         } catch (TCLAP::ArgException &e) {
            std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl; 
            exit(1);
         }
         //curLowerBound_ = lp->LowerBound();

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


         //if( ! (standardReparametrization_ == "anisotropic" || standardReparametrization_ == "uniform") ) {
         //   throw std::runtime_error("standard repararametrization mode must be {anisotropic|uniform}, is " + standardReparametrization_);
         //}
         //if( ! (roundingReparametrization_ == "anisotropic" || roundingReparametrization_ == "uniform") ) {
         //   throw std::runtime_error("rounding repararametrization mode must be {anisotropic|uniform}, is " + roundingReparametrization_);
         //}

         //spdlog::get("logger")->info() << "Initial number of factors = " << lp->GetNumberOfFactors();
         spdlog::get("logger")->info("Initial lower bound before optimizing = {}", curLowerBound_);
         beginTime_ = std::chrono::steady_clock::now();


         LpControl ret;
         ret.repam = roundingReparametrization_;
         ret.computePrimal = true;
         ret.computeLowerBound = true;
         return ret;
      }

      // LpControl says what was last command to solver, return type gives next
      //template<typename SOLVER>
      LpControl visit(const LpControl c, const REAL lowerBound, const REAL primalBound)
      {
         const INDEX timeElapsed = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - beginTime_).count();
         auto logger = spdlog::get("logger");

         // first output based on what lp solver did in last iteration
         if(c.computePrimal == false && c.computeLowerBound == false) {
            // output nothing
         } else {
            logger->info("iteration = {}{}{}, time elapsed = {}", curIter_, 
                  (c.computeLowerBound ? fmt::format(", lower bound = {}", lowerBound) : ""), 
                  (c.computePrimal ? fmt::format(", upper bound = {}", primalBound) : ""), 
                  timeElapsed);
               //<< (c.computeLowerBound ? fmt::format(", lower bound = {0}", lowerBound) : "")
               //<< (c.computePrimal ? fmt::format(", upper bound = {0}", primalBound) : "")
               //<< ", time elapsed = " << timeElapsed << " milliseconds";
         }
         if(c.end == true) {
            auto endTime = std::chrono::steady_clock::now();
            logger->info("Optimization took {} milliseconds and {} iterations.", std::chrono::duration_cast<std::chrono::milliseconds>(endTime - beginTime_).count(), curIter_);
         } else if(c.error == true) {
            assert(false); // this case is currently not handled
         }

         curIter_++;
         remainingIter_--;

         LpControl ret;

         if(c.computePrimal) {
            prevLowerBound_ = curLowerBound_;
            curLowerBound_ = lowerBound;
         }
         // check if optimization has to be terminated
         if(remainingIter_ == 0) {
            logger->info("One iteration remaining");
            ret.end = true;
            return ret;
         } 
         if(primalBound <= lowerBound + eps) {
            assert(primalBound + eps >= lowerBound);
            logger->info("Primal cost {} greater equal lower bound {}", primalBound, lowerBound);
            ret.end = true;
            return ret;
         }
         if(timeout_ != std::numeric_limits<REAL>::max() && timeElapsed/1000 >= timeout_) {
            logger->info("Timeout reached after {} seconds", timeElapsed);
            remainingIter_ = std::min(INDEX(1),remainingIter_);
         }
         if(maxMemory_ > 0) {
            const INDEX memoryUsed = memory_used()/(1024*1024);
            if(maxMemory_ < memoryUsed) {
               remainingIter_ = std::min(INDEX(1),remainingIter_);
               logger->info("Solver uses {} MB memory, aborting optimization", memoryUsed);
            }
         }
         if(minDualImprovement_ > 0 && curLowerBound_ - prevLowerBound_ < minDualImprovement_) {
            logger->info("Dual improvement smaller than {}", minDualImprovement_);
            remainingIter_ = std::min(INDEX(1),remainingIter_);
         }

         if(remainingIter_ == 1) {
            ret.computePrimal = true;
            ret.computeLowerBound = true;
            ret.repam = roundingReparametrization_;
            return ret;
         }


         // determine next steps of solver
         if(curIter_ % primalComputationInterval_ == 0 && curIter_ % lowerBoundComputationInterval_ == 0) {
            ret.computePrimal = true;
            ret.computeLowerBound = true;
            ret.repam = roundingReparametrization_;
         } else if(curIter_ % primalComputationInterval_ == 0) {
            ret.computePrimal = true;
            ret.repam = roundingReparametrization_;
         } else if(curIter_ % lowerBoundComputationInterval_ == 0) {
            ret.computeLowerBound = true;
            ret.repam = standardReparametrization_;
         } else {
            ret.repam = standardReparametrization_;
         }
         return ret;
      }
      
      using TimeType = decltype(std::chrono::steady_clock::now());
      TimeType GetBeginTime() const { return beginTime_; }
      //`REAL GetLowerBound() const { return curLowerBound_; }
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
      LPReparametrizationMode standardReparametrization_;
      LPReparametrizationMode roundingReparametrization_;
      bool protocolateConsole_;

      // internal state of visitor
      INDEX remainingIter_;
      INDEX curIter_ = 0;
      REAL prevLowerBound_ = -std::numeric_limits<REAL>::max();
      REAL curLowerBound_ = -std::numeric_limits<REAL>::max();
      TimeType beginTime_;

      // primal
      //REAL bestPrimalCost_ = std::numeric_limits<REAL>::infinity();
      //REAL currentPrimalCost_ = std::numeric_limits<REAL>::infinity();

      //REAL bestLowerBound_ = -std::numeric_limits<REAL>::infinity();
      //REAL currentLowerBound_ = -std::numeric_limits<REAL>::infinity();

      //PrimalSolutionStorage currentPrimal_, bestPrimal_;
      INDEX primalTime_;
   };

   //template<class SOLVER>
   class StandardTighteningVisitor : public StandardVisitor //<SOLVER>
   {
      using BaseVisitorType = StandardVisitor; //<SOLVER>;
      public:
      StandardTighteningVisitor(TCLAP::CmdLine& cmd)
         :
            BaseVisitorType(cmd),
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
         cmd.xorAdd(tightenConstraintsMaxArg_,tightenConstraintsPercentageArg_); // do zrobienia: this means that exactly one must be chosen. We want at most one to be chosen
      }

      LpControl begin() // called, after problem is constructed. 
      {
         try {
            tighten_ = tightenArg_.getValue();
            tightenIteration_ = tightenIterationArg_.getValue();
            tightenInterval_ = tightenIntervalArg_.getValue();
            if(tightenConstraintsPercentageArg_.isSet()) {
               assert(false); // not supported due to not having lp information. Possibly change this
               tightenConstraintsPercentage_ = tightenConstraintsPercentageArg_.getValue();
              // tightenConstraintsMax_ = INDEX(tightenConstraintsPercentage_ * lp->GetNumberOfFactors());
            } else if(tightenConstraintsMaxArg_.isSet()) {
               tightenConstraintsMax_ = tightenConstraintsMaxArg_.getValue();
            } else if(tightenArg_.isSet()) {
               throw std::runtime_error("must set number of constraints to add");
            }
            tightenMinDualIncrease_ = tightenMinDualIncreaseArg_.getValue();
            tightenMinDualDecreaseFactor_ = tightenMinDualDecreaseFactorArg_.getValue();
         } catch (TCLAP::ArgException &e) {
            std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl; 
            exit(1);
         }

         return BaseVisitorType::begin();
      }

      // the default
      //template<LPVisitorReturnType LP_STATE>
      LpControl visit(const LpControl c, const REAL lowerBound, const REAL primalBound)
      {
         auto ret = BaseVisitorType::visit(c, lowerBound, primalBound);
         auto logger = spdlog::get("logger");
         // do zrobienia: introduce tighten reparametrization


         if(this->GetIter() == tightenIteration_ || this->GetIter() >= lastTightenIteration_ + tightenInterval_) {
            ret.tighten = true;
            ret.tightenConstraints = tightenConstraintsMax_;
            lastTightenIteration_ = this->GetIter();
         }
         if(this->prevLowerBound_ >= lowerBound - tightenMinDualIncrease_) {
            ret.tighten = true;
            ret.tightenConstraints = tightenConstraintsMax_;
            lastTightenIteration_ = this->GetIter();
         }

         //if(c.end) {
         //   logger->info() << "Tightening took " << tightenTime_ << " milliseconds";
         //}

         return ret;
      }
      /*
      INDEX Tighten()
      {
         auto tightenBeginTime = std::chrono::steady_clock::now();
         const INDEX constraintsAdded = BaseVisitorType::pd_.Tighten(tightenMinDualIncrease_, tightenConstraintsMax_);
         auto tightenEndTime = std::chrono::steady_clock::now();
         tightenTime_ += std::chrono::duration_cast<std::chrono::milliseconds>(tightenEndTime - tightenBeginTime).count();
         return constraintsAdded;
      }
      */

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
