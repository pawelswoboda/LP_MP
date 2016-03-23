#ifndef LP_MP_TIKZ_VISITOR_HXX
#define LP_MP_TIKZ_VISITOR_HXX

#include "standard_visitor.hxx"
#include <iostream>
#include <fstream>

#include "solvers/evaluation/evaluate.hxx" // do zrobienia: for extract file name. Put this in special file

namespace LP_MP {

struct IterationStatistics {
   INDEX iteration_;
   INDEX timeElapsed_; // in milliseconds
   REAL lowerBound_;
   REAL upperBound_;
};

// visitor deriving from StandardVisitor and adding option to output tikz graphs for runtime and convergence plots 
// rename to PlotsVisitor
// this visitor generates tikz plots and a file with runtime information
template<class PROBLEM_DECOMPOSITION>
class TikzVisitor : public StandardVisitor<PROBLEM_DECOMPOSITION> {

public:
   TikzVisitor(TCLAP::CmdLine& cmd, PROBLEM_DECOMPOSITION& pd) 
      :
         StandardVisitor<PROBLEM_DECOMPOSITION>(cmd,pd),
         pd_(pd),
         runtimePlotFileArg_("","runtimePlotFile","file into which to protocolate runtime vs. primal/dual energy",false,"","file name",cmd),
         iterationPlotFileArg_("","iterationPlotFile","file into which to protocolate iteration vs. primal/dual energy",false,"","file name",cmd),
         upperBoundColorArg_("","upperBoundColor","color for upper bound",false,"red","color",cmd),
         lowerBoundColorArg_("","lowerBoundColor","color for lower bound",false,"blue","color",cmd),
         upperBoundLegendArg_("","upperBoundLegend","legend entry for upper bound",false,"","text",cmd),
         lowerBoundLegendArg_("","lowerBoundLegend","legend entry for lower bound",false,"","text",cmd),
         statisticsFileArg_("","statisticsFile","file name into which to write iteration/runtime statistic",false,"","file name",cmd)
   {}

   LPVisitorReturnType begin(const LP* lp) // called, after problem is constructed. 
   {
      try {
         runtimePlotFile_ = runtimePlotFileArg_.getValue();
         iterationPlotFile_ = iterationPlotFileArg_.getValue();
         upperBoundColor_ = upperBoundColorArg_.getValue();
         lowerBoundColor_ = lowerBoundColorArg_.getValue();
         upperBoundLegend_ = upperBoundLegendArg_.getValue();
         lowerBoundLegend_ = lowerBoundLegendArg_.getValue();
         statisticsFile_ = statisticsFileArg_.getValue();
      } catch (TCLAP::ArgException &e) {
         std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl; 
         exit(1);
      }
      return StandardVisitor<PROBLEM_DECOMPOSITION>::begin(lp);
   }

   template<LPVisitorReturnType LP_STATE>
   LPVisitorReturnType visit(const LP* lp)
   {
      auto ret_state = this->template StandardVisitor<PROBLEM_DECOMPOSITION>::template visit<LP_STATE>(lp);
      switch(LP_STATE) {
         case LPVisitorReturnType::ReparametrizeAndComputePrimal:
            {
            const REAL lowerBound = this->GetLowerBound();
            const REAL upperBound = lp->BestPrimalBound();
            const INDEX timeElapsed = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - StandardVisitor<PROBLEM_DECOMPOSITION>::GetBeginTime()).count();
            const INDEX curIter = StandardVisitor<PROBLEM_DECOMPOSITION>::GetIter();
            iterationStatistics_.push_back({curIter,timeElapsed,lowerBound,upperBound});
            break;
            }
         case LPVisitorReturnType::Break:
            {
            PlotBounds(iterationPlotFile_, lowerBoundColor_, upperBoundColor_, lowerBoundLegend_, upperBoundLegend_,
                  iterationStatistics_, [](const IterationStatistics& stat) -> REAL { return stat.iteration_; });
            PlotBounds(runtimePlotFile_, lowerBoundColor_, upperBoundColor_, lowerBoundLegend_, upperBoundLegend_,
                  iterationStatistics_, [](const IterationStatistics& stat) -> REAL { return REAL(stat.timeElapsed_)/1000.0; });
            WriteBounds(statisticsFile_,ExtractFilename(pd_.GetInputFileName()),PROBLEM_DECOMPOSITION::FactorMessageConnection::name, iterationStatistics_);
            break;
            }
      }
      return ret_state;
   }

   // write bounds into file
   void WriteBounds(const std::string& filename, const std::string datasetName, const std::string algorithmName, const std::vector<IterationStatistics>& iterStats)
   {
      if(!filename.empty()) {
         std::ofstream bounds({filename}, std::ios_base::app);
         if(!bounds.is_open()) { 
            throw std::runtime_error("Could not open output file for bounds " + filename);
         }
         bounds << "dataset name " << datasetName << "\n";
         bounds << "Algorithm " << algorithmName << "\n";
         for(auto it : iterStats) {
            bounds << it.iteration_ << ", " << it.timeElapsed_ << ", " << it.lowerBound_ << ", " << it.upperBound_ << "\n";
         }
         bounds.close();
      }
   }

   template<typename LAMBDA>
   void PlotBounds(const std::string filename, 
         const std::string& dualColor, const std::string& primalColor, 
         const std::string& lowerBoundLegend, const std::string& upperBoundLegend,
         const std::vector<IterationStatistics>& iterStats, LAMBDA horizontalAxisGetter)
   {
      if(!filename.empty()) {
         std::ofstream plot({filename}, std::ios_base::app);
         if(!plot.is_open()) { 
            throw std::runtime_error("Could not open output file for tikz plot " + filename);
         }
         if(iterStats.size() < 1001) {
            plot << "\\addplot[thick," << dualColor;
            if(lowerBoundLegend.empty()) {
               plot << ", forget plot";
            }
            plot << "] plot coordinates {\n";
         } else {
            plot << "\\addplot[thick," << dualColor << ",each nth point=" << iterationStatistics_.size()/1000 << ", unbounded coords=discard] plot coordinates {\n";
         }

         for(auto it : iterStats) {
            plot << "(" << horizontalAxisGetter(it) << "," << it.lowerBound_ << ")\n";
         }
         plot << "};\n";
         if(!lowerBoundLegend.empty()) {
            plot << "\\addlegendentry{" << lowerBoundLegend << "}\n\n";
         }

         if(iterStats.back().upperBound_ != std::numeric_limits<REAL>::max()) {
            plot << "\\addplot[thick, " << primalColor;
            if(upperBoundLegend.empty()) {
               plot << ", forget plot";
            }
            plot << "] plot coordinates {\n";
            for(INDEX i=0; i<iterStats.size(); ++i) {
               //if(iterStats[i].upperBound_ != std::numeric_limits<REAL>::max()) { // upper bound need not be feasible
               // heuristic upper bound.
               if(iterStats[i].upperBound_ < 1e10) { // upper bound need not be feasible
                  if(i==0 || i==iterStats.size()-1 ||
                        iterationStatistics_[i-1].upperBound_ > iterationStatistics_[i].upperBound_ ||
                        i < iterStats[i].upperBound_ > iterStats[i+1].upperBound_) {
                     plot << "(" << horizontalAxisGetter(iterStats[i]) << "," << iterStats[i].upperBound_ << ")\n";
                  }
               }
            }
            plot << "};\n";
            //if(!upperBoundLegend.empty()) {
            //   plot << "\\addlegendentry{" << upperBoundLegend << "}\n\n";
            //}
         }
         plot.close();
      }
   }

private:
   TCLAP::ValueArg<std::string> runtimePlotFileArg_;
   TCLAP::ValueArg<std::string> iterationPlotFileArg_;
   TCLAP::ValueArg<std::string> upperBoundColorArg_;
   TCLAP::ValueArg<std::string> lowerBoundColorArg_;
   TCLAP::ValueArg<std::string> upperBoundLegendArg_;
   TCLAP::ValueArg<std::string> lowerBoundLegendArg_;
   TCLAP::ValueArg<std::string> statisticsFileArg_;

   std::string runtimePlotFile_;
   std::string iterationPlotFile_;
   std::string upperBoundColor_;
   std::string lowerBoundColor_;
   std::string upperBoundLegend_;
   std::string lowerBoundLegend_;
   std::string statisticsFile_;

   std::vector<IterationStatistics> iterationStatistics_;
   PROBLEM_DECOMPOSITION& pd_;
};

} // end namespace LP_MP

#endif // LP_MP_TIKZ_VISITOR_HXX
