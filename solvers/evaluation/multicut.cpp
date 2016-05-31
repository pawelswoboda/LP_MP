#include "evaluate.hxx"
#include "visitors/sqlite_visitor.hxx"
#include "solvers/multicut/multicut.h"

using namespace LP_MP;

// algorithms considered
using FMC_MULTICUT_ODD_CYCLE = FMC_MULTICUT<MessageSendingType::SRMP>;
using FMC_MULTICUT_ODD_CYCLE_ODD_WHEEL = FMC_ODD_WHEEL_MULTICUT<MessageSendingType::SRMP>;

// OpenGM input
static auto OPENGM_MULTICUT_ODD_CYCLE_INPUT = MulticutOpenGmInput::ParseProblem<FMC_MULTICUT_ODD_CYCLE>;
static auto OPENGM_MULTICUT_ODD_CYCLE_ODD_WHEEL_INPUT = MulticutOpenGmInput::ParseProblem<FMC_MULTICUT_ODD_CYCLE_ODD_WHEEL>;


std::string knott_150_prefix = "../../../solvers/multicut/knott-3d-150/";
std::vector<std::string> knott_150_dataset = {
   {knott_150_prefix + "gm_knott_3d_032.h5"},
   {knott_150_prefix + "gm_knott_3d_033.h5"},
   {knott_150_prefix + "gm_knott_3d_034.h5"},
   {knott_150_prefix + "gm_knott_3d_035.h5"},
   {knott_150_prefix + "gm_knott_3d_036.h5"},
   {knott_150_prefix + "gm_knott_3d_037.h5"},
   {knott_150_prefix + "gm_knott_3d_039.h5"}
};

std::string knott_300_prefix = "../../../solvers/multicut/knott-3d-300/";
std::vector<std::string> knott_300_dataset = {
   {knott_300_prefix + "gm_knott_3d_072.h5"},
   {knott_300_prefix + "gm_knott_3d_073.h5"},
   {knott_300_prefix + "gm_knott_3d_074.h5"},
   {knott_300_prefix + "gm_knott_3d_075.h5"},
   {knott_300_prefix + "gm_knott_3d_076.h5"},
   {knott_300_prefix + "gm_knott_3d_077.h5"},
   {knott_300_prefix + "gm_knott_3d_078.h5"},
   {knott_300_prefix + "gm_knott_3d_079.h5"}
};

std::string knott_450_prefix = "../../../solvers/multicut/knott-3d-450/";
std::vector<std::string> knott_450_dataset = {
   {knott_450_prefix + "gm_knott_3d_096.h5"},
   {knott_450_prefix + "gm_knott_3d_097.h5"},
   {knott_450_prefix + "gm_knott_3d_098.h5"},
   {knott_450_prefix + "gm_knott_3d_099.h5"},
   {knott_450_prefix + "gm_knott_3d_100.h5"},
   {knott_450_prefix + "gm_knott_3d_101.h5"},
   {knott_450_prefix + "gm_knott_3d_102.h5"},
   {knott_450_prefix + "gm_knott_3d_103.h5"}
};

std::string knott_550_prefix = "../../../solvers/multicut/knott-3d-550/";
std::vector<std::string> knott_550_dataset = {
   {knott_550_prefix + "gm_knott_3d_112.h5"},
   {knott_550_prefix + "gm_knott_3d_113.h5"},
   {knott_550_prefix + "gm_knott_3d_114.h5"},
   {knott_550_prefix + "gm_knott_3d_115.h5"},
   {knott_550_prefix + "gm_knott_3d_116.h5"},
   {knott_550_prefix + "gm_knott_3d_117.h5"},
   {knott_550_prefix + "gm_knott_3d_118.h5"},
   {knott_550_prefix + "gm_knott_3d_119.h5"}
};

std::string modularity_clustering_prefix = "../../../solvers/multicut/modularity-clustering/";
std::vector<std::string> modularity_clustering_dataset = {
   {modularity_clustering_prefix + "adjnoun.h5"},
   {modularity_clustering_prefix + "dolphins.h5"},
   {modularity_clustering_prefix + "football.h5"},
   {modularity_clustering_prefix + "karate.h5"},
   {modularity_clustering_prefix + "lesmis.h5"},
   {modularity_clustering_prefix + "polbooks.h5"}
};


int main()
{
   // emulate command line options
   std::vector<std::string> options = {
      {"--maxIter"}, {"1"},
      {"--timeout"}, {"3600"},
      //{"--minDualImprovement"}, {"0.00001"},
      {"--lowerBoundComputationInterval"}, {"1"},
      {"--primalComputationInterval"}, {"20"},
      {"--standardReparametrization"}, {"uniform"},
      {"--roundingReparametrization"}, {"uniform"},
      //{"--overwriteDbRecord"}, // do zrobienia: possibly deactivate this. Then we do not overwrite
      {"--databaseFile"}, {"multicut.db"},
      {"--tighten"},
      {"--tightenIteration"}, {"1"},
      {"--tightenInterval"}, {"30"},
      {"--tightenConstraintsPercentage"}, {"0.01"},
      {"--tightenMinDualIncrease"}, {"1"},
      {"--tightenMinDualDecreaseFactor"}, {"0.5"}
   };

   using FMC_MULTICUT_ODD_CYCLE_VISITOR = SqliteVisitor<ProblemDecomposition<FMC_MULTICUT_ODD_CYCLE>,StandardTighteningVisitor>;
   RunSolver<FMC_MULTICUT_ODD_CYCLE,FMC_MULTICUT_ODD_CYCLE_VISITOR>(OPENGM_MULTICUT_ODD_CYCLE_INPUT,knott_150_dataset,options,"knott-3d-150","MPMC-OC");
   RunSolver<FMC_MULTICUT_ODD_CYCLE,FMC_MULTICUT_ODD_CYCLE_VISITOR>(OPENGM_MULTICUT_ODD_CYCLE_INPUT,knott_300_dataset,options,"knott-3d-300","MPMC-OC");
   RunSolver<FMC_MULTICUT_ODD_CYCLE,FMC_MULTICUT_ODD_CYCLE_VISITOR>(OPENGM_MULTICUT_ODD_CYCLE_INPUT,knott_450_dataset,options,"knott-3d-450","MPMC-OC");
   RunSolver<FMC_MULTICUT_ODD_CYCLE,FMC_MULTICUT_ODD_CYCLE_VISITOR>(OPENGM_MULTICUT_ODD_CYCLE_INPUT,knott_550_dataset,options,"knott-3d-550","MPMC-OC");
   RunSolver<FMC_MULTICUT_ODD_CYCLE,FMC_MULTICUT_ODD_CYCLE_VISITOR>(OPENGM_MULTICUT_ODD_CYCLE_INPUT,modularity_clustering_dataset,options,"modularity_clustering","MPMC-OC");

   using FMC_MULTICUT_ODD_CYCLE_ODD_WHEEL_VISITOR = SqliteVisitor<ProblemDecomposition<FMC_MULTICUT_ODD_CYCLE_ODD_WHEEL>,StandardTighteningVisitor>;
   RunSolver<FMC_MULTICUT_ODD_CYCLE_ODD_WHEEL,FMC_MULTICUT_ODD_CYCLE_ODD_WHEEL_VISITOR>(OPENGM_MULTICUT_ODD_CYCLE_ODD_WHEEL_INPUT,knott_150_dataset,options,"knott-3d-150","MPMC-OCOW");
   RunSolver<FMC_MULTICUT_ODD_CYCLE_ODD_WHEEL,FMC_MULTICUT_ODD_CYCLE_ODD_WHEEL_VISITOR>(OPENGM_MULTICUT_ODD_CYCLE_ODD_WHEEL_INPUT,knott_300_dataset,options,"knott-3d-300","MPMC-OCOW");
   RunSolver<FMC_MULTICUT_ODD_CYCLE_ODD_WHEEL,FMC_MULTICUT_ODD_CYCLE_ODD_WHEEL_VISITOR>(OPENGM_MULTICUT_ODD_CYCLE_ODD_WHEEL_INPUT,knott_450_dataset,options,"knott-3d-450","MPMC-OCOW");
   RunSolver<FMC_MULTICUT_ODD_CYCLE_ODD_WHEEL,FMC_MULTICUT_ODD_CYCLE_ODD_WHEEL_VISITOR>(OPENGM_MULTICUT_ODD_CYCLE_ODD_WHEEL_INPUT,knott_550_dataset,options,"knott-3d-550","MPMC-OCOW");
   RunSolver<FMC_MULTICUT_ODD_CYCLE_ODD_WHEEL,FMC_MULTICUT_ODD_CYCLE_ODD_WHEEL_VISITOR>(OPENGM_MULTICUT_ODD_CYCLE_ODD_WHEEL_INPUT,modularity_clustering_dataset,options,"modularity_clustering","MPMC-OCOW");


   //FinishTikzFiles(graphMatchingDatasets);
}

