#include "evaluate.hxx"
#include "visitors/sqlite_visitor.hxx"
#include "solvers/multicut/multicut.h"

using namespace LP_MP;

// algorithms considered
using FMC_MULTICUT_ODD_CYCLE = FMC_MULTICUT<MessageSendingType::SRMP>;
using FMC_MULTICUT_ODD_CYCLE_ODD_WHEEL = FMC_ODD_WHEEL_MULTICUT<MessageSendingType::SRMP>;



std::string knott_150_prefix = "multicut_datasets/knott-3d-150/";
std::vector<std::string> knott_150_dataset = {
   {knott_150_prefix + "gm_knott_3d_032.h5"},
   {knott_150_prefix + "gm_knott_3d_033.h5"},
   {knott_150_prefix + "gm_knott_3d_034.h5"},
   {knott_150_prefix + "gm_knott_3d_035.h5"},
   {knott_150_prefix + "gm_knott_3d_036.h5"},
   {knott_150_prefix + "gm_knott_3d_037.h5"},
   {knott_150_prefix + "gm_knott_3d_039.h5"}
};

std::string knott_300_prefix = "multicut_datasets/knott-3d-300/";
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

std::string knott_450_prefix = "multicut_datasets/knott-3d-450/";
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

std::string knott_550_prefix = "multicut_datasets/knott-3d-550/";
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

std::string modularity_clustering_prefix = "multicut_datasets/modularity-clustering/";
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
      {"--maxIter"}, {"10000"},
      {"--timeout"}, {"3600"}, // one hour
      {"--minDualImprovement"}, {"0.001"},
      {"--minDualImprovementInterval"}, {"50"},
      {"--lowerBoundComputationInterval"}, {"10"},
      {"--primalComputationInterval"}, {"20"},
      {"--standardReparametrization"}, {"uniform"},
      {"--roundingReparametrization"}, {"uniform"},
      {"--overwriteDbRecord"}, // do zrobienia: possibly deactivate this. Then we do not overwrite
      {"--databaseFile"}, {"multicut.db"},
      {"--tighten"},
      {"--tightenReparametrization"}, {"uniform"},
      {"--tightenIteration"}, {"1"},
      {"--tightenInterval"}, {"5"},
      {"--tightenConstraintsPercentage"}, {"0.1"}
   };

   {
      using FMC = FMC_MULTICUT_ODD_CYCLE;
      using VisitorType = SqliteVisitor<StandardTighteningVisitor>;
      using SolverType = ProblemConstructorRoundingSolver<FMC>;
      static auto Input = MulticutOpenGmInput::ParseProblem<FMC>;
      RunSolver<FMC,VisitorSolver<SolverType,VisitorType>>(Input,knott_150_dataset,options,"knott-3d-150","MPMC-C");
      RunSolver<FMC,VisitorSolver<SolverType,VisitorType>>(Input,knott_300_dataset,options,"knott-3d-300","MPMC-C");
      //RunSolver<FMC,VisitorType,SolverType>(Input,knott_450_dataset,options,"knott-3d-450","MPMC-C");
      //RunSolver<FMC,VisitorType,SolverType>(Input,knott_550_dataset,options,"knott-3d-550","MPMC-C");
      RunSolver<FMC,VisitorSolver<SolverType,VisitorType>>(Input,modularity_clustering_dataset,options,"modularity_clustering","MPMC-C");
   }

   {
      using FMC = FMC_MULTICUT_ODD_CYCLE_ODD_WHEEL;
      using VisitorType = SqliteVisitor<StandardTighteningVisitor>;
      using SolverType = ProblemConstructorRoundingSolver<FMC>;
      static auto Input = MulticutOpenGmInput::ParseProblem<FMC>;
      RunSolver<FMC,VisitorSolver<SolverType,VisitorType>>(Input,knott_150_dataset,options,"knott-3d-150","MPMC-COW");
      RunSolver<FMC,VisitorSolver<SolverType,VisitorType>>(Input,knott_300_dataset,options,"knott-3d-300","MPMC-COW");
      //RunSolver<FMC,VisitorType,SolverType>(Input,knott_450_dataset,options,"knott-3d-450","MPMC-COW");
      //RunSolver<FMC,VisitorType,SolverType>(Input,knott_550_dataset,options,"knott-3d-550","MPMC-COW");
      RunSolver<FMC,VisitorSolver<SolverType,VisitorType>>(Input,modularity_clustering_dataset,options,"modularity_clustering","MPMC-COW");
   }
   return 0;
}

