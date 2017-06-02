#include "solvers/evaluation/evaluate.hxx"
#include "visitors/sqlite_visitor.hxx"
#include "solvers/discrete_tomography/discrete_tomography.h"
#include "discrete_tomography_eval_problems.h"

using namespace LP_MP;

int main()
{
   std::vector<std::string> options = {
      {"--maxIter"}, {"1000"},
      //{"--minDualImprovement"}, {"0.0001"},
      //{"--minDualImprovementInterval"}, {"50"},
      {"--standardReparametrization"},{"anisotropic"},
      {"--roundingReparametrization"},{"uniform"},
      {"--lowerBoundComputationInterval"}, {"1"},
      {"--primalComputationInterval"},{"2"},
      {"--overwriteDbRecord"},
      {"--databaseFile"},{"discrete_tomography.db"}
   };    

   using VisitorType = SqliteVisitor<StandardTighteningVisitor>;
   using FMC = FMC_DT;
   using LP_type = LP_sat<LP>;
   using SOLVER = MpRoundingSolver<Solver<FMC,LP_type,VisitorType>>;

   RunSolver<FMC_DT, SOLVER>(DiscreteTomographyTextInput::ParseProblem<Solver<FMC,LP_type,VisitorType>>,get_discrete_tomography_instances(2,1),options,"2-1","MP");
   return 0;
}
