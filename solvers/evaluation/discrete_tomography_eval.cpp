#include "evaluate.hxx"
#include "visitors/sqlite_visitor.hxx"
#include "solvers/discrete_tomography/discrete_tomography.h"

using namespace LP_MP;

std::string dt_synthetic_prefix = "discrete_tomography_datasets/discrete_tomography_synthetic/mp";

std::string dt_synthetic_prefix_2_1 = "/2/1/";
std::vector<std::string> dt_synthetic_2_1 = {
   {dt_synthetic_prefix_2_1 + "0.10_0.10_100.txt"},
   {dt_synthetic_prefix_2_1 + "0.10_0.10_10.txt"},
   {dt_synthetic_prefix_2_1 + "0.10_0.10_11.txt"},
   {dt_synthetic_prefix_2_1 + "0.10_0.10_12.txt"},
   {dt_synthetic_prefix_2_1 + "0.10_0.10_13.txt"},
   {dt_synthetic_prefix_2_1 + "0.10_0.10_14.txt"},
   {dt_synthetic_prefix_2_1 + "0.10_0.10_15.txt"},
   {dt_synthetic_prefix_2_1 + "0.10_0.10_16.txt"},
   {dt_synthetic_prefix_2_1 + "0.10_0.10_17.txt"},
   {dt_synthetic_prefix_2_1 + "0.10_0.10_18.txt"},
   {dt_synthetic_prefix_2_1 + "0.10_0.10_19.txt"},
   {dt_synthetic_prefix_2_1 + "0.10_0.10_1.txt"},
   {dt_synthetic_prefix_2_1 + "0.10_0.10_20.txt"},
   {dt_synthetic_prefix_2_1 + "0.10_0.10_21.txt"},
   {dt_synthetic_prefix_2_1 + "0.10_0.10_22.txt"},
   {dt_synthetic_prefix_2_1 + "0.10_0.10_23.txt"},
   {dt_synthetic_prefix_2_1 + "0.10_0.10_24.txt"},
   {dt_synthetic_prefix_2_1 + "0.10_0.10_25.txt"},
   {dt_synthetic_prefix_2_1 + "0.10_0.10_26.txt"},
   {dt_synthetic_prefix_2_1 + "0.10_0.10_27.txt"},
   {dt_synthetic_prefix_2_1 + "0.10_0.10_28.txt"},
   {dt_synthetic_prefix_2_1 + "0.10_0.10_29.txt"},
   {dt_synthetic_prefix_2_1 + "0.10_0.10_2.txt"},
   {dt_synthetic_prefix_2_1 + "0.10_0.10_30.txt"},
   {dt_synthetic_prefix_2_1 + "0.10_0.10_31.txt"},
   {dt_synthetic_prefix_2_1 + "0.10_0.10_32.txt"},
   {dt_synthetic_prefix_2_1 + "0.10_0.10_33.txt"},
   {dt_synthetic_prefix_2_1 + "0.10_0.10_34.txt"},
   {dt_synthetic_prefix_2_1 + "0.10_0.10_35.txt"},
   {dt_synthetic_prefix_2_1 + "0.10_0.10_36.txt"},
   {dt_synthetic_prefix_2_1 + "0.10_0.10_37.txt"},
   {dt_synthetic_prefix_2_1 + "0.10_0.10_38.txt"},
   {dt_synthetic_prefix_2_1 + "0.10_0.10_39.txt"},
   {dt_synthetic_prefix_2_1 + "0.10_0.10_3.txt"},
   {dt_synthetic_prefix_2_1 + "0.10_0.10_40.txt"},
   {dt_synthetic_prefix_2_1 + "0.10_0.10_41.txt"},
   {dt_synthetic_prefix_2_1 + "0.10_0.10_42.txt"},
   {dt_synthetic_prefix_2_1 + "0.10_0.10_43.txt"},
   {dt_synthetic_prefix_2_1 + "0.10_0.10_44.txt"},
   {dt_synthetic_prefix_2_1 + "0.10_0.10_45.txt"},
   {dt_synthetic_prefix_2_1 + "0.10_0.10_46.txt"},
   {dt_synthetic_prefix_2_1 + "0.10_0.10_47.txt"},
   {dt_synthetic_prefix_2_1 + "0.10_0.10_48.txt"},
   {dt_synthetic_prefix_2_1 + "0.10_0.10_49.txt"},
   {dt_synthetic_prefix_2_1 + "0.10_0.10_4.txt"},
   {dt_synthetic_prefix_2_1 + "0.10_0.10_50.txt"},
   {dt_synthetic_prefix_2_1 + "0.10_0.10_51.txt"},
   {dt_synthetic_prefix_2_1 + "0.10_0.10_52.txt"},
   {dt_synthetic_prefix_2_1 + "0.10_0.10_53.txt"},
   {dt_synthetic_prefix_2_1 + "0.10_0.10_54.txt"},
   {dt_synthetic_prefix_2_1 + "0.10_0.10_55.txt"},
   {dt_synthetic_prefix_2_1 + "0.10_0.10_56.txt"},
   {dt_synthetic_prefix_2_1 + "0.10_0.10_57.txt"},
   {dt_synthetic_prefix_2_1 + "0.10_0.10_58.txt"},
   {dt_synthetic_prefix_2_1 + "0.10_0.10_59.txt"},
   {dt_synthetic_prefix_2_1 + "0.10_0.10_5.txt"},
   {dt_synthetic_prefix_2_1 + "0.10_0.10_60.txt"},
   {dt_synthetic_prefix_2_1 + "0.10_0.10_61.txt"},
   {dt_synthetic_prefix_2_1 + "0.10_0.10_62.txt"},
   {dt_synthetic_prefix_2_1 + "0.10_0.10_63.txt"},
   {dt_synthetic_prefix_2_1 + "0.10_0.10_64.txt"},
   {dt_synthetic_prefix_2_1 + "0.10_0.10_65.txt"},
   {dt_synthetic_prefix_2_1 + "0.10_0.10_66.txt"},
   {dt_synthetic_prefix_2_1 + "0.10_0.10_67.txt"},
   {dt_synthetic_prefix_2_1 + "0.10_0.10_68.txt"},
   {dt_synthetic_prefix_2_1 + "0.10_0.10_69.txt"},
   {dt_synthetic_prefix_2_1 + "0.10_0.10_6.txt"},
   {dt_synthetic_prefix_2_1 + "0.10_0.10_70.txt"},
   {dt_synthetic_prefix_2_1 + "0.10_0.10_71.txt"},
   {dt_synthetic_prefix_2_1 + "0.10_0.10_72.txt"},
   {dt_synthetic_prefix_2_1 + "0.10_0.10_73.txt"},
   {dt_synthetic_prefix_2_1 + "0.10_0.10_74.txt"},
   {dt_synthetic_prefix_2_1 + "0.10_0.10_75.txt"},
   {dt_synthetic_prefix_2_1 + "0.10_0.10_76.txt"},
   {dt_synthetic_prefix_2_1 + "0.10_0.10_77.txt"},
   {dt_synthetic_prefix_2_1 + "0.10_0.10_78.txt"},
   {dt_synthetic_prefix_2_1 + "0.10_0.10_79.txt"},
   {dt_synthetic_prefix_2_1 + "0.10_0.10_7.txt"},
   {dt_synthetic_prefix_2_1 + "0.10_0.10_80.txt"},
   {dt_synthetic_prefix_2_1 + "0.10_0.10_81.txt"},
   {dt_synthetic_prefix_2_1 + "0.10_0.10_82.txt"},
   {dt_synthetic_prefix_2_1 + "0.10_0.10_83.txt"},
   {dt_synthetic_prefix_2_1 + "0.10_0.10_84.txt"},
   {dt_synthetic_prefix_2_1 + "0.10_0.10_85.txt"},
   {dt_synthetic_prefix_2_1 + "0.10_0.10_86.txt"},
   {dt_synthetic_prefix_2_1 + "0.10_0.10_87.txt"},
   {dt_synthetic_prefix_2_1 + "0.10_0.10_88.txt"},
   {dt_synthetic_prefix_2_1 + "0.10_0.10_89.txt"},
   {dt_synthetic_prefix_2_1 + "0.10_0.10_8.txt"},
   {dt_synthetic_prefix_2_1 + "0.10_0.10_90.txt"},
   {dt_synthetic_prefix_2_1 + "0.10_0.10_91.txt"},
   {dt_synthetic_prefix_2_1 + "0.10_0.10_92.txt"},
   {dt_synthetic_prefix_2_1 + "0.10_0.10_93.txt"},
   {dt_synthetic_prefix_2_1 + "0.10_0.10_94.txt"},
   {dt_synthetic_prefix_2_1 + "0.10_0.10_95.txt"},
   {dt_synthetic_prefix_2_1 + "0.10_0.10_96.txt"},
   {dt_synthetic_prefix_2_1 + "0.10_0.10_97.txt"},
   {dt_synthetic_prefix_2_1 + "0.10_0.10_98.txt"},
   {dt_synthetic_prefix_2_1 + "0.10_0.10_99.txt"},
   {dt_synthetic_prefix_2_1 + "0.10_0.10_9.txt"}
};

int main()
{
   // emulate command line options
   std::vector<std::string> options = {
      {"--maxIter"}, {"1000"},
      {"--timeout"}, {"3600"}, // one hour
      {"--minDualImprovement"}, {"0.00001"},
      {"--minDualImprovementInterval"}, {"5"},
      {"--lowerBoundComputationInterval"}, {"10"},
      {"--primalComputationInterval"}, {"100000"},
      //{"--overwriteDbRecord"}, // do zrobienia: possibly deactivate this. Then we do not overwrite
      {"--databaseFile"}, {"discrete_tomography.db"}
   };

   using VisitorType = SqliteVisitor<StandardTighteningVisitor>;

   RunSolver<FMC_DT,VisitorType,Solver<FMC_DT>>(DiscreteTomographyTextInput::ParseProblem,dt_synthetic_2_1,options,"2 projections, sparsity 1","MP");
   return 0;
}

