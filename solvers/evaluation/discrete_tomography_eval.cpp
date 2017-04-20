#include "evaluate.hxx"
#include "visitors/postgresql_visitor.hxx"
#include "solvers/discrete_tomography/discrete_tomography.h"
#include "discrete_tomography_eval_problems.h"
//#include "lp_interface/lp_cplex.hxx"
//#include "ConicBundle/cb_interface.hxx"

using namespace LP_MP;

std::mutex UpdateIter;

template<class FMC,class SOLVER,class VisitorType,class I1,class I2,class J1,class J2,class O,class P>
void RunExperiment(std::vector<std::thread>& threads,I1 InstBegin,I2 InstEnd,J1 DatasetBegin,J2 DatasetEnd,O options,std::string algo,P parser){
  for(INDEX i=0;i<threads.size();i++){
    threads[i] =  std::thread([&](){
        INDEX idx = i;
        while(true){
          std::unique_lock<std::mutex> UpdateIterGuard(UpdateIter);
          if( InstBegin != InstEnd ){
            std::string inst = *InstBegin;
            auto sp = *DatasetBegin;
            std::vector<std::string> instv = {inst};

            InstBegin++;
            DatasetBegin++;
            UpdateIterGuard.unlock();
	     
            RunSolver<FMC,VisitorSolver<SOLVER,VisitorType>>
              (parser,instv,options,sp,algo);
          } else { break; }
        }
      });
  }
   
  for(INDEX i=0;i<threads.size();i++){
    if(threads[i].joinable()){
      threads[i].join();
    }
  }
}

int main()
{
  // run synthetic problems on all projection and sparsity levels 
  std::vector<std::string> projections = {{"2"},{"4"}};//,{"6"}};
  std::vector<std::vector<std::string>> sparsities = {
    dt_synthetic_sparsity_1, 
    dt_synthetic_sparsity_2, 
    //dt_synthetic_sparsity_3, 
    //dt_synthetic_sparsity_4, 
    //dt_synthetic_sparsity_5, 
    //dt_synthetic_sparsity_6, 
    //dt_synthetic_sparsity_7, 
    //dt_synthetic_sparsity_8, 
    //dt_synthetic_sparsity_9
  }; 

  int sparsity;
  std::vector<std::string> AllInst;
  std::vector<std::string> AllInstDatasets;
  for(auto p : projections) {
    sparsity = 1;
    for(auto s : sparsities) {
      for(auto i : s) {
        AllInst.push_back("discrete_tomography_datasets/discrete_tomography_synthetic/mp/" + p + "/" + i);
        AllInstDatasets.push_back(p + " projections, sparsity " + std::to_string(sparsity));
      }
      sparsity++;
    }    
  }
   
  /*
  for(auto p : projections) {
    for(auto i : sheep_logan_64){
      AllInst.push_back("discrete_tomography_datasets/discrete_tomography_synthetic/mp/" + p + "/" + i);
      AllInstDatasets.push_back(p + " sheep_logan_64");
    }
  }
  */
  
  using VisitorType = PostgresqlVisitor<StandardTighteningVisitor>;
  std::vector<std::thread> threads(1);
 
  { // MessagePassing_Cplex
    std::vector<std::string> options = {
      {"--maxIter"}, {"10000"},
      {"--timeout"}, {"3600"},
      {"--minDualImprovement"}, {"0.0001"},
      {"--minDualImprovementInterval"}, {"50"},
      {"--standardReparametrization"},{"uniform"},
      {"--lowerBoundComputationInterval"}, {"10"},
      {"--primalComputationInterval"},{"50000"},
      {"--databaseName"}, {"discrete_tomography"},
      {"--databaseAuth"},{"jkuske:tomo"},
      {"--databaseIp"},{"129.206.113.168"}
    };    
    
    using FMC = FMC_DT;
    using SOLVER = Solver<FMC>;
    RunExperiment<FMC,SOLVER,VisitorType>(threads,AllInst.begin(),AllInst.end(),
                                          AllInstDatasets.begin(),AllInstDatasets.end(),
                                          options,"FastMessagePassing_Uniform",
                                          DiscreteTomographyTextInput::ParseProblem<FMC>);
    
  }
    
  // run sheep logan 
  //for(auto& i : sheep_logan) {
  //   i = "discrete_tomography_datasets/discrete_tomography_synthetic/mp/" + i;
  //}
  //RunSolver<FMC_DT,VisitorSolver<Solver<FMC_DT>,VisitorType>>(DiscreteTomographyTextInput::ParseProblem<FMC_DT>,sheep_logan,options,"sheep logan","MP");

  return 0;
}

