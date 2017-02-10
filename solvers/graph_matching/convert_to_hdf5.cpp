#include "solvers/graph_matching/graph_matching.h"
#include "opengm/opengm.hxx"
#include "opengm/graphicalmodel/graphicalmodel.hxx"
#include "opengm/graphicalmodel/graphicalmodel_hdf5.hxx"
#include "opengm/operations/adder.hxx"
#include "assert.h"

int main(int argc, char** argv)
{
   assert(argc == 3); // first arg input in torresani format, second is output in opengm format

   auto gm_input = TorresaniEtAlInput::ParseFile(std::string(argv[1]));

   typedef double                                                               ValueType;          // type used for values
   typedef size_t                                                               IndexType;          // type used for indexing nodes and factors (default : size_t)
   typedef size_t                                                               LabelType;          // type used for labels (default : size_t)
   typedef opengm::Adder                                                        OpType;             // operation used to combine terms
   typedef opengm::ExplicitFunction<ValueType,IndexType,LabelType>              ExplicitFunction;   // shortcut for explicite function
   typedef opengm::meta::TypeListGenerator<ExplicitFunction>::type              FunctionTypeList;   // list of all function the model cal use (this trick avoids virtual methods) - here only one
   typedef opengm::DiscreteSpace<IndexType, LabelType>                          SpaceType;          // type used to define the feasible statespace
   typedef opengm::GraphicalModel<ValueType,OpType,FunctionTypeList,SpaceType>  Model;              // type of the model
   typedef Model::FunctionIdentifier                                            FunctionIdentifier; // type of the function identifier


   //*******************
   //** Code
   //*******************

   std::cout << "Start building the model ... "<<std::endl;

   auto left_unaries = TorresaniEtAlInput::build_left_unaries(gm_input, 1.0);
   std::vector<LabelType> numbersOfLabels;
   for(const auto& v : left_unaries) {
      numbersOfLabels.push_back(v.size());
   }

   // Build empty Model
   Model gm(SpaceType(numbersOfLabels.begin(), numbersOfLabels.end()));

   // Add 1st order functions and factors to the model
   for(IndexType variable = 0; variable < gm.numberOfVariables(); ++variable) {
      // construct 1st order function
      const LabelType shape[] = {gm.numberOfLabels(variable)};
      ExplicitFunction f(shape, shape + 1);
      for(LabelType state = 0; state < gm.numberOfLabels(variable); ++state) { 
         f(&state) = left_unaries[variable][state];
         //f(state) = ValueType(rand()) / RAND_MAX; // only works for ExpliciteFunction
      }
      // add function
      FunctionIdentifier id = gm.addFunction(f);
      // add factor
      IndexType variableIndex[] = {variable};
      gm.addFactor(id, variableIndex, variableIndex + 1);
   }
   // add 2nd order function and factors to the model
   auto quadr_pot = BuildLeftPairwisePotentials(gm_input, 1.0);
   for(auto& q : quadr_pot) {
      IndexType vars[] = {q.first.first, q.first.second};
      LabelType shape[] = {numbersOfLabels[vars[0]],numbersOfLabels[vars[1]]};
      LabelType state[] = {0,0};
      ExplicitFunction f(shape, shape + 2);
      for(state[0] = 0; state[0] < gm.numberOfLabels(vars[0]); ++state[0]){
         for(state[1] = 0; state[1] < gm.numberOfLabels(vars[1]); ++state[1]) {
            //f(state) = q.second[state[0]*gm.numberOfLabels(vars[1]) + state[1]];
            f(state) = q.second[state[0] + state[1]*gm.numberOfLabels(vars[0])];
         }
      }
      FunctionIdentifier fid = gm.addFunction(f);
      // add factor
      gm.addFactor(fid, vars, vars + 2);
   }

   // View some model information
   std::cout << "The model has " << gm.numberOfVariables() << " variables."<<std::endl;
   for(size_t i=0; i<gm.numberOfVariables(); ++i){
      std::cout << " * Variable " << i << " has "<< gm.numberOfLabels(i) << " labels."<<std::endl; 
   } 
   std::cout << "The model has " << gm.numberOfFactors() << " factors."<<std::endl;
   for(size_t f=0; f<gm.numberOfFactors(); ++f){
      std::cout << " * Factor " << f << " has order "<< gm[f].numberOfVariables() << "."<<std::endl; 
   }

   opengm::hdf5::save(gm, argv[2],"gm");

   return 0;
}
