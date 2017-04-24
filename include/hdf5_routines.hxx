#ifndef LP_MP_HDF5_ROUTINES_HXX
#define LP_MP_HDF5_ROUTINES_HXX

#include "LP_MP.h"
#include "hdf5.h"
#include "opengm/opengm.hxx"
#include <opengm/graphicalmodel/graphicalmodel.hxx>
#include <opengm/graphicalmodel/graphicalmodel_hdf5.hxx>
#include <opengm/operations/minimizer.hxx>
#include <opengm/operations/adder.hxx>
#include <opengm/functions/explicit_function.hxx>
#include <opengm/functions/potts.hxx>
#include <opengm/functions/pottsn.hxx>
#include <opengm/functions/pottsg.hxx>
#include "opengm/functions/truncated_absolute_difference.hxx"
#include "opengm/functions/truncated_squared_difference.hxx"

namespace LP_MP {

   /*
   INDEX GetDimension(const hid_t handle, const std::string datasetName) 
   {
      auto dataset = H5Dopen(handle, datasetName.c_str(), H5P_DEFAULT);
      if(dataset < 0) {
         throw std::runtime_error("cannot open dataset " + datasetName);
      }
      hid_t filespace = H5Dget_space(dataset);
      hsize_t dimension = H5Sget_simple_extent_ndims(filespace);
      assert(dimension == 1);
      hsize_t shape[] = {0};// hsize_t[(size_t)(dimension)];
      auto status = H5Sget_simple_extent_dims(filespace, shape, NULL);
      assert(status >= 0);

      H5Dclose(dataset);
      H5Sclose(filespace);

      return shape[0];
   }

   template<typename T>
   std::vector<T> ReadVector(const hid_t handle, const std::string datasetName)
   {
      herr_t status;
      auto dataset = H5Dopen(handle, datasetName.c_str(), H5P_DEFAULT);
      if(dataset < 0) {
         throw std::runtime_error("cannot open dataset " + datasetName);
      }
      hid_t filespace = H5Dget_space(dataset);
      hid_t type = H5Dget_type(dataset);
      hid_t nativeType = H5Tget_native_type(type, H5T_DIR_DESCEND);
      //if(!H5Tequal(nativeType, hdf5Type<T>())) {
      //  H5Dclose(dataset);
      //  H5Tclose(nativeType);
      //  H5Tclose(type);
      //  H5Sclose(filespace);
      //  throw std::runtime_error("Data types not equal error.");
      //}
      

      hsize_t dimension = H5Sget_simple_extent_ndims(filespace);
      assert(dimension == 1);
      hsize_t shape[] = {0};// hsize_t[(size_t)(dimension)];
      status = H5Sget_simple_extent_dims(filespace, shape, NULL);
      assert(status >= 0);
      hid_t memspace = H5Screate_simple(dimension, &shape[0], NULL);

      std::vector<T> out(shape[0]);
      status = H5Dread(dataset, nativeType, memspace, filespace, H5P_DEFAULT, &out[0]);
      H5Dclose(dataset);
      H5Tclose(nativeType);
      H5Tclose(type);
      H5Sclose(memspace);
      H5Sclose(filespace);
      assert(status >= 0);
      return std::move(out);
   }

   template<typename MRF_CONSTRUCTOR>
   bool ParseGMObsolete(const std::string filename, MRF_CONSTRUCTOR& mrf)
   {
      //read with hdf5
      auto fileHandle = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
      if(fileHandle < 0) {
         throw std::runtime_error("Cannot open file " + filename);
      }
      // first look into dataset "gm".
      auto gmHandle = H5Gopen(fileHandle, "gm", H5P_DEFAULT);
      if(gmHandle < 0) {
         throw std::runtime_error("Cannot open HDF5 group gm in " + filename);
      }

      // Read 'numbers-of-states' to get the variable number. Note: this is an array of length equal to number of variables, with each entry having the length as its number.
      // -> numberOfVariables
      INDEX numberOfVariables = GetDimension(gmHandle, "numbers-of-states"); 
      std::cout << "Has " << numberOfVariables << " variables\n";
      //auto statesHandle = H5Gopen(gmHandle, "numbers-of-states", H5P_DEFAULT);
      static_assert(sizeof(size_t) == 8, "");
      auto cardinality = ReadVector<size_t>(gmHandle, "numbers-of-states");

      // second, read nodes and edges (unaries and pairwise) in graph in 'factors'.
      static_assert(sizeof(size_t) == 8,"HDF5 file format has 64 bits for integers");
      auto factorVec = ReadVector<size_t>(gmHandle, "factors");
      INDEX c = 0; // current index in factors
      std::vector<std::vector<INDEX>> factors;
      while(c < factorVec.size()) {
         const INDEX factorNo = factorVec[c];
         const INDEX zero = factorVec[c+1];
         assert(zero == 0);
         const INDEX factorCardinality = factorVec[c+2];
         if(factorCardinality == 1) {
            const INDEX i1 = factorVec[c+3];
            assert(i1 < cardinality.size());
            factors.push_back({factorNo,i1});
            c += 4;
         } else if(factorCardinality == 2) {
            const INDEX i1 = factorVec[c+3];
            const INDEX i2 = factorVec[c+4];
            assert(i1<i2);
            assert(i2 < cardinality.size());
            factors.push_back({factorNo,i1,i2});
            c += 5;
         } else {
            assert(false);
         }
      }

      // third, read edge costs in 'function-id-...16000'.
      auto functionHandle = H5Gopen(gmHandle, "function-id-16000", H5P_DEFAULT);
      assert(functionHandle >= 0);
      static_assert(sizeof(REAL) == 8, "HDF5 file format has 64 bits for floats");
      auto values = ReadVector<REAL>(functionHandle,"values");
      c = 0; // current index in values
      for(const auto& factor : factors) {
         if(factor.size() == 2) {
            const INDEX i1 = factor[1];
            const INDEX c1 = cardinality[i1];
            std::vector<REAL> unaryCost(c1);
            for(INDEX i=0; i<c1; ++i) {
               unaryCost[i] = values[c+i];
            }
            c += c1;
            mrf.AddUnaryFactor(i1,unaryCost);
         } else if(factor.size() == 3) {
            const INDEX i1 = factor[1];
            const INDEX i2 = factor[2];
            const INDEX c1 = cardinality[i1];
            const INDEX c2 = cardinality[i2];
            std::vector<REAL> pairwiseCost(c1*c2);
            for(INDEX i=0; i<c1; ++i) {
               for(INDEX j=0; j<c2; ++j) {
                  pairwiseCost[i + j*c1] = values[c+i+j*c1];
                  //pairwiseCost[i + j*c1] = values[c+i*c2+j]; // which one is correct?
               }
            }
            c += c1*c2;
            mrf.AddPairwiseFactor(i1,i2,pairwiseCost);
         } else {
            assert(false);
         }
      }
      assert(c == values.size());

      return true;
   }
   */

   template<typename MRF_CONSTRUCTOR>
   bool ParseGM(const std::string filename, MRF_CONSTRUCTOR& mrf)
   {
   typedef double ValueType;
   typedef size_t IndexType;
   typedef size_t LabelType;
   typedef opengm::Adder OperatorType;
   typedef opengm::Minimizer AccumulatorType;
   typedef opengm::DiscreteSpace<IndexType, LabelType> SpaceType;

   // Set functions for graphical model
   typedef opengm::meta::TypeListGenerator<
      opengm::ExplicitFunction<ValueType, IndexType, LabelType>,
      opengm::PottsFunction<ValueType, IndexType, LabelType>,
      opengm::PottsNFunction<ValueType, IndexType, LabelType>,
      opengm::PottsGFunction<ValueType, IndexType, LabelType>,
      opengm::TruncatedSquaredDifferenceFunction<ValueType, IndexType, LabelType>,
      opengm::TruncatedAbsoluteDifferenceFunction<ValueType, IndexType, LabelType>
   >::type FunctionTypeList;


   typedef opengm::GraphicalModel<
      ValueType,
      OperatorType,
      FunctionTypeList,
      SpaceType
   > GmType;
   

   GmType gm; 
   opengm::hdf5::load(gm, filename,"gm");

   for(INDEX f=0; f<gm.numberOfFactors(); ++f){

      if(gm[f].numberOfVariables()==0){
         // ignore for now
      }
      else if(gm[f].numberOfVariables()==1){
         const INDEX i = gm.variableOfFactor(f,0);
         std::vector<REAL> unaryCost(gm[f].numberOfLabels(0));         
         for(INDEX l=0; l<gm[f].numberOfLabels(0); ++l){
            unaryCost[l] = gm[f](std::array<INDEX,1>({l}).begin()); 
         } 
         mrf.AddUnaryFactor(i,unaryCost);
      } 
      else if(gm[f].numberOfVariables()==2){
         const INDEX i = gm.variableOfFactor(f,0);
         const INDEX j = gm.variableOfFactor(f,1);
         matrix<REAL> pairwiseCost(gm[f].numberOfLabels(0), gm[f].numberOfLabels(1));//gm[f].size());
         for(INDEX l1=0; l1<gm[f].numberOfLabels(0); ++l1){
            for(INDEX l2=0; l2<gm[f].numberOfLabels(1); ++l2){
               pairwiseCost(l1, l2) = gm[f](std::array<INDEX,2>({l1,l2}).begin()); 
            }
         }
         mrf.AddPairwiseFactor(i,j,pairwiseCost);
      }
      else{
         std::cout << "Factors of order higher than 2 are so far not supported !" <<std::endl;
         return 1;
      }

   }


   return true;
   }


} // end namespace LP_MP

#endif // LP_MP_HDF5_ROUTINES_HXX
