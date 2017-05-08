#ifndef LP_MP_MULTICUT_H
#define LP_MP_MULTICUT_H

#include "LP_MP.h"
#include "solver.hxx"
#include "factors_messages.hxx"
#include "multicut_factors_messages.hxx"
#include "multiway_cut_factors_messages.hxx"
#include "multicut_triplet_factor.hxx"
#include "factors/constant_factor.hxx"
#include "lifted_multicut_factors_messages.hxx"
#include "multicut_constructor.hxx"

#include "factors/simplex_factor.hxx"
#include "messages/simplex_marginalization_message.hxx"
#include "problem_constructors/mrf_problem_construction.hxx"

#include "parse_rules.h"

#include "hdf5_routines.hxx"


#include "andres/graph/graph.hxx"
#include "andres/graph/grid-graph.hxx"

#include "andres/graph/hdf5/graph.hxx"
#include "andres/graph/hdf5/grid-graph.hxx"
#include "andres/functional.hxx"

#include "andres/graph/multicut/kernighan-lin.hxx"
#include "andres/graph/multicut/greedy-additive.hxx"
#include "andres/graph/multicut-lifted/kernighan-lin.hxx"
#include "andres/graph/multicut-lifted/greedy-additive.hxx"

#include <iostream>
#include <vector>

namespace LP_MP {

// for now, we assume that the model is (re)-constructed each time the functor is called.
// however it might be more efficient to construct the model once in the beginning and
// then only change it according to the internal parametrization when the functor is called
// in this case, the functor would need to have the model as internal state
struct KlRounder {
    
    typedef andres::graph::Graph<> GraphType;

    KlRounder()
    {}

    // TODO do we have to call by value here due to using async or could we also use a call by refernce?
    std::vector<char> operator()(GraphType g, std::vector<REAL> edgeValues) {
   
      std::vector<char> labeling(g.numberOfEdges(), 0);
      if(g.numberOfEdges() > 0) {
         andres::graph::multicut::greedyAdditiveEdgeContraction(g, edgeValues, labeling);
         andres::graph::multicut::kernighanLin(g, edgeValues, labeling, labeling);
      }
      return labeling;

    }

    static std::string name() {
        return "KlRounder";
    }
};

struct LiftedKlRounder {

    typedef andres::graph::Graph<> GraphType;

    LiftedKlRounder()
    {}

    std::vector<char> operator()(
            GraphType originalGraph,
            GraphType liftedGraph,
            std::vector<REAL> edgeValues) {
        
        std::vector<char> labeling(edgeValues.size());
        if(originalGraph.numberOfEdges() > 0) {
           andres::graph::multicut_lifted::greedyAdditiveEdgeContraction(originalGraph, liftedGraph, edgeValues, labeling);
           andres::graph::multicut_lifted::kernighanLin(originalGraph, liftedGraph, edgeValues, labeling, labeling);
        }
        return labeling;
    
    }
    
    // FIXME this is a bit hacky, but we need this overload, s.t. the MulticutConstructor can
    // also be used with LiftedKlRounder
    // TODO figure out if this is actually used for Multicut Rounding in the Lifted Multicut
    // if not, change this to a dummy implementation
    std::vector<char> operator()(GraphType g, std::vector<REAL> edgeValues) {
   
      std::vector<char> labeling(g.numberOfEdges(), 0);
      if(g.numberOfEdges() > 0) {
         andres::graph::multicut::greedyAdditiveEdgeContraction(g, edgeValues, labeling);
         andres::graph::multicut::kernighanLin(g, edgeValues, labeling, labeling);
      }
      return labeling;

    }

    static std::string name() {
        return "LiftedKlRounder";
    }
};

// do zrobienia: possibly rename unary to edge factor

// having the rounder as template with default setting works for the to multicut factor messages,
// but for some reason not for the lifted and multiway-cut factor messages (there having the default value does not work,
// hence the type must be explicitly set in the cpps)

template<MessageSendingType MESSAGE_SENDING, typename ROUNDER>
struct FMC_MULTICUT {
   constexpr static const char* name = "Multicut with cycle constraints";
   //constexpr static MessageSendingType MESSAGE_SENDING = MessageSendingType::SRMP;

   using edge_factor_container = FactorContainer<multicut_edge_factor, FMC_MULTICUT, 0, true>;
   using triplet_factor_container = FactorContainer<multicut_triplet_factor, FMC_MULTICUT, 1>;
   using ConstantFactorContainer = FactorContainer<ConstantFactor, FMC_MULTICUT, 2>;

   using edge_triplet_message_0_container = MessageContainer<multicut_edge_triplet_message_0, 0, 1, variableMessageNumber, 1, FMC_MULTICUT, 0 >;
   using edge_triplet_message_1_container = MessageContainer<multicut_edge_triplet_message_1, 0, 1, variableMessageNumber, 1, FMC_MULTICUT, 1 >;
   using edge_triplet_message_2_container = MessageContainer<multicut_edge_triplet_message_2, 0, 1, variableMessageNumber, 1, FMC_MULTICUT, 2 >;

   using FactorList = meta::list< edge_factor_container, triplet_factor_container, ConstantFactorContainer>;
   using MessageList = meta::list<edge_triplet_message_0_container,edge_triplet_message_1_container,edge_triplet_message_2_container>;

   using multicut = MulticutConstructor<FMC_MULTICUT,0,1,0,1,2,2>;
   using ProblemDecompositionList = meta::list<multicut>;
   using RounderType = ROUNDER;
};

// It would be nice to be able to derive from FMC_MULTICUT. This is not possible due to deviating FMCs. Possibly parametrize above FMC with template
template<MessageSendingType MESSAGE_SENDING, class ROUNDER>
struct FMC_ODD_WHEEL_MULTICUT {
   constexpr static const char* name = "Multicut with cycle and odd wheel constraints";

   using edge_factor_container = FactorContainer<multicut_edge_factor, FMC_ODD_WHEEL_MULTICUT, 0, true>;
   using triplet_factor_container = FactorContainer<multicut_triplet_factor, FMC_ODD_WHEEL_MULTICUT, 1>;
   using odd_3_wheel_factor_container = FactorContainer<multicut_odd_3_wheel_factor, FMC_ODD_WHEEL_MULTICUT, 2>;
   using ConstantFactorContainer = FactorContainer<ConstantFactor, FMC_ODD_WHEEL_MULTICUT, 3>;
      
   using edge_triplet_message_0_container = MessageContainer<multicut_edge_triplet_message_0, 0, 1, variableMessageNumber, 1, FMC_ODD_WHEEL_MULTICUT, 0 >;
   using edge_triplet_message_1_container = MessageContainer<multicut_edge_triplet_message_1, 0, 1, variableMessageNumber, 1, FMC_ODD_WHEEL_MULTICUT, 1 >;
   using edge_triplet_message_2_container = MessageContainer<multicut_edge_triplet_message_2, 0, 1, variableMessageNumber, 1, FMC_ODD_WHEEL_MULTICUT, 2 >;

   using triplet_odd_wheel_message_012 = MessageContainer<multicut_triplet_odd_3_wheel_message_012, 1, 2, variableMessageNumber, 1, FMC_ODD_WHEEL_MULTICUT, 3>;
   using triplet_odd_wheel_message_013 = MessageContainer<multicut_triplet_odd_3_wheel_message_013, 1, 2, variableMessageNumber, 1, FMC_ODD_WHEEL_MULTICUT, 4>;
   using triplet_odd_wheel_message_023 = MessageContainer<multicut_triplet_odd_3_wheel_message_023, 1, 2, variableMessageNumber, 1, FMC_ODD_WHEEL_MULTICUT, 5>;
   using triplet_odd_wheel_message_123 = MessageContainer<multicut_triplet_odd_3_wheel_message_123, 1, 2, variableMessageNumber, 1, FMC_ODD_WHEEL_MULTICUT, 6>;

   using FactorList = meta::list< edge_factor_container, triplet_factor_container, odd_3_wheel_factor_container, ConstantFactorContainer>;
   using MessageList = meta::list<
      edge_triplet_message_0_container, edge_triplet_message_1_container, edge_triplet_message_2_container,  
      triplet_odd_wheel_message_012, triplet_odd_wheel_message_013, triplet_odd_wheel_message_023, triplet_odd_wheel_message_123 
      >;

   using multicut_c = MulticutConstructor<FMC_ODD_WHEEL_MULTICUT,0,1, 0,1,2, 3>;
   using multicut_cow = MulticutOddWheelConstructor<multicut_c,2, 3,4,5,6>;
   using ProblemDecompositionList = meta::list<multicut_cow>;
   using RounderType = ROUNDER;
};


template <typename ROUNDER>
struct FMC_LIFTED_MULTICUT {
   constexpr static const char* name = "Lifted Multicut with cycle constraints";
   constexpr static MessageSendingType MESSAGE_SENDING = MessageSendingType::SRMP;

   // no rounding performed: do it via GAEC and K&L, called from problem constructor
   using edge_factor_container = FactorContainer<multicut_edge_factor, FMC_LIFTED_MULTICUT, 0, true>;
   using triplet_factor_container = FactorContainer<multicut_triplet_factor, FMC_LIFTED_MULTICUT, 1>;
   using cut_factor_container = FactorContainer<LiftedMulticutCutFactor, FMC_LIFTED_MULTICUT, 2>;
   using ConstantFactorContainer = FactorContainer<ConstantFactor, FMC_LIFTED_MULTICUT, 3>;

   using edge_triplet_message_0_container = MessageContainer<multicut_edge_triplet_message_0, 0, 1, variableMessageNumber, 1, FMC_LIFTED_MULTICUT, 0 >;
   using edge_triplet_message_1_container = MessageContainer<multicut_edge_triplet_message_1, 0, 1, variableMessageNumber, 1, FMC_LIFTED_MULTICUT, 1 >;
   using edge_triplet_message_2_container = MessageContainer<multicut_edge_triplet_message_2, 0, 1, variableMessageNumber, 1, FMC_LIFTED_MULTICUT, 2 >;

   typedef MessageContainer<CutEdgeLiftedMulticutFactorMessage, 0, 2, variableMessageNumber, variableMessageNumber, FMC_LIFTED_MULTICUT, 3 > CutEdgeLiftedMulticutFactorMessageContainer;
   typedef MessageContainer<LiftedEdgeLiftedMulticutFactorMessage, 0, 2, variableMessageNumber, variableMessageNumber, FMC_LIFTED_MULTICUT, 4 > LiftedEdgeLiftedMulticutFactorMessageContainer;

   using FactorList = meta::list<
      edge_factor_container, 
      triplet_factor_container, 
      cut_factor_container,
      ConstantFactorContainer 
         >;
   using MessageList = meta::list<
      edge_triplet_message_0_container, edge_triplet_message_1_container, edge_triplet_message_2_container,
      CutEdgeLiftedMulticutFactorMessageContainer, 
      LiftedEdgeLiftedMulticutFactorMessageContainer
         >;

   // TODO is the (non-lifted) rounder actually used ?
   // if yes, it should be included in the lifted rounder struct
   using BaseMulticutConstructor = MulticutConstructor<FMC_LIFTED_MULTICUT,0,1,0,1,2,3>;
   using LiftedMulticutConstructor = class LiftedMulticutConstructor<BaseMulticutConstructor,2,3,4>;
   using ProblemDecompositionList = meta::list<LiftedMulticutConstructor>;
   using RounderType = ROUNDER;

};

// also only separate with violated cycles only in multiway cut
// TODO Do we need a specific rounder here (kernighan lin)? 
template <typename ROUNDER>
struct FMC_MULTIWAY_CUT {
   constexpr static const char* name = "Multiway cut with cycle and odd wheel constraints";

   // multicut
   using edge_factor_container = FactorContainer<multicut_edge_factor, FMC_MULTIWAY_CUT, 0>;
   using triplet_factor_container = FactorContainer<multicut_triplet_factor, FMC_MULTIWAY_CUT, 1>;
   using odd_3_wheel_factor_container = FactorContainer<multicut_odd_3_wheel_factor, FMC_MULTIWAY_CUT, 2>;
   using ConstantFactorContainer = FactorContainer<ConstantFactor, FMC_MULTIWAY_CUT, 3>;

   using edge_triplet_message_0_container = MessageContainer<multicut_edge_triplet_message_0, 0, 1, variableMessageNumber, 1, FMC_MULTIWAY_CUT, 0 >;
   using edge_triplet_message_1_container = MessageContainer<multicut_edge_triplet_message_1, 0, 1, variableMessageNumber, 1, FMC_MULTIWAY_CUT, 1 >;
   using edge_triplet_message_2_container = MessageContainer<multicut_edge_triplet_message_2, 0, 1, variableMessageNumber, 1, FMC_MULTIWAY_CUT, 2 >;

   using triplet_odd_wheel_message_012 = MessageContainer<multicut_triplet_odd_3_wheel_message_012, 1, 2, variableMessageNumber, 1, FMC_MULTIWAY_CUT, 3>;
   using triplet_odd_wheel_message_013 = MessageContainer<multicut_triplet_odd_3_wheel_message_013, 1, 2, variableMessageNumber, 1, FMC_MULTIWAY_CUT, 4>;
   using triplet_odd_wheel_message_023 = MessageContainer<multicut_triplet_odd_3_wheel_message_023, 1, 2, variableMessageNumber, 1, FMC_MULTIWAY_CUT, 5>;
   using triplet_odd_wheel_message_123 = MessageContainer<multicut_triplet_odd_3_wheel_message_123, 1, 2, variableMessageNumber, 1, FMC_MULTIWAY_CUT, 6>;

   // mrf
   using unary_factor_container = FactorContainer<UnarySimplexFactor, FMC_MULTIWAY_CUT, 4, true>;
   using potts_factor_container = FactorContainer<pairwise_potts_factor, FMC_MULTIWAY_CUT, 5>;

   using unary_pairwise_message_0_container = MessageContainer<UnaryPairwiseMessageLeft<MessageSendingType::SRMP,false,true>, 4, 5, variableMessageNumber, 1, FMC_MULTIWAY_CUT, 7>;
   using unary_pairwise_message_1_container = MessageContainer<UnaryPairwiseMessageRight<MessageSendingType::SRMP,false,true>, 4, 5, variableMessageNumber, 1, FMC_MULTIWAY_CUT, 8>;

   // join multicut edge and Potts factor
   using multicut_edge_potts_message_container = MessageContainer<multicut_edge_potts_message, 0, 5, atMostOneMessage, atMostOneMessage, FMC_MULTIWAY_CUT, 9>; 
   // when we tighten, additional edges may not be connected to any MRF factor. Also, before we tighten we actually

   using FactorList = meta::list< 
      edge_factor_container,
      triplet_factor_container,
      odd_3_wheel_factor_container,
      ConstantFactorContainer,

      unary_factor_container,
      potts_factor_container 
         >;
   using MessageList = meta::list<
      edge_triplet_message_0_container, edge_triplet_message_1_container, edge_triplet_message_2_container,  
      triplet_odd_wheel_message_012, triplet_odd_wheel_message_013, triplet_odd_wheel_message_023, triplet_odd_wheel_message_123,

      unary_pairwise_message_0_container, unary_pairwise_message_1_container,
      multicut_edge_potts_message_container 
      >;

   using multicut_c = MulticutConstructor<FMC_MULTIWAY_CUT,0,1, 0,1,2, 3>;
   using multicut_cow = MulticutOddWheelConstructor<multicut_c,2, 3,4,5,6>;
   using mrf = StandardMrfConstructor<FMC_MULTIWAY_CUT, 4, 5, 7, 8>;
   using multiway_cut_c = multiway_cut_constructor<FMC_MULTIWAY_CUT,0,1,9>;
   using ProblemDecompositionList = meta::list<multicut_cow, mrf, multiway_cut_c>; 
   using RounderType = ROUNDER;
};

// TODO Do we need a specific rounder here (kerninghan lin)?
template <typename ROUNDER>
struct FMC_ASYMMETRIC_MULTIWAY_CUT {
   constexpr static const char* name = "Asymmetric multiway cut with cycle and odd wheel constraints";

   // multicut
   using edge_factor_container = FactorContainer<multicut_edge_factor, FMC_ASYMMETRIC_MULTIWAY_CUT, 0>;
   using triplet_factor_container = FactorContainer<multicut_triplet_factor, FMC_ASYMMETRIC_MULTIWAY_CUT, 1>;
   using odd_3_wheel_factor_container = FactorContainer<multicut_odd_3_wheel_factor, FMC_ASYMMETRIC_MULTIWAY_CUT, 2>;
   using ConstantFactorContainer = FactorContainer<ConstantFactor, FMC_ASYMMETRIC_MULTIWAY_CUT, 3>;

   using edge_triplet_message_0_container = MessageContainer<multicut_edge_triplet_message_0, 0, 1, variableMessageNumber, 1, FMC_ASYMMETRIC_MULTIWAY_CUT, 0 >;
   using edge_triplet_message_1_container = MessageContainer<multicut_edge_triplet_message_1, 0, 1, variableMessageNumber, 1, FMC_ASYMMETRIC_MULTIWAY_CUT, 1 >;
   using edge_triplet_message_2_container = MessageContainer<multicut_edge_triplet_message_2, 0, 1, variableMessageNumber, 1, FMC_ASYMMETRIC_MULTIWAY_CUT, 2 >;

   using triplet_odd_wheel_message_012 = MessageContainer<multicut_triplet_odd_3_wheel_message_012, 1, 2, variableMessageNumber, 1, FMC_ASYMMETRIC_MULTIWAY_CUT, 3>;
   using triplet_odd_wheel_message_013 = MessageContainer<multicut_triplet_odd_3_wheel_message_013, 1, 2, variableMessageNumber, 1, FMC_ASYMMETRIC_MULTIWAY_CUT, 4>;
   using triplet_odd_wheel_message_023 = MessageContainer<multicut_triplet_odd_3_wheel_message_023, 1, 2, variableMessageNumber, 1, FMC_ASYMMETRIC_MULTIWAY_CUT, 5>;
   using triplet_odd_wheel_message_123 = MessageContainer<multicut_triplet_odd_3_wheel_message_123, 1, 2, variableMessageNumber, 1, FMC_ASYMMETRIC_MULTIWAY_CUT, 6>;

   // mrf
   using unary_factor_container = FactorContainer<UnarySimplexFactor, FMC_ASYMMETRIC_MULTIWAY_CUT, 4, true>;
   using potts_factor_container = FactorContainer<amwc_pairwise_potts_factor, FMC_ASYMMETRIC_MULTIWAY_CUT, 5>;

   using unary_pairwise_message_0_container = MessageContainer<UnaryPairwiseMessageLeft<MessageSendingType::SRMP,false,true>, 4, 5, variableMessageNumber, 1, FMC_ASYMMETRIC_MULTIWAY_CUT, 7>;
   using unary_pairwise_message_1_container = MessageContainer<UnaryPairwiseMessageRight<MessageSendingType::SRMP,false,true>, 4, 5, variableMessageNumber, 1, FMC_ASYMMETRIC_MULTIWAY_CUT, 8>;

   // join multicut edge and Potts factor
   using multicut_edge_potts_message_container = MessageContainer<multicut_edge_potts_message, 0, 5, atMostOneMessage, atMostOneMessage, FMC_ASYMMETRIC_MULTIWAY_CUT, 9>; 
   // when we tighten, additional edges may not be connected to any MRF factor. Also, before we tighten we actually

   using FactorList = meta::list< 
      edge_factor_container,
      triplet_factor_container,
      odd_3_wheel_factor_container,
      ConstantFactorContainer,

      unary_factor_container,
      potts_factor_container 
         >;
   using MessageList = meta::list<
      edge_triplet_message_0_container, edge_triplet_message_1_container, edge_triplet_message_2_container,  
      triplet_odd_wheel_message_012, triplet_odd_wheel_message_013, triplet_odd_wheel_message_023, triplet_odd_wheel_message_123,

      unary_pairwise_message_0_container, unary_pairwise_message_1_container,
      multicut_edge_potts_message_container 
      >;

   using multicut_c = MulticutConstructor<FMC_ASYMMETRIC_MULTIWAY_CUT,0,1, 0,1,2, 3>;
   using multicut_cow = MulticutOddWheelConstructor<multicut_c,2, 3,4,5,6>;
   using mrf = StandardMrfConstructor<FMC_ASYMMETRIC_MULTIWAY_CUT, 4, 5, 7, 8>;
   using multiway_cut_c = multiway_cut_constructor<FMC_ASYMMETRIC_MULTIWAY_CUT,0,1,9>;
   using ProblemDecompositionList = meta::list<multicut_cow, mrf, multiway_cut_c>; 
   using RounderType = ROUNDER;
};


namespace MulticutOpenGmInput {

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
      /*
      if(!H5Tequal(nativeType, hdf5Type<T>())) {
        H5Dclose(dataset);
        H5Tclose(nativeType);
        H5Tclose(type);
        H5Sclose(filespace);
        throw std::runtime_error("Data types not equal error.");
      }
      */

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
      return out;
   }

   template<typename SOLVER>
   bool ParseProblemOld(const std::string filename, SOLVER& pd)
   {
      //read with hdf5
      auto fileHandle = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
      if(fileHandle < 0) {
         throw std::runtime_error("cannot open file " + filename);
      }
      // first look into dataset "gm".
      auto gmHandle = H5Gopen(fileHandle, "gm", H5P_DEFAULT);
      if(gmHandle < 0) {
         throw std::runtime_error("Could not open HDF5 group gm in " + filename);
      }

      // Read 'numbers-of-states' to get the variable number. Note: this is an array of length equal to number of variables, with each entry having the length as its number.
      // -> numberOfVariables
      INDEX numberOfVariables = GetDimension(gmHandle, "numbers-of-states"); 
      std::cout << "Has " << numberOfVariables << " variables\n";

      // second, read edges in graph in 'factors' in chunks of 5. (${factorNo},1,2,${i1},${i2}), where (i1,i2) is the edge
      // note: we must double check the datatype here
      static_assert(sizeof(size_t) == 8,"HDF5 file format has 64 bits for integers");
      auto factor = ReadVector<size_t>(gmHandle, "factors");
      std::vector<std::tuple<INDEX,INDEX,REAL>> edges;
      for(INDEX i=0; i<factor.size()/5; ++i) {
         const INDEX factorNo = factor[5*i];
         const INDEX one = factor[5*i+1];
         const INDEX two = factor[5*i+2];
         const INDEX i1 = factor[5*i+3];
         const INDEX i2 = factor[5*i+4];
         assert(i1<i2);
         edges.push_back(std::make_tuple(i1,i2,0.0));
      }

      // third, read edge costs in 'function-id-...16006' in chunks of two. (${x1},${x2}) -> the edge cost is x1-x2 plus some global offset, if x1 != 0
      auto functionHandle = H5Gopen(gmHandle, "function-id-16006", H5P_DEFAULT);
      assert(functionHandle >= 0);
      static_assert(sizeof(REAL) == 8, "HDF5 file format has 64 bits for floats");
      auto edgeCosts = ReadVector<REAL>(functionHandle,"values");
      assert(edgeCosts.size()/2 == edges.size());
      auto& mc = pd.template GetProblemConstructor<0>();
      for(INDEX i=0; i<edges.size(); ++i) {
         std::get<2>(edges[i]) = -edgeCosts[2*i] + edgeCosts[2*i+1]; // do zrobienia: use constant factor and add edgeCosts[2*i] to it
         mc.AddToConstant(edgeCosts[2*i]);
         //std::get<2>(edges[i]) = edgeCosts[2*i] - edgeCosts[2*i+1];
      }
      // note: theoretically, this could be much more complicated, if some factors are shared etc. Then this simple approach above will not work.
      
      // fourth, sort the edges by edges. There may be multiple edges (e.g. image-seg dataset). Merge such edges
      std::sort(edges.begin(), edges.end(), [](auto a, auto b) { return std::get<0>(a) == std::get<0>(b) ? std::get<1>(a) < std::get<1>(b) : std::get<0>(a) < std::get<0>(b); });


      for(INDEX i=0; i<edges.size(); ++i) {
         INDEX i1 = std::get<0>(edges[i]); 
         INDEX i2 = std::get<1>(edges[i]); 
         REAL cost = std::get<2>(edges[i]);

         while( i+1 < edges.size() && std::get<0>(edges[i+1]) == i1 && std::get<1>(edges[i+1]) == i2) {
            cost += std::get<2>(edges[i+1]);
            ++i;
         }
         if(i1 > i2) {
            assert(false);
            std::swap(i1,i2);
         }
         //std::cout << "add edge (" << i1 << "," << i2 << ") with cost " << cost << "\n";
         mc.AddUnaryFactor(i1, i2, cost);
      }

      return true;
   }

   // transform a Potts problem into a multiway cut one
   template<typename SOLVER>
   bool ParseProblem(const std::string filename, SOLVER& s)
   {
      auto& mc = s.template GetProblemConstructor<0>();
      return ParseMulticutOpenGM(filename, mc);
   }

   template<typename SOLVER>
   bool ParsePottsProblem(const std::string filename, SOLVER& pd)
   {
      auto& mrf = pd.template GetProblemConstructor<1>();
      const bool ret = ParsePottsGM(filename, mrf);
      assert(ret);

      /*
      if(ret) {
         auto& multicut = pd.template GetProblemConstructor<0>();
         for(INDEX i=0; i<mrf.GetNumberOfPairwiseFactors(); ++i) {
            auto vars = mrf.GetPairwiseVariables(i);
            auto* f = multicut.AddUnaryFactor(std::get<0>(vars), std::get<1>(vars), 0.0);
            auto* m = new typename SOLVER::FMC::multicut_edge_potts_message_container(f, mrf.GetPairwiseFactor(i)); 
            pd.GetLP().AddMessage(m); 
         }
      }
      */

      return ret;
   } 
}

namespace MulticutTextInput {

   using Parsing::opt_whitespace;
   using Parsing::mand_whitespace;
   using Parsing::positive_integer;
   using Parsing::real_number;

   struct init_line : pegtl::seq< opt_whitespace, pegtl::string<'M','U','L','T','I','C','U','T'>, opt_whitespace > {};
   struct numberOfVariables_line : pegtl::seq< opt_whitespace, positive_integer, opt_whitespace > {};
   struct edge_line : pegtl::seq< opt_whitespace, positive_integer, opt_whitespace, positive_integer, opt_whitespace, real_number, opt_whitespace > {};

   struct grammar : pegtl::must<
                    init_line, pegtl::eol,
                    numberOfVariables_line, pegtl::eol,
                    pegtl::star<edge_line, pegtl::eol>,
                    pegtl::opt<edge_line>,
                    pegtl::star<pegtl::sor<mand_whitespace, pegtl::eol>>,
                    pegtl::eof> {};

   struct lifted_line : pegtl::seq<opt_whitespace, pegtl::string<'L','I','F','T','E','D'>, opt_whitespace> {};
   struct lifted_edge_line : pegtl::seq< opt_whitespace, positive_integer, opt_whitespace, positive_integer, opt_whitespace, real_number, opt_whitespace > {};

   struct LiftedMulticutGrammar : pegtl::must<
                    init_line, pegtl::eol,
                    numberOfVariables_line, pegtl::eol,
                    pegtl::star<edge_line, pegtl::eol>,
                    pegtl::opt<edge_line>,
                    pegtl::star<pegtl::sor<mand_whitespace, pegtl::eol>>,
                    // now come lifted edges //
                    lifted_line, pegtl::eol,
                    pegtl::star<lifted_edge_line, pegtl::eol>,
                    pegtl::opt<lifted_edge_line>,
                    pegtl::star<pegtl::sor<mand_whitespace, pegtl::eol>>,
                    pegtl::eof> {};


   // do zrobienia: remove this structure
   struct MulticutInput {
      INDEX numberOfVariables_;
      std::vector<std::tuple<INDEX,INDEX,REAL>> edges_;
   };

   template<typename SOLVER, typename Rule >
      struct action
      : pegtl::nothing< Rule > {};


   template<typename SOLVER> struct action<SOLVER, positive_integer > {
      template<typename INPUT>
      static void apply(const INPUT & in, SOLVER&, std::stack<SIGNED_INDEX>& integer_stack, std::stack<REAL>& real_stack, MulticutInput &)
      {
         integer_stack.push(std::stoul(in.string())); 
      }
   };
   template<typename SOLVER> struct action<SOLVER, real_number > {
      template<typename INPUT>
      static void apply(const INPUT & in, SOLVER&, std::stack<SIGNED_INDEX>& integer_stack, std::stack<REAL>& real_stack, MulticutInput &)
      {
         real_stack.push(std::stod(in.string())); 
      }
   };
   template<typename SOLVER> struct action<SOLVER, numberOfVariables_line > {
      template<typename INPUT>
      static void apply(const INPUT & in, SOLVER&, std::stack<SIGNED_INDEX>& integer_stack, std::stack<REAL>& real_stack, MulticutInput & mcInput)
      {
         assert(integer_stack.size() == 1);
         mcInput.numberOfVariables_ = integer_stack.top();
         integer_stack.pop();
      }
   };
   template<typename SOLVER> struct action<SOLVER, edge_line > {
      template<typename INPUT>
      static void apply(const INPUT & in, SOLVER& pd, std::stack<SIGNED_INDEX>& integer_stack, std::stack<REAL>& real_stack, MulticutInput & mcInput)
      {
         assert(integer_stack.size() == 2);
         assert(real_stack.size() == 1);
         INDEX i2 = integer_stack.top();
         integer_stack.pop();
         INDEX i1 = integer_stack.top();
         integer_stack.pop();
         const REAL cost = real_stack.top();
         real_stack.pop();
         if(i1 > i2) {
            std::swap(i1,i2);
         }
         assert(i1 < i2);
         assert(i2 < mcInput.numberOfVariables_);
         auto& mc = pd.template GetProblemConstructor<0>();
         mc.AddUnaryFactor( i1,i2,cost );
      }
   };
   template<typename SOLVER> struct action<SOLVER, pegtl::eof> {
      template<typename INPUT>
      static void apply(const INPUT & in, SOLVER& pd, std::stack<SIGNED_INDEX>& integer_stack, std::stack<REAL>& real_stack, MulticutInput& mcInput)
      {}
   };

   template<typename SOLVER>
      struct actionSpecialization {
         template<typename RULE> struct type : public action<SOLVER,RULE> {};
      };


   template<typename SOLVER> struct action<SOLVER, lifted_edge_line > {
      template<typename INPUT>
      static void apply(const INPUT & in, SOLVER& pd, std::stack<SIGNED_INDEX>& integer_stack, std::stack<REAL>& real_stack, MulticutInput & mcInput)
      {
         assert(integer_stack.size() == 2);
         assert(real_stack.size() == 1);
         INDEX i2 = integer_stack.top();
         integer_stack.pop();
         INDEX i1 = integer_stack.top();
         integer_stack.pop();
         const REAL cost = real_stack.top();
         real_stack.pop();
         if(i1 > i2) {
            std::swap(i1,i2);
         }
         assert(i1 < i2);
         assert(i2 < mcInput.numberOfVariables_);
         auto& mc = pd.template GetProblemConstructor<0>();
         mc.AddLiftedUnaryFactor( i1,i2,cost );
      }
   };

      
   template<typename SOLVER>
   bool ParseProblem(const std::string filename, SOLVER& pd)
   {
      std::stack<SIGNED_INDEX> integer_stack;
      std::stack<REAL> real_stack;
      MulticutInput mcInput;
      std::cout << "parsing " << filename << "\n";

      pegtl::file_parser problem(filename);

      return problem.parse< grammar, actionSpecialization<SOLVER>::template type >(pd, integer_stack, real_stack, mcInput);
   }

   template<typename SOLVER>
   bool ParseLiftedProblem(const std::string filename, SOLVER& pd)
   {
      std::stack<SIGNED_INDEX> integer_stack;
      std::stack<REAL> real_stack;
      MulticutInput mcInput;
      std::cout << "parsing " << filename << "\n";

      pegtl::file_parser problem(filename);
      return problem.parse< LiftedMulticutGrammar, actionSpecialization<SOLVER>::template type >(pd, integer_stack, real_stack, mcInput);
   }


} // end namespace MulticutTextInput

// HDF5 input as in data of Andres et al.
namespace MulticutH5Input {

   template<typename SOLVER, bool GRID_GRAPH=false>
   bool ParseLiftedProblem(const std::string filename, SOLVER& pd)
   {
      auto& mc = pd.template GetProblemConstructor<0>();
      
      using orig_graph_type = typename std::conditional<GRID_GRAPH, andres::graph::GridGraph<2>, andres::graph::Graph<>>::type;
      orig_graph_type originalGraph;
      andres::graph::Graph<> liftedGraph;
      std::vector<REAL> edgeValues;

      auto fileHandle = andres::graph::hdf5::openFile(filename);
      andres::graph::hdf5::load(fileHandle, "graph", originalGraph);
      andres::graph::hdf5::load(fileHandle, "graph-lifted", liftedGraph);

      std::vector<size_t> shape;
      andres::graph::hdf5::load(fileHandle, "edge-cut-probabilities", shape, edgeValues);
      andres::graph::hdf5::closeFile(fileHandle);

      // transform to energy cost
      std::transform(edgeValues.begin(), edgeValues.end(), edgeValues.begin(), andres::NegativeLogProbabilityRatio<REAL,REAL>());
      assert(edgeValues.size() == liftedGraph.numberOfEdges());

      // FIXME
      // FIXME This looks wrong! Local and lifted edges should be different! 
      // FIXME
      for(std::size_t e=0; e<liftedGraph.numberOfEdges(); ++e) {
         auto i = liftedGraph.vertexOfEdge(e,0);
         auto j = liftedGraph.vertexOfEdge(e,1);
         if(originalGraph.findEdge(i,j).first) {
            mc.AddUnaryFactor( i, j, edgeValues[e] );
         }
      }
      for(std::size_t e=0; e<liftedGraph.numberOfEdges(); ++e) {
         auto i = liftedGraph.vertexOfEdge(e,0);
         auto j = liftedGraph.vertexOfEdge(e,1);
         if(!originalGraph.findEdge(i,j).first) {
            mc.AddLiftedUnaryFactor( i, j, edgeValues[e] );
         }
      }
      return true;
   }
} // end namespace MulticutH5Input

} // end namespace LP_MP

#endif // LP_MP_MULTICUT_H
