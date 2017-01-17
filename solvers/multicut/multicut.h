#ifndef LP_MP_MULTICUT_HXX
#define LP_MP_MULTICUT_HXX

#include "LP_MP.h"
#include "solver.hxx"
#include "factors_messages.hxx"
#include "multicut_unary_factor.hxx"
#include "multicut_triplet_factor.hxx"
#include "multicut_odd_wheel.hxx"
#include "factors/constant_factor.hxx"
#include "lifted_multicut_factors_messages.hxx"
#include "multicut_constructor.hxx"

#include "parse_rules.h"

#include "hdf5_routines.hxx"

#include "andres/graph/graph.hxx"
#include "andres/graph/grid-graph.hxx"
#include "andres/graph/hdf5/graph.hxx"
#include "andres/graph/hdf5/grid-graph.hxx"
#include "andres/functional.hxx"
#include "andres/graph/multicut-lifted/kernighan-lin.hxx"

#include <iostream>
#include <vector>

//#include "boost/hana.hpp"

namespace LP_MP {

//namespace hana = boost::hana;

// do zrobienia: possibly rename unary to edge factor

template<MessageSendingType MESSAGE_SENDING>
struct FMC_MULTICUT {
   constexpr static const char* name = "Multicut with cycle constraints";
   //constexpr static MessageSendingType MESSAGE_SENDING = MessageSendingType::SRMP;

   typedef FactorContainer<MulticutUnaryFactor, FMC_MULTICUT, 0, true> MulticutUnaryFactorContainer;
   typedef FactorContainer<MulticutTripletFactor, FMC_MULTICUT, 1> MulticutTripletFactorContainer;
   using ConstantFactorContainer = FactorContainer<ConstantFactor, FMC_MULTICUT, 2>;

   typedef MessageContainer<MulticutUnaryTripletMessage<MESSAGE_SENDING>, 0, 1, variableMessageNumber, 3, FMC_MULTICUT, 0 > MulticutUnaryTripletMessageContainer;

   using FactorList = meta::list< MulticutUnaryFactorContainer, MulticutTripletFactorContainer, ConstantFactorContainer>;
   using MessageList = meta::list<MulticutUnaryTripletMessageContainer>;

   using multicut = MulticutConstructor<FMC_MULTICUT,0,1,0,2>;
   using ProblemDecompositionList = meta::list<multicut>;
};

// It would be nice to be able to derive from FMC_MULTICUT. This is not possible due to deviating FMCs. Possibly parametrize above FMC with template
template<MessageSendingType MESSAGE_SENDING>
struct FMC_ODD_WHEEL_MULTICUT {
   constexpr static const char* name = "Multicut with cycle and odd wheel constraints";

   typedef FactorContainer<MulticutUnaryFactor, FMC_ODD_WHEEL_MULTICUT, 0, true> MulticutUnaryFactorContainer;
   typedef FactorContainer<MulticutTripletFactor, FMC_ODD_WHEEL_MULTICUT, 1> MulticutTripletFactorContainer;
   typedef FactorContainer<MulticutTripletPlusSpokeFactor, FMC_ODD_WHEEL_MULTICUT, 2> MulticutTripletPlusSpokeFactorContainer;
   using ConstantFactorContainer = FactorContainer<ConstantFactor, FMC_ODD_WHEEL_MULTICUT, 3>;
      
   typedef MessageContainer<MulticutUnaryTripletMessage<MESSAGE_SENDING>, 0, 1, variableMessageNumber, 3, FMC_ODD_WHEEL_MULTICUT, 0 > MulticutUnaryTripletMessageContainer;
   typedef MessageContainer<MulticutTripletPlusSpokeMessage, 1, 2, variableMessageNumber, variableMessageNumber, FMC_ODD_WHEEL_MULTICUT, 1> MulticutTripletPlusSpokeMessageContainer;
   typedef MessageContainer<MulticutTripletPlusSpokeCoverMessage, 1, 2, variableMessageNumber, 1, FMC_ODD_WHEEL_MULTICUT, 2> MulticutTripletPlusSpokeCoverMessageContainer;

   using FactorList = meta::list< 
      MulticutUnaryFactorContainer,
      MulticutTripletFactorContainer,
      MulticutTripletPlusSpokeFactorContainer,
      ConstantFactorContainer 
         >;
   using MessageList = meta::list<
      MulticutUnaryTripletMessageContainer, 
      MulticutTripletPlusSpokeMessageContainer, // one unfortunately cannot use a fixed size container on the right, as one spoke factor might participate in more than one odd wheel constraint
      MulticutTripletPlusSpokeCoverMessageContainer
      >;

   using multicut_c = MulticutConstructor<FMC_ODD_WHEEL_MULTICUT,0,1,0,3>;
   using multicut_cow = MulticutOddWheelConstructor<multicut_c,2,1,2>;
   using ProblemDecompositionList = meta::list<multicut_cow>;
};

struct FMC_LIFTED_MULTICUT {
   constexpr static const char* name = "Lifted Multicut with cycle constraints";
   constexpr static MessageSendingType MESSAGE_SENDING = MessageSendingType::SRMP;

   // no rounding performed: do it via GAEC and K&L, called from problem constructor
   typedef FactorContainer<MulticutUnaryFactor, FMC_LIFTED_MULTICUT, 0> MulticutUnaryFactorContainer;
   typedef FactorContainer<MulticutTripletFactor, FMC_LIFTED_MULTICUT, 1> MulticutTripletFactorContainer;
   typedef FactorContainer<LiftedMulticutCutFactor, FMC_LIFTED_MULTICUT, 2> LiftedMulticutCutFactorContainer;
   using ConstantFactorContainer = FactorContainer<ConstantFactor, FMC_LIFTED_MULTICUT, 3>;

   typedef MessageContainer<MulticutUnaryTripletMessage<MESSAGE_SENDING>, 0, 1, variableMessageNumber, 3, FMC_LIFTED_MULTICUT, 0 > MulticutUnaryTripletMessageContainer;
   typedef MessageContainer<CutEdgeLiftedMulticutFactorMessage, 0, 2, variableMessageNumber, variableMessageNumber, FMC_LIFTED_MULTICUT, 1 > CutEdgeLiftedMulticutFactorMessageContainer;
   typedef MessageContainer<LiftedEdgeLiftedMulticutFactorMessage, 0, 2, variableMessageNumber, variableMessageNumber, FMC_LIFTED_MULTICUT, 2 > LiftedEdgeLiftedMulticutFactorMessageContainer;

   using FactorList = meta::list<
      MulticutUnaryFactorContainer, 
      MulticutTripletFactorContainer, 
      LiftedMulticutCutFactorContainer,
      ConstantFactorContainer 
         >;
   using MessageList = meta::list<
      MulticutUnaryTripletMessageContainer, 
      CutEdgeLiftedMulticutFactorMessageContainer, 
      LiftedEdgeLiftedMulticutFactorMessageContainer
         >;

   using BaseMulticutConstructor = MulticutConstructor<FMC_LIFTED_MULTICUT,0,1,0,2>;
   using LiftedMulticutConstructor = LiftedMulticutConstructor<BaseMulticutConstructor,2,1,2>;
   using ProblemDecompositionList = meta::list<LiftedMulticutConstructor>;


//   using MessageListHana = hana::tuple<
//      MulticutUnaryTripletMessageContainer, 
//      CutEdgeLiftedMulticutFactorMessageContainer, 
//      LiftedEdgeLiftedMulticutFactorMessageContainer
//         >;

//   using FactorListHana = hana::tuple<
//      MulticutUnaryFactorContainer, 
//      MulticutTripletFactorContainer, 
//      LiftedMulticutCutFactorContainer 
//         >;
   //using MessageListHana = hana::tuple<
   //   MulticutUnaryTripletMessageContainer, 
   //   CutEdgeLiftedMulticutFactorMessageContainer, 
   //   LiftedEdgeLiftedMulticutFactorMessageContainer
   //      >;
//   using ProblemDecompositionListHana = hana::tuple<LiftedMulticutConstructor>;

   /*
   hana::tuple_t<
      MulticutUnaryFactorContainer, 
      MulticutTripletFactorContainer, 
      LiftedMulticutCutFactorContainer 
         > FactorListHana;
   hana::tuple_t<
      MulticutUnaryTripletMessageContainer, 
      CutEdgeLiftedMulticutFactorMessageContainer, 
      LiftedEdgeLiftedMulticutFactorMessageContainer
         > MessageListHana;
   hana::tuple_t<LiftedMulticutConstructor> ProblemDecompositionListHana;
   */

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

   template<typename FMC>
   bool ParseProblem(const std::string filename, Solver<FMC>& pd)
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

      for(INDEX i=0; i<edges.size(); ++i) {
         INDEX i1 = std::get<0>(edges[i]); 
         INDEX i2 = std::get<1>(edges[i]); 
         if(i1 > i2) {
            std::swap(i1,i2);
         }
         const REAL cost = std::get<2>(edges[i]);
         //std::cout << "add edge (" << i1 << "," << i2 << ") with cost " << cost << "\n";
         mc.AddUnaryFactor(i1, i2, cost);
      }

      return true;
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

   template<typename FMC, typename Rule >
      struct action
      : pegtl::nothing< Rule > {};


   template<typename FMC> struct action<FMC, positive_integer > {
      template<typename INPUT>
      static void apply(const INPUT & in, Solver<FMC>&, std::stack<SIGNED_INDEX>& integer_stack, std::stack<REAL>& real_stack, MulticutInput &)
      {
         integer_stack.push(std::stoul(in.string())); 
      }
   };
   template<typename FMC> struct action<FMC, real_number > {
      template<typename INPUT>
      static void apply(const INPUT & in, Solver<FMC>&, std::stack<SIGNED_INDEX>& integer_stack, std::stack<REAL>& real_stack, MulticutInput &)
      {
         real_stack.push(std::stod(in.string())); 
      }
   };
   template<typename FMC> struct action<FMC, numberOfVariables_line > {
      template<typename INPUT>
      static void apply(const INPUT & in, Solver<FMC>&, std::stack<SIGNED_INDEX>& integer_stack, std::stack<REAL>& real_stack, MulticutInput & mcInput)
      {
         assert(integer_stack.size() == 1);
         mcInput.numberOfVariables_ = integer_stack.top();
         integer_stack.pop();
      }
   };
   template<typename FMC> struct action<FMC, edge_line > {
      template<typename INPUT>
      static void apply(const INPUT & in, Solver<FMC>& pd, std::stack<SIGNED_INDEX>& integer_stack, std::stack<REAL>& real_stack, MulticutInput & mcInput)
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
   template<typename FMC> struct action<FMC, pegtl::eof> {
      template<typename INPUT>
      static void apply(const INPUT & in, Solver<FMC>& pd, std::stack<SIGNED_INDEX>& integer_stack, std::stack<REAL>& real_stack, MulticutInput& mcInput)
      {}
   };

   template<typename FMC>
      struct actionSpecialization {
         template<typename RULE> struct type : public action<FMC,RULE> {};
      };


   template<typename FMC> struct action<FMC, lifted_edge_line > {
      template<typename INPUT>
      static void apply(const INPUT & in, Solver<FMC>& pd, std::stack<SIGNED_INDEX>& integer_stack, std::stack<REAL>& real_stack, MulticutInput & mcInput)
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

      
   template<typename FMC>
   bool ParseProblem(const std::string filename, Solver<FMC>& pd)
   {
      std::stack<SIGNED_INDEX> integer_stack;
      std::stack<REAL> real_stack;
      MulticutInput mcInput;
      std::cout << "parsing " << filename << "\n";

      pegtl::file_parser problem(filename);

      return problem.parse< grammar, actionSpecialization<FMC>::template type >(pd, integer_stack, real_stack, mcInput);
   }

   template<typename FMC>
   bool ParseLiftedProblem(const std::string filename, Solver<FMC>& pd)
   {
      std::stack<SIGNED_INDEX> integer_stack;
      std::stack<REAL> real_stack;
      MulticutInput mcInput;
      std::cout << "parsing " << filename << "\n";

      pegtl::file_parser problem(filename);
      return problem.parse< LiftedMulticutGrammar, actionSpecialization<FMC>::template type >(pd, integer_stack, real_stack, mcInput);
   }


} // end namespace MulticutTextInput

// HDF5 input as in data of Andres et al.
namespace MulticutH5Input {

   template<typename SOLVER>
   bool ParseLiftedProblem(const std::string filename, SOLVER& pd)
   {
      auto& mc = pd.template GetProblemConstructor<0>();
      
      andres::graph::GridGraph<2> originalGraph;
      andres::graph::Graph<> liftedGraph;
      std::vector<REAL> edgeValues;

      auto fileHandle = andres::graph::hdf5::openFile(filename);
      andres::graph::hdf5::load(fileHandle, "graph", originalGraph);
      andres::graph::hdf5::load(fileHandle, "graph-lifted", liftedGraph);

      std::vector<size_t> shape;
      andres::graph::hdf5::load(fileHandle, "edge-cut-probabilities", shape, edgeValues);
      andres::graph::hdf5::closeFile(fileHandle);

      // transform to energy cost
      {
         std::transform(edgeValues.begin(), edgeValues.end(), edgeValues.begin(), andres::NegativeLogProbabilityRatio<REAL,REAL>());
      }

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

#endif // LP_MP_MULTICUT_HXX
