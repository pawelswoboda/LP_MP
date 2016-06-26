#ifndef LP_MP_MULTICUT_HXX
#define LP_MP_MULTICUT_HXX

#include "LP_MP.h"
#include "problem_decomposition.hxx"
#include "factors_messages.hxx"
#include "multicut_unary_factor.hxx"
#include "multicut_triplet_factor.hxx"
#include "multicut_global_factor.hxx"
#include "multicut_odd_wheel_factor.hxx"
#include "multicut_odd_wheel.hxx"
#include "lifted_multicut_factors_messages.hxx"
#include "multicut_constructor.hxx"

#include "parse_rules.h"

#include "hdf5_routines.hxx"

#include <iostream>
#include <vector>

namespace LP_MP {

// do zrobienia: possibly rename unary to edge factor

//template<MessageSendingType MESSAGE_SENDING>
struct FMC_MULTICUT {
   constexpr static const char* name = "Multicut with odd cycle constraints";
   constexpr static MessageSendingType MESSAGE_SENDING = MessageSendingType::SRMP;

   typedef FactorContainer<MulticutUnaryFactor, FixedSizeExplicitRepamStorage<MulticutUnaryFactor::size()>::type, FMC_MULTICUT, 0, true> MulticutUnaryFactorContainer;
   typedef FactorContainer<MulticutTripletFactor, FixedSizeExplicitRepamStorage<MulticutTripletFactor::size()>::type, FMC_MULTICUT, 1> MulticutTripletFactorContainer;
   typedef FactorContainer<MulticutGlobalFactor, MulticutGlobalRepamStorage, FMC_MULTICUT, 2> MulticutGlobalFactorContainer;

   typedef MessageContainer<MulticutUnaryTripletMessage<MESSAGE_SENDING>, 0, 1, variableMessageNumber, 3, MulticutUnaryTripletMessage<MESSAGE_SENDING>::size(), FMC_MULTICUT, 0 > MulticutUnaryTripletMessageContainer;
   typedef MessageContainer<MulticutUnaryGlobalMessage, 0, 2, 1, variableMessageNumber, 0, FMC_MULTICUT, 1> MulticutUnaryGlobalMessageContainer;

   using FactorList = meta::list< MulticutUnaryFactorContainer, MulticutTripletFactorContainer, MulticutGlobalFactorContainer >;
   using MessageList = meta::list<MulticutUnaryTripletMessageContainer, MulticutUnaryGlobalMessageContainer>;

   using multicut = MulticutConstructor<FMC_MULTICUT,0,1,2,0,1>;
   using ProblemDecompositionList = meta::list<multicut>;
};

// It would be nice to be able to derive from FMC_MULTICUT. This is not possible due to deviating FMCs
template<MessageSendingType MESSAGE_SENDING>
struct FMC_ODD_WHEEL_MULTICUT {
   constexpr static const char* name = "Multicut with odd cycle and wheel constraints";

   typedef FactorContainer<MulticutUnaryFactor, FixedSizeExplicitRepamStorage<MulticutUnaryFactor::size()>::type, FMC_ODD_WHEEL_MULTICUT, 0, true> MulticutUnaryFactorContainer;
   typedef FactorContainer<MulticutTripletFactor, FixedSizeExplicitRepamStorage<MulticutTripletFactor::size()>::type, FMC_ODD_WHEEL_MULTICUT, 1> MulticutTripletFactorContainer;
   typedef FactorContainer<MulticutGlobalFactor, MulticutGlobalRepamStorage, FMC_ODD_WHEEL_MULTICUT, 2> MulticutGlobalFactorContainer;

   typedef FactorContainer<MulticutTripletPlusSpokeFactor, FixedSizeExplicitRepamStorage<MulticutTripletPlusSpokeFactor::size()>::type, FMC_ODD_WHEEL_MULTICUT, 3> MulticutTripletPlusSpokeFactorContainer;
      
   typedef MessageContainer<MulticutUnaryTripletMessage<MESSAGE_SENDING>, 0, 1, variableMessageNumber, 3, MulticutUnaryTripletMessage<MESSAGE_SENDING>::size(), FMC_ODD_WHEEL_MULTICUT, 0 > MulticutUnaryTripletMessageContainer;
   typedef MessageContainer<MulticutUnaryGlobalMessage, 0, 2, 1, variableMessageNumber, 0, FMC_ODD_WHEEL_MULTICUT, 1> MulticutUnaryGlobalMessageContainer;
   
   typedef MessageContainer<MulticutTripletPlusSpokeMessage, 1, 3, variableMessageNumber, variableMessageNumber, MulticutTripletPlusSpokeMessage::size(), FMC_ODD_WHEEL_MULTICUT, 2> MulticutTripletPlusSpokeMessageContainer;
   typedef MessageContainer<MulticutTripletPlusSpokeCoverMessage, 1, 3, variableMessageNumber, 1, MulticutTripletPlusSpokeCoverMessage::size(), FMC_ODD_WHEEL_MULTICUT, 3> MulticutTripletPlusSpokeCoverMessageContainer;

   using FactorList = meta::list< MulticutUnaryFactorContainer, MulticutTripletFactorContainer, MulticutGlobalFactorContainer, MulticutTripletPlusSpokeFactorContainer >;
   using MessageList = meta::list<
      MulticutUnaryTripletMessageContainer, 
      MulticutUnaryGlobalMessageContainer, 
      MulticutTripletPlusSpokeMessageContainer, // one unfortunately cannot use a fixed size container on the right, as one spoke factor might participate in more than one odd wheel constraint
      MulticutTripletPlusSpokeCoverMessageContainer
      >;

   using multicut = MulticutOddWheelConstructor<FMC_ODD_WHEEL_MULTICUT,0,1,2,0,1,3,2,3>;
   using ProblemDecompositionList = meta::list<multicut>;
};

struct FMC_LIFTED_MULTICUT_ODD_CYCLE {
   constexpr static const char* name = "Lifted Multicut";
   constexpr static MessageSendingType MESSAGE_SENDING = MessageSendingType::SRMP;

   typedef FactorContainer<MulticutUnaryFactor, FixedSizeExplicitRepamStorage<MulticutUnaryFactor::size()>::type, FMC_LIFTED_MULTICUT_ODD_CYCLE, 0, true> MulticutUnaryFactorContainer;
   typedef FactorContainer<MulticutTripletFactor, FixedSizeExplicitRepamStorage<MulticutTripletFactor::size()>::type, FMC_LIFTED_MULTICUT_ODD_CYCLE, 1> MulticutTripletFactorContainer;
   typedef FactorContainer<MulticutGlobalFactor, MulticutGlobalRepamStorage, FMC_LIFTED_MULTICUT_ODD_CYCLE, 2> MulticutGlobalFactorContainer;
   
   typedef MessageContainer<MulticutUnaryTripletMessage<MESSAGE_SENDING>, 0, 1, variableMessageNumber, 3, MulticutUnaryTripletMessage<MESSAGE_SENDING>::size(), FMC_LIFTED_MULTICUT_ODD_CYCLE, 0 > MulticutUnaryTripletMessageContainer;
   typedef MessageContainer<MulticutUnaryGlobalMessage, 0, 2, 1, variableMessageNumber, 0, FMC_LIFTED_MULTICUT_ODD_CYCLE, 1> MulticutUnaryGlobalMessageContainer;

   typedef FactorContainer<LiftedMulticutCutFactor, ExplicitRepamStorage, FMC_LIFTED_MULTICUT_ODD_CYCLE, 3> LiftedMulticutCutFactorContainer;
   typedef MessageContainer<CutEdgeLiftedMulticutFactorMessage, 0, 3, variableMessageNumber, variableMessageNumber, CutEdgeLiftedMulticutFactorMessage::size(), FMC_LIFTED_MULTICUT_ODD_CYCLE, 2 > CutEdgeLiftedMulticutFactorMessageContainer;
   typedef MessageContainer<LiftedEdgeLiftedMulticutFactorMessage, 0, 3, variableMessageNumber, variableMessageNumber, LiftedEdgeLiftedMulticutFactorMessage::size(), FMC_LIFTED_MULTICUT_ODD_CYCLE, 3 > LiftedEdgeLiftedMulticutFactorMessageContainer;

   using FactorList = meta::list<MulticutUnaryFactorContainer, MulticutTripletFactorContainer, MulticutGlobalFactorContainer, LiftedMulticutCutFactorContainer >;
   using MessageList = meta::list<MulticutUnaryTripletMessageContainer, MulticutUnaryGlobalMessageContainer, CutEdgeLiftedMulticutFactorMessageContainer, LiftedEdgeLiftedMulticutFactorMessageContainer>;

   using BaseMulticutConstructor = MulticutConstructor<FMC_LIFTED_MULTICUT_ODD_CYCLE,0,1,2,0,1>;
   using LiftedMulticutConstructor = LiftedMulticutConstructor<BaseMulticutConstructor,3,2,3>;
   using ProblemDecompositionList = meta::list<LiftedMulticutConstructor>;
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
   bool ParseProblem(const std::string filename, ProblemDecomposition<FMC>& pd)
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
      for(INDEX i=0; i<edges.size(); ++i) {
         std::get<2>(edges[i]) = -edgeCosts[2*i] + edgeCosts[2*i+1]; // do zrobienia: or the other way around
         //std::get<2>(edges[i]) = edgeCosts[2*i] - edgeCosts[2*i+1];
      }
      // note: theoretically, this could be much more complicated, if some factors are shared etc. Then this simple approach above will not work.

      auto& mc = pd.template GetProblemConstructor<0>();
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
      static void apply(const pegtl::input & in, ProblemDecomposition<FMC>&, std::stack<SIGNED_INDEX>& integer_stack, std::stack<REAL>& real_stack, MulticutInput &)
      {
         integer_stack.push(std::stoul(in.string())); 
      }
   };
   template<typename FMC> struct action<FMC, real_number > {
      static void apply(const pegtl::input & in, ProblemDecomposition<FMC>&, std::stack<SIGNED_INDEX>& integer_stack, std::stack<REAL>& real_stack, MulticutInput &)
      {
         real_stack.push(std::stod(in.string())); 
      }
   };
   template<typename FMC> struct action<FMC, numberOfVariables_line > {
      static void apply(const pegtl::input & in, ProblemDecomposition<FMC>&, std::stack<SIGNED_INDEX>& integer_stack, std::stack<REAL>& real_stack, MulticutInput & mcInput)
      {
         assert(integer_stack.size() == 1);
         mcInput.numberOfVariables_ = integer_stack.top();
         integer_stack.pop();
      }
   };
   template<typename FMC> struct action<FMC, edge_line > {
      static void apply(const pegtl::input & in, ProblemDecomposition<FMC>& pd, std::stack<SIGNED_INDEX>& integer_stack, std::stack<REAL>& real_stack, MulticutInput & mcInput)
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
      static void apply(const pegtl::input & in, ProblemDecomposition<FMC>& pd, std::stack<SIGNED_INDEX>& integer_stack, std::stack<REAL>& real_stack, MulticutInput& mcInput)
      {}
   };

   template<typename FMC>
      struct actionSpecialization {
         template<typename RULE> struct type : public action<FMC,RULE> {};
      };


   template<typename FMC> struct action<FMC, lifted_edge_line > {
      static void apply(const pegtl::input & in, ProblemDecomposition<FMC>& pd, std::stack<SIGNED_INDEX>& integer_stack, std::stack<REAL>& real_stack, MulticutInput & mcInput)
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
   bool ParseProblem(const std::string filename, ProblemDecomposition<FMC>& pd)
   {
      std::stack<SIGNED_INDEX> integer_stack;
      std::stack<REAL> real_stack;
      MulticutInput mcInput;
      std::cout << "parsing " << filename << "\n";

      pegtl::file_parser problem(filename);

      return problem.parse< grammar, actionSpecialization<FMC>::template type >(pd, integer_stack, real_stack, mcInput);
   }

   template<typename FMC>
   bool ParseLiftedProblem(const std::string filename, ProblemDecomposition<FMC>& pd)
   {
      std::stack<SIGNED_INDEX> integer_stack;
      std::stack<REAL> real_stack;
      MulticutInput mcInput;
      std::cout << "parsing " << filename << "\n";

      pegtl::file_parser problem(filename);
      return problem.parse< LiftedMulticutGrammar, actionSpecialization<FMC>::template type >(pd, integer_stack, real_stack, mcInput);
   }


} // end namespace MulticutTextInput

} // end namespace LP_MP

#endif // LP_MP_MULTICUT_HXX
