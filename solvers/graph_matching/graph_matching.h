#ifndef LP_MP_GRAPH_MATCHING_H
#define LP_MP_GRAPH_MATCHING_H

//#include "problem_decomposition.hxx"
#include "factors_messages.hxx"
#include "LP_MP.h"
#include "solver.hxx"
#include "factors/simplex_factor.hxx"
#include "const_array_types.h"
#include "messages/simplex_marginalization_message.hxx"
#include "messages/equality_message.hxx"
#include "problem_constructors/mrf_problem_construction.hxx"
//#include "factors/min_cost_flow_factor_lemon.hxx"
#include "factors/min_cost_flow_factor_cs2.hxx"
//#include "factors/minimum_cost_flow_labeling.hxx"
#include "../cosegmentation/assignment_via_min_cost_flow_constructor.hxx" // move file to problem_constructors
#include "../cosegmentation/assignment_via_message_passing_problem_constructor.hxx" // move file to problem_constructors

#include "problem_constructors/cycle_inequalities.hxx"

#include "parse_rules.h"

#include <vector>

// this file contains definitions of various graph matching solvers and grammars.
//
// solvers:
// FMC_MP implements graph matching with the uniqueness constraints implemented via messages.
// FMC_MCF implements graph matching with a global min cost flow factor.
// FMC_GM amounts to TRWS with infinity on diagonals
// iFMC_${MODEL}_T implements tightening version of all three solvers using violated cycle tightening of Sontag
//
// input grammars:
// TorresaniEtAlInput contains the grammar used by the dual decomposition algorithm of Torresani, Kolmogorov and Rother.
// UAIInput contains the grammar in uai MRF format plus constraints section.

using namespace LP_MP;

enum class PairwiseConstruction {Left,Right,BothSides}; // Indicates whether pairwise potentials should be built on {left|right|both} side(s) of assignment graph.

typedef UnaryLoop<> UnaryLoopType;
typedef PairwiseLoop<0> LeftLoopType;
typedef PairwiseLoop<1> RightLoopType;

typedef SimplexMarginalizationMessage<UnaryLoopType,LeftLoopType,true,false,false,true> LeftMargMessage;
typedef SimplexMarginalizationMessage<UnaryLoopType,RightLoopType,true,false,false,true> RightMargMessage;

typedef PairwiseTripletLoop<0,1> PairwiseTripletLoopType12;
typedef PairwiseTripletLoop<0,2> PairwiseTripletLoopType13;
typedef PairwiseTripletLoop<1,2> PairwiseTripletLoopType23;
typedef SimplexMarginalizationMessage<UnaryLoopType,PairwiseTripletLoopType12,true,false,false,true> PairwiseTriplet12Message;
typedef SimplexMarginalizationMessage<UnaryLoopType,PairwiseTripletLoopType13,true,false,false,true> PairwiseTriplet13Message;
typedef SimplexMarginalizationMessage<UnaryLoopType,PairwiseTripletLoopType23,true,false,false,true> PairwiseTriplet23Message;

typedef SimplexFactor<> Simplex;

// graph matching with assignment via message passing
template<PairwiseConstruction PAIRWISE_CONSTRUCTION = PairwiseConstruction::Left>
struct FMC_MP {
   using FMC_MP_PARAM = FMC_MP<PAIRWISE_CONSTRUCTION>;
   constexpr static const char* name =
      PAIRWISE_CONSTRUCTION == PairwiseConstruction::Left ? "AMP-O"
      : (PAIRWISE_CONSTRUCTION == PairwiseConstruction::Right ? "AMP-I"
      : (PAIRWISE_CONSTRUCTION == PairwiseConstruction::BothSides ? "AMP-B"
      : "unknown variant"));
      
   typedef FactorContainer<Simplex, ExplicitRepamStorage, FMC_MP_PARAM, 0, true, true > UnaryFactor; // set to true if labeling by unaries is desired
   typedef FactorContainer<Simplex, ExplicitRepamStorage, FMC_MP_PARAM, 1, false, false > PairwiseFactor;
   //typedef FactorContainer<MinimumCostFlowLabelingFactor, MinimumCostFlowLabelingRepamStorage, FMC_MP_PARAM, 2, true, false> McfLabelingFactor;

   typedef MessageContainer<EqualityMessage, 0, 0, variableMessageNumber, variableMessageNumber, 1, FMC_MP_PARAM, 0 > AssignmentConstraintMessage;
   typedef MessageContainer<LeftMargMessage, 0, 1, variableMessageNumber, 1, variableMessageSize, FMC_MP_PARAM, 1 > UnaryPairwiseMessageLeft;
   typedef MessageContainer<RightMargMessage, 0, 1, variableMessageNumber, 1, variableMessageSize, FMC_MP_PARAM, 2 > UnaryPairwiseMessageRight;
   //typedef MessageContainer<SimplexMinimumCostFlowLabelingMessage, 0, 2, 1, variableMessageNumber, 0, FMC_MP_PARAM, 3 > UnaryMcfLabelingMessage;

   using FactorList = meta::list< UnaryFactor, PairwiseFactor>;//, McfLabelingFactor >;
   using MessageList = meta::list< AssignmentConstraintMessage, UnaryPairwiseMessageLeft, UnaryPairwiseMessageRight>;//, UnaryMcfLabelingMessage >;

   using assignment = AssignmentViaMessagePassingProblemConstructor<FMC_MP_PARAM,0,0>;
   using mrfLeft = StandardMrfConstructor<FMC_MP_PARAM,0,1,1,2>;
   using mrfRight = StandardMrfConstructor<FMC_MP_PARAM,0,1,1,2>;
   //using mcfLabeling = MinimumCostFlowLabelingConstructor<FMC_MP_PARAM,0,2,3>;
   using ProblemDecompositionList = meta::list<assignment, mrfLeft, mrfRight>;//, mcfLabeling>;
};

// graph matching with assignment via message passing + tightening triplets
template<PairwiseConstruction PAIRWISE_CONSTRUCTION = PairwiseConstruction::Left>
struct FMC_MP_T {
   using FMC_MP_PARAM = FMC_MP_T<PAIRWISE_CONSTRUCTION>;
   constexpr static const char* name =
      PAIRWISE_CONSTRUCTION == PairwiseConstruction::Left ? "AMP-O-T"
      : (PAIRWISE_CONSTRUCTION == PairwiseConstruction::Right ? "AMP-I-T"
      : (PAIRWISE_CONSTRUCTION == PairwiseConstruction::BothSides ? "AMP-B-T"
      : "unknown variant"));
      
   typedef FactorContainer<Simplex, ExplicitRepamStorage, FMC_MP_PARAM, 0, false, true > UnaryFactor; // set to true if labeling by unaries is desired
   typedef FactorContainer<Simplex, ExplicitRepamStorage, FMC_MP_PARAM, 1, false, false > PairwiseFactor;
   //typedef FactorContainer<MinimumCostFlowLabelingFactor, MinimumCostFlowLabelingRepamStorage, FMC_MP_PARAM, 2, true, false> McfLabelingFactor;

   typedef MessageContainer<EqualityMessage, 0, 0, variableMessageNumber, variableMessageNumber, 1, FMC_MP_PARAM, 0 > AssignmentConstraintMessage;
   typedef MessageContainer<LeftMargMessage, 0, 1, variableMessageNumber, 1, variableMessageSize, FMC_MP_PARAM, 1 > UnaryPairwiseMessageLeft;
   typedef MessageContainer<RightMargMessage, 0, 1, variableMessageNumber, 1, variableMessageSize, FMC_MP_PARAM, 2 > UnaryPairwiseMessageRight;
   //typedef MessageContainer<SimplexMinimumCostFlowLabelingMessage, 0, 2, 1, variableMessageNumber, 0, FMC_MP_PARAM, 3 > UnaryMcfLabelingMessage;

   typedef FactorContainer<Simplex, ExplicitRepamStorage, FMC_MP_PARAM, 2 > EmptyTripletFactor;
   typedef MessageContainer<PairwiseTriplet12Message, 1, 2, variableMessageNumber, 1, variableMessageSize, FMC_MP_PARAM, 4> PairwiseTriplet12MessageContainer;
   typedef MessageContainer<PairwiseTriplet13Message, 1, 2, variableMessageNumber, 1, variableMessageSize, FMC_MP_PARAM, 5> PairwiseTriplet13MessageContainer;
   typedef MessageContainer<PairwiseTriplet23Message, 1, 2, variableMessageNumber, 1, variableMessageSize, FMC_MP_PARAM, 6> PairwiseTriplet23MessageContainer;

   using FactorList = meta::list< UnaryFactor, PairwiseFactor, EmptyTripletFactor>;//McfLabelingFactor, EmptyTripletFactor >;
   using MessageList = meta::list< 
      AssignmentConstraintMessage,
      UnaryPairwiseMessageLeft,
      UnaryPairwiseMessageRight,
      //UnaryMcfLabelingMessage,
      PairwiseTriplet12MessageContainer, 
      PairwiseTriplet13MessageContainer, 
      PairwiseTriplet23MessageContainer 
         >;

   using assignment = AssignmentViaMessagePassingProblemConstructor<FMC_MP_PARAM,0,0>;
   using mrf = StandardMrfConstructor<FMC_MP_PARAM,0,1,1,2>;
   using tighteningMrf = TighteningMRFProblemConstructor<mrf,3,4,5,6>;
   using mrfLeft = tighteningMrf;
   using mrfRight = tighteningMrf;
   //using mcfLabeling = MinimumCostFlowLabelingConstructor<FMC_MP_PARAM,0,2,3>;
   using ProblemDecompositionList = meta::list<assignment, mrfLeft, mrfRight>;//, mcfLabeling>;
};



// graph matching with assignment via minimum cost flow solver

// first good option: construct pairwise potentials on both sides, only send messages from unary to assignment (no receiving) and adjust all factors.
// another good algorithm can be obtained by just using the left side, only receiving messages from unaries (no sending) and adjusting all factors.
// For both algorithms: In maximal perturbation problem capacities have to be set [-1,-inf] for edges = 1 and [0,inf] for edges = 0

template<PairwiseConstruction PAIRWISE_CONSTRUCTION = PairwiseConstruction::Left>
struct FMC_MCF {
   using FMC_MCF_PARAM = FMC_MCF<PAIRWISE_CONSTRUCTION>;
   constexpr static const char* name =
      PAIRWISE_CONSTRUCTION == PairwiseConstruction::Left ? "AMCF-O"
      : (PAIRWISE_CONSTRUCTION == PairwiseConstruction::Right ? "AMCF-I"
      : (PAIRWISE_CONSTRUCTION == PairwiseConstruction::BothSides ? "AMCF-B"
      : "unknown variant"));
      
   constexpr static INDEX McfCoveringFactor = PAIRWISE_CONSTRUCTION == PairwiseConstruction::BothSides ? 2 : 1;

   typedef FactorContainer<MinCostFlowFactorCS2, MinCostFlowReparametrizationStorageCS2, FMC_MCF_PARAM, 0, true, true> MinCostFlowAssignmentFactor;
   typedef FactorContainer<Simplex, ExplicitRepamStorage, FMC_MCF_PARAM, 1 > UnaryFactor;
   typedef FactorContainer<Simplex, ExplicitRepamStorage, FMC_MCF_PARAM, 2 > PairwiseFactor;

   typedef MessageContainer<LeftMargMessage, 1, 2, variableMessageNumber, 1, variableMessageSize, FMC_MCF_PARAM, 0 > UnaryPairwiseMessageLeft;
   typedef MessageContainer<RightMargMessage, 1, 2, variableMessageNumber, 1, variableMessageSize, FMC_MCF_PARAM, 1 > UnaryPairwiseMessageRight;
   typedef MessageContainer<UnaryToAssignmentMessage<McfCoveringFactor>, 1, 0, 1, variableMessageNumber, variableMessageSize, FMC_MCF_PARAM, 2> UnaryToAssignmentMessageContainer;

   using FactorList = meta::list<MinCostFlowAssignmentFactor, UnaryFactor, PairwiseFactor>;
   using MessageList = meta::list<
       UnaryPairwiseMessageLeft,  
       UnaryPairwiseMessageRight, 
       UnaryToAssignmentMessageContainer
      >;

   using assignment = AssignmentViaMinCostFlowConstructor<FMC_MCF_PARAM,0>;
   using mrf = StandardMrfConstructor<FMC_MCF_PARAM,1,2,0,1>;
   using mrfLeft = mrf;
   using mrfRight = mrf;
   using ProblemDecompositionList = meta::list<assignment,mrfLeft,mrfRight>;
};

// + tightening
template<PairwiseConstruction PAIRWISE_CONSTRUCTION = PairwiseConstruction::Left>
struct FMC_MCF_T {
   using FMC_MCF_PARAM = FMC_MCF_T<PAIRWISE_CONSTRUCTION>;
   constexpr static const char* name =
      PAIRWISE_CONSTRUCTION == PairwiseConstruction::Left ? "AMCF-O-T"
      : (PAIRWISE_CONSTRUCTION == PairwiseConstruction::Right ? "AMCF-I-T"
      : (PAIRWISE_CONSTRUCTION == PairwiseConstruction::BothSides ? "AMCF-B-T"
      : "unknown variant"));
      
   constexpr static INDEX McfCoveringFactor = PAIRWISE_CONSTRUCTION == PairwiseConstruction::BothSides ? 2 : 1;

   typedef FactorContainer<MinCostFlowFactorCS2, MinCostFlowReparametrizationStorageCS2, FMC_MCF_PARAM, 0, true, true> MinCostFlowAssignmentFactor;
   typedef FactorContainer<Simplex, ExplicitRepamStorage, FMC_MCF_PARAM, 1 > UnaryFactor;
   typedef FactorContainer<Simplex, ExplicitRepamStorage, FMC_MCF_PARAM, 2 > PairwiseFactor;
   typedef FactorContainer<Simplex, ExplicitRepamStorage, FMC_MCF_PARAM, 3 > EmptyTripletFactor;

   typedef MessageContainer<LeftMargMessage, 1, 2, variableMessageNumber, 1, variableMessageSize, FMC_MCF_PARAM, 0 > UnaryPairwiseMessageLeft;
   typedef MessageContainer<RightMargMessage, 1, 2, variableMessageNumber, 1, variableMessageSize, FMC_MCF_PARAM, 1 > UnaryPairwiseMessageRight;
   typedef MessageContainer<UnaryToAssignmentMessage<McfCoveringFactor>, 1, 0, 1, variableMessageNumber, variableMessageSize, FMC_MCF_PARAM, 2> UnaryToAssignmentMessageContainer;

   typedef MessageContainer<PairwiseTriplet12Message, 2, 3, variableMessageNumber, 1, variableMessageSize, FMC_MCF_PARAM, 3> PairwiseTriplet12MessageContainer;
   typedef MessageContainer<PairwiseTriplet13Message, 2, 3, variableMessageNumber, 1, variableMessageSize, FMC_MCF_PARAM, 4> PairwiseTriplet13MessageContainer;
   typedef MessageContainer<PairwiseTriplet23Message, 2, 3, variableMessageNumber, 1, variableMessageSize, FMC_MCF_PARAM, 5> PairwiseTriplet23MessageContainer;

   using FactorList = meta::list<MinCostFlowAssignmentFactor, UnaryFactor, PairwiseFactor, EmptyTripletFactor>;
   using MessageList = meta::list<
       UnaryPairwiseMessageLeft,  
       UnaryPairwiseMessageRight, 
       UnaryToAssignmentMessageContainer, 
       PairwiseTriplet12MessageContainer, 
       PairwiseTriplet13MessageContainer, 
       PairwiseTriplet23MessageContainer 
      >;

   using assignment = AssignmentViaMinCostFlowConstructor<FMC_MCF_PARAM,0>;
   using mrf = StandardMrfConstructor<FMC_MCF_PARAM,1,2,0,1>;
   using tighteningMrf = TighteningMRFProblemConstructor<mrf,3,3,4,5>;
   using mrfLeft = tighteningMrf;
   using mrfRight = tighteningMrf;
   using ProblemDecompositionList = meta::list<assignment,mrfLeft,mrfRight>;
};


// naive version where the assignment is enforced through inf on pairwise diagonals. One has to insert all possible diagonals then.
// this results in a dense standard pairwise graphical model
template<PairwiseConstruction PAIRWISE_CONSTRUCTION = PairwiseConstruction::Left> // note: both sides makes no sense here
struct FMC_GM {
   using FMC_GM_PARAM = FMC_GM<PAIRWISE_CONSTRUCTION>;
   constexpr static const char* name =
      PAIRWISE_CONSTRUCTION == PairwiseConstruction::Left ? "GM-O"
      : (PAIRWISE_CONSTRUCTION == PairwiseConstruction::Right ? "GM-I"
      : "unknown variant");
      
   typedef FactorContainer<Simplex, ExplicitRepamStorage, FMC_GM_PARAM, 0, true, true > UnaryFactor; // make true, if primal rounding similar to TRW-S is required
   typedef FactorContainer<Simplex, ExplicitRepamStorage, FMC_GM_PARAM, 1, false, false > PairwiseFactor;
   //typedef FactorContainer<MinimumCostFlowLabelingFactor, MinimumCostFlowLabelingRepamStorage, FMC_GM_PARAM, 2, true, false> McfLabelingFactor;

   typedef MessageContainer<LeftMargMessage, 0, 1, variableMessageNumber, 1, variableMessageSize, FMC_GM_PARAM, 0 > UnaryPairwiseMessageLeft;
   typedef MessageContainer<RightMargMessage, 0, 1, variableMessageNumber, 1, variableMessageSize, FMC_GM_PARAM, 1 > UnaryPairwiseMessageRight;
   //typedef MessageContainer<SimplexMinimumCostFlowLabelingMessage, 0, 2, 1, variableMessageNumber, 0, FMC_GM_PARAM, 2 > UnaryMcfLabelingMessage;

   using FactorList = meta::list< UnaryFactor, PairwiseFactor>;//, McfLabelingFactor >;
   using MessageList = meta::list< UnaryPairwiseMessageLeft, UnaryPairwiseMessageRight>;//, UnaryMcfLabelingMessage >;

   using mrf = StandardMrfConstructor<FMC_GM_PARAM,0,1,0,1>;
   //using mcfLabeling = MinimumCostFlowLabelingConstructor<FMC_GM_PARAM,0,2,2>;
   using ProblemDecompositionList = meta::list<mrf>;//,mcfLabeling>;
};

// + tightening triplets
template<PairwiseConstruction PAIRWISE_CONSTRUCTION = PairwiseConstruction::Left> // note: both sides makes no sense here
struct FMC_GM_T {
   using FMC_GM_PARAM = FMC_GM_T<PAIRWISE_CONSTRUCTION>;
   constexpr static const char* name =
      PAIRWISE_CONSTRUCTION == PairwiseConstruction::Left ? "GM-O-T"
      : (PAIRWISE_CONSTRUCTION == PairwiseConstruction::Right ? "GM-I-T"
      : "unknown variant");
      
   typedef FactorContainer<Simplex, ExplicitRepamStorage, FMC_GM_PARAM, 0, true, false > UnaryFactor; // make true, if primal rounding similar to TRW-S is required
   typedef FactorContainer<Simplex, ExplicitRepamStorage, FMC_GM_PARAM, 1, false, false > PairwiseFactor;
   //typedef FactorContainer<MinimumCostFlowLabelingFactor, MinimumCostFlowLabelingRepamStorage, FMC_GM_PARAM, 2, true, false> McfLabelingFactor;

   typedef MessageContainer<LeftMargMessage, 0, 1, variableMessageNumber, 1, variableMessageSize, FMC_GM_PARAM, 0 > UnaryPairwiseMessageLeft;
   typedef MessageContainer<RightMargMessage, 0, 1, variableMessageNumber, 1, variableMessageSize, FMC_GM_PARAM, 1 > UnaryPairwiseMessageRight;
   //typedef MessageContainer<SimplexMinimumCostFlowLabelingMessage, 0, 2, 1, variableMessageNumber, 0, FMC_GM_PARAM, 2 > UnaryMcfLabelingMessage;

   // tightening
   typedef FactorContainer<Simplex, ExplicitRepamStorage, FMC_GM_PARAM, 2 > EmptyTripletFactor;
   typedef MessageContainer<PairwiseTriplet12Message, 1, 3, variableMessageNumber, 1, variableMessageSize, FMC_GM_PARAM, 2> PairwiseTriplet12MessageContainer;
   typedef MessageContainer<PairwiseTriplet13Message, 1, 3, variableMessageNumber, 1, variableMessageSize, FMC_GM_PARAM, 3> PairwiseTriplet13MessageContainer;
   typedef MessageContainer<PairwiseTriplet23Message, 1, 3, variableMessageNumber, 1, variableMessageSize, FMC_GM_PARAM, 4> PairwiseTriplet23MessageContainer;

   using FactorList = meta::list< UnaryFactor, PairwiseFactor, EmptyTripletFactor >;
   using MessageList = meta::list<
      UnaryPairwiseMessageLeft,
      UnaryPairwiseMessageRight,
      //UnaryMcfLabelingMessage,
      PairwiseTriplet12MessageContainer, 
      PairwiseTriplet13MessageContainer, 
      PairwiseTriplet23MessageContainer 
         >;

   using mrf = StandardMrfConstructor<FMC_GM_PARAM,0,1,0,1>;
   using tighteningMrf = TighteningMRFProblemConstructor<mrf,2,2,3,4>;
   //using mcfLabeling = MinimumCostFlowLabelingConstructor<FMC_GM_PARAM,0,2,2>;
   using ProblemDecompositionList = meta::list<tighteningMrf>;//,mcfLabeling>;
};


// helper function for extracting types from FMC
template<template <PairwiseConstruction> class FMC, PairwiseConstruction PC>
constexpr PairwiseConstruction FmcConstruction(FMC<PC>)
{
   return PC;
}

template<template <PairwiseConstruction> class FMC1, template <PairwiseConstruction> class FMC2, PairwiseConstruction PC>
constexpr bool FmcTypeCheck(FMC2<PC>)
{
   if(std::is_same<FMC1<PC>,FMC2<PC>>::value) return true;
   else return false;
}


// grammar for reading in files in the format of the Dual Decomposition algorithm of Torresani, Kolmogorov and Rother
namespace TorresaniEtAlInput {
   // here we collect information about graph matching problems. This structure is then given to the respective solver
   struct GraphMatchingInput {
      struct Assignment { INDEX leftNode_, rightNode_; REAL cost_; };
      std::vector<Assignment> assignment_;
      std::vector<std::vector<INDEX> >  leftGraph_, rightGraph_;
      // meaning of tuple: (leftNode1, leftNode2) matched to (rightNode1, rightNode2) with given cost
      std::map<std::tuple<INDEX,INDEX,INDEX,INDEX>, REAL> pairwisePotentials_;
   };

   // first two integers are number of left nodes, number of right nodes, then comes number of assignments, and then number of quadratic terms
   struct init_line : pegtl::seq< opt_whitespace, pegtl::string<'p'>, mand_whitespace, positive_integer, mand_whitespace, positive_integer, mand_whitespace, positive_integer, mand_whitespace, positive_integer, opt_whitespace >{};
   // numbers mean: assignment number (consecutive), then comes left node number, right node number, cost
   struct assignment_line : pegtl::seq< opt_whitespace, pegtl::string<'a'>, mand_whitespace, positive_integer, mand_whitespace, positive_integer, mand_whitespace, positive_integer, mand_whitespace, real_number, opt_whitespace>{};
   // numbers mean: number of left assignment, number of right assignment, cost
   struct quadratic_pot_line : pegtl::seq<opt_whitespace, pegtl::string<'e'>, mand_whitespace, positive_integer, mand_whitespace, positive_integer, mand_whitespace, real_number, opt_whitespace>{};

   // do zrobienia: templatize and put into parse_rules.hxx <- rename
   using comment_line = comment_line<pegtl::string<'c'>>;

   // artifacts from dual decomposition file format
   struct neighbor_line : pegtl::seq< pegtl::sor<pegtl::string<'n','0'>, pegtl::string<'n','1'>>, pegtl::until< pegtl::eol>> {};
   struct coordinate_line : pegtl::seq< pegtl::sor<pegtl::string<'i','0'>, pegtl::string<'i','1'>>, pegtl::until< pegtl::eol>> {};

   // better way to cope with comment lines? On each line there may be a comment
   struct grammar : pegtl::must<
                    pegtl::star<comment_line>,
                    init_line,pegtl::eol,
                    pegtl::star< pegtl::sor<
                       pegtl::seq<quadratic_pot_line,pegtl::eol>,
                       pegtl::seq<assignment_line,pegtl::eol>,
                       comment_line,
                       neighbor_line,
                       coordinate_line,
                       pegtl::seq<opt_whitespace, pegtl::eol>, 
                       opt_whitespace
                    > >, 
                    pegtl::eof
                    > {};

   template< typename FMC, typename Rule >
      struct action
      : pegtl::nothing< Rule > {};

   template<typename FMC> struct action< FMC, positive_integer > {
      static void apply(const pegtl::input & in, Solver<FMC>&, std::stack<SIGNED_INDEX>& integer_stack, std::stack<REAL>&, GraphMatchingInput&)
      { 
         integer_stack.push(std::stoul(in.string())); 
      };
   };

   template<typename FMC> struct action< FMC, real_number > {
      static void apply(const pegtl::input & in, Solver<FMC>&, std::stack<SIGNED_INDEX>&, std::stack<REAL>& real_stack, GraphMatchingInput&)
      {
         real_stack.push(std::stod(in.string()));
      }
   };

   template<typename FMC> struct action< FMC, init_line > {
      static void apply(const pegtl::input& in, Solver<FMC>&, std::stack<SIGNED_INDEX>& integer_stack, std::stack<REAL>& real_stack, GraphMatchingInput& gmInput)
      {
         const INDEX quadraticNo = integer_stack.top(); integer_stack.pop();
         const INDEX assignmentNo = integer_stack.top(); integer_stack.pop();
         const INDEX rightFactorNo = integer_stack.top(); integer_stack.pop();
         const INDEX leftFactorNo = integer_stack.top(); integer_stack.pop();

         gmInput.leftGraph_.resize(leftFactorNo);
         gmInput.rightGraph_.resize(rightFactorNo);
         gmInput.assignment_.reserve(assignmentNo);
         //gmInput.pairwisePotentials_.reserve(quadraticNo);
      }
   };

   template<typename FMC> struct action< FMC, assignment_line > {
      static void apply(const pegtl::input& in, Solver<FMC>&, std::stack<SIGNED_INDEX>& integer_stack, std::stack<REAL>& real_stack, GraphMatchingInput& gmInput)
      {
         const INDEX rightNode = integer_stack.top();
         integer_stack.pop();
         const INDEX leftNode = integer_stack.top();
         integer_stack.pop();
         const INDEX assignmentNo = integer_stack.top();
         integer_stack.pop();
         const REAL cost = real_stack.top();
         real_stack.pop();

         assert(assignmentNo == gmInput.assignment_.size());
         gmInput.assignment_.push_back({leftNode,rightNode,cost});

         assert(leftNode < gmInput.leftGraph_.size());
         gmInput.leftGraph_[leftNode].push_back(rightNode);

         assert(rightNode < gmInput.rightGraph_.size());
         gmInput.rightGraph_[rightNode].push_back(leftNode);
      }
   };
   template<typename FMC> struct action< FMC, quadratic_pot_line > {
      static void apply(const pegtl::input & in, Solver<FMC>&, std::stack<SIGNED_INDEX>& integer_stack, std::stack<REAL>& real_stack, GraphMatchingInput& gmInput)
      {
         const INDEX assignment2 = integer_stack.top();
         integer_stack.pop();
         const INDEX assignment1 = integer_stack.top();
         integer_stack.pop();
         const REAL cost = real_stack.top();
         real_stack.pop();

         const INDEX leftNode1 = gmInput.assignment_[assignment1].leftNode_;
         const INDEX leftNode2 = gmInput.assignment_[assignment2].leftNode_;
         const INDEX rightNode1 = gmInput.assignment_[assignment1].rightNode_;
         const INDEX rightNode2 = gmInput.assignment_[assignment2].rightNode_;

         assert(gmInput.pairwisePotentials_.find( std::make_tuple(leftNode1,leftNode2,rightNode1,rightNode2) ) == gmInput.pairwisePotentials_.end() );
         gmInput.pairwisePotentials_.insert( std::make_pair(std::make_tuple(leftNode1,leftNode2,rightNode1,rightNode2),cost) );
      }
   };


   template<typename MAP>
   void AddQuadraticPotential(MAP& q, INDEX node1, INDEX node2, INDEX oppositeNode1, INDEX oppositeNode2, REAL cost, const std::vector<std::vector<INDEX>>& graph)
   {
      INDEX index1 = std::find(graph[node1].begin(), graph[node1].end(), oppositeNode1) - graph[node1].begin();
      INDEX index2 = std::find(graph[node2].begin(), graph[node2].end(), oppositeNode2) - graph[node2].begin();
      assert(index1 < graph[node1].size());
      assert(index2 < graph[node2].size());

      // always assume that node1 < node2 in q
      //assert(node1 != node2);
      if(node1 == node2) {
         std::cout << "This value is not useful, due to matching constraint\n";
         cost = 1e10;
      }
      if(node1>node2) {
         std::swap(node1,node2);
         std::swap(index1,index2);
      }

      const INDEX leftDim = graph[node1].size()+1;
      const INDEX rightDim = graph[node2].size()+1;

      // initialize empty pairwise cost with infinity on diagonal
      if(q.find(std::make_pair(node1,node2)) == q.end()) {
         std::vector<REAL> cost(leftDim * rightDim, 0.0);
         
         for(INDEX i1=0; i1<graph[node1].size(); ++i1) {
            for(INDEX i2=0; i2<graph[node2].size(); ++i2) {
               if(graph[node1][i1] == graph[node2][i2]) {
                  assert(i1 + i2*leftDim < cost.size());
                  cost[i1 + i2*leftDim] = std::numeric_limits<REAL>::max();
                  cost[i1 + i2*leftDim] = 10000.0; // do zrobienia: make max again. However tightening only works with finite numbers
               }
            }
         }
         q.insert(std::make_pair(std::make_pair(node1,node2), std::move(cost)));
      }

      // now add specific cost
      assert(q.find(std::make_pair(node1,node2)) != q.end());
      std::vector<REAL>& costVec = (*q.find(std::make_pair(node1,node2))).second;
      assert(costVec[index1 + index2*leftDim] == 0.0);
      costVec[index1 + index2*leftDim] = cost;
   }

   std::map<std::pair<INDEX,INDEX>, std::vector<REAL>> BuildLeftPairwisePotentials(const GraphMatchingInput& gmInput, const REAL weight = 1.0)
   {
      std::map<std::pair<INDEX,INDEX>, std::vector<REAL>> q;
      for(auto i : gmInput.pairwisePotentials_) {
         const INDEX leftNode1 = std::get<0>(i.first);
         const INDEX leftNode2 = std::get<1>(i.first);
         const INDEX rightNode1 = std::get<2>(i.first);
         const INDEX rightNode2 = std::get<3>(i.first);
         const REAL cost = i.second;

         const INDEX rightIndex1 = std::find(gmInput.leftGraph_[leftNode1].begin(), gmInput.leftGraph_[leftNode1].end(), rightNode1) - gmInput.leftGraph_[leftNode1].begin();
         const INDEX rightIndex2 = std::find(gmInput.leftGraph_[leftNode2].begin(), gmInput.leftGraph_[leftNode2].end(), rightNode2) - gmInput.leftGraph_[leftNode2].begin();

         AddQuadraticPotential(q, leftNode1,leftNode2, rightIndex1,rightIndex2, weight*cost, gmInput.leftGraph_);
      }
      return q;
   }
   std::map<std::pair<INDEX,INDEX>, std::vector<REAL>> BuildRightPairwisePotentials(const GraphMatchingInput& gmInput, const REAL weight = 1.0)
   {
      std::map<std::pair<INDEX,INDEX>, std::vector<REAL>> q;
      for(auto i : gmInput.pairwisePotentials_) {
         const INDEX leftNode1 = std::get<0>(i.first);
         const INDEX leftNode2 = std::get<1>(i.first);
         const INDEX rightNode1 = std::get<2>(i.first);
         const INDEX rightNode2 = std::get<3>(i.first);
         const REAL cost = i.second;

         const INDEX leftIndex1 = std::find(gmInput.rightGraph_[rightNode1].begin(), gmInput.rightGraph_[rightNode1].end(), leftNode1) - gmInput.rightGraph_[rightNode1].begin();
         const INDEX leftIndex2 = std::find(gmInput.rightGraph_[rightNode2].begin(), gmInput.rightGraph_[rightNode2].end(), leftNode2) - gmInput.rightGraph_[rightNode2].begin();

         AddQuadraticPotential(q, rightNode1,rightNode2, leftIndex1,leftIndex2, weight*cost, gmInput.rightGraph_);
      }
      return q;
   }

   

   // specializations for FMC_MP, FMC_MCF and FMC_GM. Here construction of the model takes place
   template<typename FMC> struct action<FMC, typename std::enable_if<FmcTypeCheck<FMC_MP>(FMC{}) || FmcTypeCheck<FMC_MP_T>(FMC{}), pegtl::eof>::type > {
      static void apply(const pegtl::input& in, Solver<FMC>& pd, std::stack<SIGNED_INDEX>& integer_stack, std::stack<REAL>& real_stack, GraphMatchingInput& gmInput)
      {
         std::cout << "Parsed problem, now construct mp version\n";

         auto& assignment = pd.template GetProblemConstructor<0>();
         for(auto& a : gmInput.assignment_) {
            assignment.AddAssignment(a.leftNode_,a.rightNode_,a.cost_);
         }
         assignment.AddSlack();
         assignment.Construct(pd);

         constexpr PairwiseConstruction pc = FmcConstruction(FMC{});
         auto& mrfLeft = pd.template GetProblemConstructor<1>();
         auto& mrfRight = pd.template GetProblemConstructor<2>();
         // first register unary factors from assignment problem
         for(INDEX i=0; i<assignment.GetNumberOfLeftFactors(); i++) {
            mrfLeft.RegisterUnaryFactor(i,assignment.GetLeftFactor(i));
         }
         for(INDEX i=0; i<assignment.GetNumberOfRightFactors(); i++) {
            //mrfRight.RegisterUnaryFactor(assignment.GetNumberOfLeftFactors()+i,assignment.GetRightFactor(i));
            mrfRight.RegisterUnaryFactor(i,assignment.GetRightFactor(i));
         }
         constexpr REAL pairwiseWeight = pc == PairwiseConstruction::BothSides ? 0.5 : 1.0;

         if(pc == PairwiseConstruction::BothSides || pc == PairwiseConstruction::Left) {
            std::map<std::pair<INDEX,INDEX>, std::vector<REAL>> leftQuadraticPot = BuildLeftPairwisePotentials(gmInput,pairwiseWeight);
            for(auto& q : leftQuadraticPot) {
               //auto p = mrf.AddPairwiseFactor(q.second);
               //mrf.LinkUnaryPairwiseFactor(assignment.GetLeftFactor(q.first.first), p, assignment.GetLeftFactor(q.first.second));
               auto p = mrfLeft.AddPairwiseFactor(q.first.first, q.first.second, q.second);
            }
         }

         if(pc == PairwiseConstruction::BothSides || pc == PairwiseConstruction::Right) {
            std::map<std::pair<INDEX,INDEX>, std::vector<REAL>> rightQuadraticPot = BuildRightPairwisePotentials(gmInput,pairwiseWeight);
            for(auto& q : rightQuadraticPot) {
               //auto p = mrf.AddPairwiseFactor(q.second);
               //mrf.LinkUnaryPairwiseFactor(assignment.GetRightFactor(q.first.first), p, assignment.GetRightFactor(q.first.second));
               auto p = mrfRight.AddPairwiseFactor(q.first.first, q.first.second, q.second);
            }
         }

         mrfLeft.Construct(pd);
         mrfRight.Construct(pd);
         
         // connect left and right mp to mcf labeling factor
         // construct labeling factor
         //auto& mcf = pd.template GetProblemConstructor<3>();
         //mcf.ConstructLinearAssignmentGraph(gmInput.leftGraph_,true);
         //mcf.LinkLeftUnaries(mrfLeft);
         //mcf.LinkRightUnaries(mrfRight,gmInput.leftGraph_,true);
      }
   };

   template<typename FMC> struct action<FMC, typename std::enable_if<FmcTypeCheck<FMC_MCF>(FMC{}) || FmcTypeCheck<FMC_MCF_T>(FMC{}),pegtl::eof>::type> {
      static void
      apply(const pegtl::input& in, Solver<FMC>& pd, std::stack<SIGNED_INDEX>& integer_stack, std::stack<REAL>& real_stack, GraphMatchingInput& gmInput)
      {
         std::cout << "Parsed problem, now construct mcf version\n";

         constexpr PairwiseConstruction pc = FmcConstruction(FMC{});
         constexpr REAL pairwiseWeight = pc == PairwiseConstruction::BothSides ? 0.5 : 1.0;

         auto& assignment = pd.template GetProblemConstructor<0>();
         for(auto a : gmInput.assignment_) {
            assignment.AddAssignment(a.leftNode_,a.rightNode_,a.cost_);
         }
         assignment.AddSlack();
         assignment.Construct(pd);

         auto& mrfLeft = pd.template GetProblemConstructor<1>();
         auto& mrfRight = pd.template GetProblemConstructor<2>();

         // first build unaries
         // left side
         if(pc == PairwiseConstruction::BothSides || pc == PairwiseConstruction::Left) {
            for(INDEX i=0; i<assignment.GetNoLeftNodes()-1; ++i) { // last node is for non-assignment
               auto edgeList = assignment.GetEdgeList(i);
               auto unary = mrfLeft.AddUnaryFactor(i, std::vector<REAL>(edgeList.size(),0.0));
               // connect unary with mcf factor
               UnaryToAssignmentMessage<FMC::McfCoveringFactor> msg(edgeList);
               typename FMC::UnaryToAssignmentMessageContainer* msgContainer = new typename FMC::UnaryToAssignmentMessageContainer(msg, unary, assignment.GetMinCostFlowFactorContainer(), edgeList.size());
               pd.GetLP().AddMessage(msgContainer);
            }
         }
         // right side
         if(pc == PairwiseConstruction::BothSides || pc == PairwiseConstruction::Right) {
            for(INDEX i=0; i<assignment.GetNoRightNodes()-1; ++i) {
               auto& edgeList = assignment.GetEdgeList(i);
               auto unary = mrfRight.AddUnaryFactor(i, std::vector<REAL>(edgeList.size(),0.0));

               UnaryToAssignmentMessage<FMC::McfCoveringFactor> msg(edgeList);
               typename FMC::UnaryToAssignmentMessageContainer* msgContainer = new typename FMC::UnaryToAssignmentMessageContainer(msg, unary, assignment.GetMinCostFlowFactorContainer(), edgeList.size());
               pd.GetLP().AddMessage(msgContainer);
            }
         }
         // now construct pairwise potentials
         // on left side
         if(pc == PairwiseConstruction::BothSides || pc == PairwiseConstruction::Left) {
            std::map<std::pair<INDEX,INDEX>, std::vector<REAL>> leftQuadraticPot = BuildLeftPairwisePotentials(gmInput,pairwiseWeight);
            for(auto it = leftQuadraticPot.begin(); it!=leftQuadraticPot.end(); ++it) {
               const INDEX leftVar = it->first.first;
               const INDEX rightVar = it->first.second;
               auto p = mrfLeft.AddPairwiseFactor(leftVar,rightVar,it->second);
            }
         }
         // on right side
         if(pc == PairwiseConstruction::BothSides || pc == PairwiseConstruction::Right) {
            std::map<std::pair<INDEX,INDEX>, std::vector<REAL>> rightQuadraticPot = BuildRightPairwisePotentials(gmInput,pairwiseWeight);
            for(auto it = rightQuadraticPot.begin(); it!=rightQuadraticPot.end(); ++it) {
               const INDEX leftVar = it->first.first;
               const INDEX rightVar = it->first.second;
               auto p = mrfRight.AddPairwiseFactor(leftVar,rightVar,it->second);
            }
         }

         // finish mrf construction
         mrfLeft.Construct(pd);
         mrfRight.Construct(pd);
      }
   };

   template<typename FMC> struct action< FMC, typename std::enable_if<FmcTypeCheck<FMC_GM>(FMC{}) || FmcTypeCheck<FMC_GM_T>(FMC{}), pegtl::eof>::type > {
      static void apply(const pegtl::input& in, Solver<FMC>& pd, std::stack<SIGNED_INDEX>& integer_stack, std::stack<REAL>& real_stack, GraphMatchingInput& gmInput)
      {
         constexpr PairwiseConstruction pc = FmcConstruction(FMC {});
         static_assert(pc == PairwiseConstruction::Left || pc == PairwiseConstruction::Right,"");
         
         std::cout << "Parsed problem, now construct gm version\n";

         auto& mrf = pd.template GetProblemConstructor<0>();
         if(pc == PairwiseConstruction::Left) {
            std::vector<std::vector<REAL>> leftCost(gmInput.leftGraph_.size()); // unary costs
            for(INDEX i=0; i<gmInput.leftGraph_.size(); ++i) {
               leftCost[i].resize(gmInput.leftGraph_[i].size()+1); // +1 for non-assignment
            }
            for(auto a : gmInput.assignment_) {
               const INDEX var = a.leftNode_;
               const INDEX label = std::find(gmInput.leftGraph_[var].begin(), gmInput.leftGraph_[var].end(), a.rightNode_) - gmInput.leftGraph_[var].begin();
               leftCost[var][label] = a.cost_;
            }
            for(INDEX i=0; i<leftCost.size(); ++i) {
               mrf.AddUnaryFactor(i, leftCost[i]);
            }
         }
         if(pc == PairwiseConstruction::Right) {
            std::vector<std::vector<REAL>> rightCost(gmInput.rightGraph_.size()); // unary costs
            for(INDEX i=0; i<gmInput.rightGraph_.size(); ++i) {
               rightCost[i].resize(gmInput.rightGraph_[i].size()+1); // +1 for non-assignment
            }
            for(auto a : gmInput.assignment_) {
               const INDEX var = a.rightNode_;
               const INDEX label = std::find(gmInput.rightGraph_[var].begin(), gmInput.rightGraph_[var].end(), a.leftNode_) - gmInput.rightGraph_[var].begin();
               rightCost[var][label] = a.cost_;
            }
            for(INDEX i=0; i<rightCost.size(); ++i) {
               mrf.AddUnaryFactor(i, rightCost[i]);
            }
         }


         if(pc == PairwiseConstruction::Left) {
            std::map<std::pair<INDEX,INDEX>, std::vector<REAL>> leftQuadraticPot = BuildLeftPairwisePotentials(gmInput);
            for(auto& q : leftQuadraticPot) {
               //auto p = mrf.AddPairwiseFactor(q.second);
               //mrf.LinkUnaryPairwiseFactor(mrf.GetUnaryFactor(q.first.first), p, mrf.GetUnaryFactor(q.first.second));
               auto p = mrf.AddPairwiseFactor(q.first.first, q.first.second, q.second);
            }
         }
         if(pc == PairwiseConstruction::Right) {
            std::map<std::pair<INDEX,INDEX>, std::vector<REAL>> rightQuadraticPot = BuildRightPairwisePotentials(gmInput);
            for(auto& q : rightQuadraticPot) {
               //auto p = mrf.AddPairwiseFactor(q.second);
               //mrf.LinkUnaryPairwiseFactor(mrf.GetUnaryFactor(q.first.first), p, mrf.GetUnaryFactor(q.first.second));
               auto p = mrf.AddPairwiseFactor(q.first.first, q.first.second, q.second);
            }
         }

         mrf.Construct(pd);
      }
   };

   // the action class templates for the grammar must only depend on one parameter, hence we take out FMC
   template<typename FMC>
      struct actionSpecialization {
         template<typename RULE> struct type : public action<FMC,RULE> {};
      };

   template<typename FMC>
   bool ParseProblem(const std::string filename, Solver<FMC>& pd)
   {
      std::stack<SIGNED_INDEX> integer_stack;
      std::stack<REAL> real_stack;
      GraphMatchingInput gmInput;

      pegtl::file_parser problem(filename);
      std::cout << "parsing " << filename << "\n";

      //return problem.parse< grammar, actionSpecialization<FMC>::template type >(pd, integer_stack, real_stack, mcInput);
      return problem.parse< grammar, actionSpecialization<FMC>::template type >(pd, integer_stack, real_stack, gmInput);
   }
} // end TorresaniEtAlInput

// the UAI mrf format followed by custom constraints describing an underlying assignment problem.
namespace UAIInput {

   struct init_line : pegtl::seq< opt_whitespace, pegtl::string<'M','A','R','K','O','V'>, opt_whitespace > {};
   struct numberOfVariables_line : pegtl::seq< opt_whitespace, positive_integer, opt_whitespace > {};
   // vector of integers denoting how many labels each variable has
   struct cardinality_line : pegtl::seq< opt_whitespace, positive_integer, pegtl::star< opt_whitespace, positive_integer>, opt_whitespace> {};
   struct numberOfCliques : pegtl::seq< opt_whitespace, positive_integer, opt_whitespace> {};
   // first is the number of variables in the clique, then the actual variables.
   struct cliqueScope_line : pegtl::seq< opt_whitespace, positive_integer, pegtl::plus< opt_whitespace, positive_integer>, opt_whitespace, pegtl::eol> {};
   struct cliqueScopes : pegtl::star<cliqueScope_line> {};
   // a function table is begun by number of entries and then a list of real numbers. Here we record all the values in the real stack
   struct functionTables : pegtl::seq<pegtl::star<pegtl::sor<mand_whitespace, pegtl::eol>>, real_number, pegtl::plus<pegtl::star<pegtl::sor<mand_whitespace, pegtl::eol>>, real_number>> {};

   struct constraints : pegtl::seq< opt_whitespace, pegtl::string<'c','o','n','s','t','r','a','i','n','t','s'>, opt_whitespace, pegtl::eol> {};
   struct variable_label_pair : pegtl::seq< pegtl::string<'('>, pegtl::opt<pegtl::string<'v'>>, positive_integer, pegtl::string<','>, pegtl::opt<pegtl::string<'l'>>, positive_integer, pegtl::string<')'>> {};
   struct constraint_line : pegtl::seq< variable_label_pair, pegtl::star<opt_whitespace, pegtl::string<'+'>, opt_whitespace, variable_label_pair>, pegtl::string<'<','='>, opt_whitespace, positive_integer> {};

   struct grammar : pegtl::seq<
                    init_line, pegtl::eol,
                    numberOfVariables_line, pegtl::eol,
                    cardinality_line, pegtl::eol,
                    numberOfCliques, pegtl::eol,
                    cliqueScopes,
                    functionTables,
                    pegtl::star<pegtl::sor<mand_whitespace,pegtl::eol>>,
                    constraints, 
                    pegtl::star<constraint_line,pegtl::eol>,
                    pegtl::opt<constraint_line>,
                    pegtl::star<pegtl::sor<mand_whitespace,pegtl::eol>>,
                    pegtl::eof> {};


   struct GraphMatchingInput {
      INDEX numberOfVariables_;
      INDEX numberOfCliques_;
      std::vector<INDEX> cardinality_;
      std::vector<std::vector<INDEX>> cliqueScope_;
      std::vector<std::vector<REAL>> functionTable_;
      std::vector<std::vector<std::pair<INDEX,INDEX>>> constraint_;
   };

   template< typename FMC, typename Rule >
      struct action
      : pegtl::nothing< Rule > {};


   template<typename FMC> struct action< FMC, positive_integer > {
      static void apply(const pegtl::input & in, Solver<FMC>&, std::stack<SIGNED_INDEX>& integer_stack, std::stack<REAL>&, GraphMatchingInput &)
      { 
         integer_stack.push(std::stoul(in.string())); 
      }
   };
   template<typename FMC> struct action< FMC, real_number > {
      static void apply(const pegtl::input & in, Solver<FMC>&, std::stack<SIGNED_INDEX>&, std::stack<REAL>& real_stack, GraphMatchingInput&)
      { 
         real_stack.push(std::stod(in.string())); 
      }
   };
   template<typename FMC> struct action< FMC, numberOfVariables_line > {
      static void apply(const pegtl::input & in, Solver<FMC>&, std::stack<SIGNED_INDEX>& integer_stack, std::stack<REAL>& real_stack, GraphMatchingInput& gmInput)
      {
         gmInput.numberOfVariables_ = integer_stack.top();
         integer_stack.pop();
      }
   };
   template<typename FMC> struct action< FMC, numberOfCliques > {
      static void apply(const pegtl::input & in, Solver<FMC>&, std::stack<SIGNED_INDEX>& integer_stack, std::stack<REAL>& real_stack, GraphMatchingInput& gmInput)
      {
         gmInput.numberOfCliques_ = integer_stack.top();
         integer_stack.pop();
      }
   };
   template<typename FMC> struct action< FMC, cardinality_line > {
      static void apply(const pegtl::input & in, Solver<FMC>&, std::stack<SIGNED_INDEX>& integer_stack, std::stack<REAL>& real_stack, GraphMatchingInput& gmInput)
      {
         assert(integer_stack.size() == gmInput.numberOfVariables_);
         while(!integer_stack.empty()) {
            gmInput.cardinality_.push_back(integer_stack.top());
            integer_stack.pop();
         }
         std::reverse(gmInput.cardinality_.begin(), gmInput.cardinality_.end());
      }
   };
   template<typename FMC> struct action< FMC, cliqueScope_line > {
      static void apply(const pegtl::input & in, Solver<FMC>&, std::stack<SIGNED_INDEX>& integer_stack, std::stack<REAL>& real_stack, GraphMatchingInput& gmInput)
      {
         gmInput.cliqueScope_.push_back(std::vector<INDEX>(0));
         while(integer_stack.size() > 1) {
            gmInput.cliqueScope_.back().push_back(integer_stack.top());
            integer_stack.pop();
         }
         std::reverse(gmInput.cliqueScope_.back().begin(), gmInput.cliqueScope_.back().end());
         const INDEX cliqueSize = integer_stack.top();
         integer_stack.pop();
         assert(gmInput.cliqueScope_.back().size() == cliqueSize);
      }
   };
   // reconstruct the function tables
   template<typename FMC> struct action< FMC, functionTables > {
      static void apply(const pegtl::input & in, Solver<FMC>&, std::stack<SIGNED_INDEX>& integer_stack, std::stack<REAL>& real_stack, GraphMatchingInput& gmInput)
      {
         // first read in real stack into array
         std::vector<REAL> values;
         while(!real_stack.empty()) {
            values.push_back(real_stack.top());
            real_stack.pop();
         }
         std::reverse(values.begin(), values.end());

         // now iterate over function tables. first read in number of values of function table, then add them to the function table.
         INDEX functionTableIdx = 0; // at which index does current function table start?
         while( functionTableIdx<values.size() ) {
            // assert that this functionTableSize is actually a real variable
            const INDEX functionTableSize = values[functionTableIdx];
            assert(functionTableIdx+functionTableSize < values.size());
            const INDEX potentialNumber = gmInput.functionTable_.size();
            const INDEX cardinality = gmInput.cliqueScope_[potentialNumber].size();

            gmInput.functionTable_.push_back(std::vector<REAL>(functionTableSize,0.0));
            if(cardinality == 1) { // unary factor
               assert(functionTableSize == gmInput.cardinality_[ gmInput.cliqueScope_[potentialNumber][0] ]);
               for(INDEX label=0; label<functionTableSize; ++label) {
                  gmInput.functionTable_.back().operator[](label) = values[functionTableIdx+1+label];
                  //std::cout << values[functionTableIdx+1+label] << ", ";
               }
               //std::cout << "\n";
            } else if(cardinality == 2) { // pairwise factor
               const INDEX var1 = gmInput.cliqueScope_[potentialNumber][0];
               const INDEX var2 = gmInput.cliqueScope_[potentialNumber][1];
               assert(functionTableSize == gmInput.cardinality_[var1] * gmInput.cardinality_[var2]);
               for(INDEX label1=0; label1<gmInput.cardinality_[var1]; ++label1) {
                  for(INDEX label2=0; label2<gmInput.cardinality_[var2]; ++label2) {
                     // note: we must transpose the matrix, that we have read in
                     gmInput.functionTable_.back().operator[](label1 + label2*gmInput.cardinality_[var1]) = values[functionTableIdx + 1 + label2 + label1*gmInput.cardinality_[var2]];
                  }
               }
            } else {
               assert(false);
               throw std::runtime_error("Only unary and pairwise potentials supported now");
            }
            functionTableIdx += functionTableSize+1;
         }
      }
   };
   template<typename FMC> struct action< FMC, constraint_line > {
      static void apply(const pegtl::input & in, Solver<FMC>&, std::stack<SIGNED_INDEX>& integer_stack, std::stack<REAL>& real_stack, GraphMatchingInput& gmInput)
      {
         //std::cout << "constraint: " << in.string() << "\n";
         const INDEX sum = integer_stack.top();
         integer_stack.pop();
         assert(integer_stack.size() % 2 == 0);
         assert(sum == 1);
         // now read in variable label pairs.
         gmInput.constraint_.push_back(std::vector<std::pair<INDEX,INDEX>>(0));
         while(!integer_stack.empty()) {
            const INDEX label = integer_stack.top(); // do zrobienia: in constraints, label start at 1, but variables start at 0
            integer_stack.pop();
            const INDEX var = integer_stack.top();
            integer_stack.pop();
            gmInput.constraint_.back().push_back(std::make_pair(var,label));
         }
         std::reverse(gmInput.constraint_.back().begin(), gmInput.constraint_.back().begin());
      }
   };

   template<typename MRF_CONSTRUCTOR>
      void BuildMrf(MRF_CONSTRUCTOR& mrf, GraphMatchingInput& gmInput)
      {
         // first input the unaries, as pairwise potentials need them to be able to link to them
         for(INDEX i=0; i<gmInput.numberOfCliques_; ++i) {
            if(gmInput.cliqueScope_[i].size() == 1) {
               const INDEX var = gmInput.cliqueScope_[i][0];
               assert(gmInput.functionTable_[i].size() == gmInput.cardinality_[var]);
               mrf.AddUnaryFactor(var,gmInput.functionTable_[i]);
            } else if(gmInput.cliqueScope_[i].size() > 2) {
               throw std::runtime_error("only pairwise models are accepted currently");
            }
         }
         // now the pairwise potentials. 
         for(INDEX i=0; i<gmInput.numberOfCliques_; ++i) {
            if(gmInput.cliqueScope_[i].size() == 2) {
               const INDEX var1 = gmInput.cliqueScope_[i][0];
               const INDEX var2 = gmInput.cliqueScope_[i][1];
               assert(var1<var2);
               assert(gmInput.functionTable_[i].size() == gmInput.cardinality_[var1]*gmInput.cardinality_[var2]);
               mrf.AddPairwiseFactor(var1,var2,gmInput.functionTable_[i]); // or do we have to transpose the values?
            }
         }
      }

   std::vector<std::vector<INDEX>> BuildGraph(GraphMatchingInput& gmInput)
   {
      std::vector<std::vector<INDEX>> graph(gmInput.numberOfVariables_);
      const INDEX noConstraints = gmInput.constraint_.size();
      for(INDEX i=0; i<gmInput.numberOfVariables_; ++i) {
         graph[i].resize(gmInput.cardinality_[i],noConstraints); // this value signifies that the edge does not participate in any constraint
      }
      for(INDEX c=0; c<gmInput.constraint_.size(); c++) {
         //std::cout << "constraint ";
         for(INDEX j=0; j<gmInput.constraint_[c].size(); ++j) {
            const INDEX var = std::get<0>(gmInput.constraint_[c][j]);
            assert(var < gmInput.numberOfVariables_);
            const INDEX label = std::get<1>(gmInput.constraint_[c][j]); // do zrobienia: this is due to the format: labels start at 1, variables start at 0
            //assert(label < mrf.GetNumberOfLabels(var));
            //std::cout << "(" << var << "," << label << ")+";
            assert(graph[var][label] == noConstraints); // assignment problem only allows node to be present in one constraint
            graph[var][label] = c;
         }
         //std::cout << "\n";
      }
      return graph;
   }

   template<typename FMC> struct action<FMC, typename std::enable_if<FmcTypeCheck<FMC_GM>(FMC{}) || FmcTypeCheck<FMC_GM_T>(FMC{}),pegtl::eof>::type> {
      static void apply(const pegtl::input & in, Solver<FMC>& pd, std::stack<SIGNED_INDEX>& integer_stack, std::stack<REAL>& real_stack, GraphMatchingInput& gmInput)
      {
         static_assert(FmcConstruction(FMC{}) == PairwiseConstruction::Left,"");
         std::cout << "construct gm version\n";
         auto& mrf = pd.template GetProblemConstructor<0>();
         BuildMrf(mrf, gmInput);

         // add additional empty pairwise potentials with infty on diagonals for each constraint, if not already there.
         std::map<std::tuple<INDEX,INDEX>, std::vector<REAL>> pairwisePot;
         for(INDEX c=0; c<gmInput.constraint_.size(); ++c) {
            for(INDEX c1=0; c1<gmInput.constraint_[c].size(); ++c1) {
               for(INDEX c2=0; c2<c1; ++c2) {
                  INDEX var1 = gmInput.constraint_[c][c1].first;
                  INDEX label1 = gmInput.constraint_[c][c1].second;
                  INDEX var2 = gmInput.constraint_[c][c2].first;
                  INDEX label2 = gmInput.constraint_[c][c2].second;
                  assert(var1 != var2);
                  if(var1 > var2) {
                     std::swap(var1,var2);
                     std::swap(label1,label2);
                  }
                  assert(var2 < gmInput.numberOfVariables_);
                  assert(label1 < gmInput.cardinality_[var1]);
                  assert(label2 < gmInput.cardinality_[var2]);

                  if(!mrf.HasPairwiseFactor(var1,var2)) {
                     const INDEX potentialSize = gmInput.cardinality_[var1] * gmInput.cardinality_[var2];
                     if(pairwisePot.find(std::make_tuple(var1,var2)) == pairwisePot.end()) {
                        pairwisePot.insert( std::make_pair(std::make_pair(var1,var2), std::vector<REAL>(potentialSize, 0.0)) );
                     } 
                     auto it = pairwisePot.find(std::make_pair(var1,var2));
                     assert(it != pairwisePot.end());
                     assert(it->second.size() == potentialSize);
                     assert(label1 + label2*gmInput.cardinality_[var1] < it->second.size());
                     it->second.operator[](label1 + label2*gmInput.cardinality_[var1]) = 3e11;
                  } else {
                     const INDEX factorId = mrf.GetPairwiseFactorId(var1,var2);
                     const REAL val = mrf.GetPairwiseValue(factorId,label1,label2);
                     assert(val > 1000000);
                  }
               }
            }
         }
         for(auto it = pairwisePot.cbegin(); it!=pairwisePot.cend(); ++it) {
            const INDEX var1 = std::get<0>(it->first);
            const INDEX var2 = std::get<1>(it->first);
            const std::vector<REAL> pot = it->second;
            mrf.AddPairwiseFactor(var1,var2,pot);
         }
  
         // construct labeling factor
         //auto& mcf = pd.template GetProblemConstructor<1>();
         //mcf.ConstructLinearAssignmentGraph(BuildGraph(gmInput));
         //mcf.LinkUnaries(mrf);
      }
   };

   template<typename FMC> struct action<FMC, typename std::enable_if<FmcTypeCheck<FMC_MP>(FMC{}) || FmcTypeCheck<FMC_MP_T>(FMC{}),pegtl::eof>::type> {
      static void apply(const pegtl::input & in, Solver<FMC>& pd, std::stack<SIGNED_INDEX>& integer_stack, std::stack<REAL>& real_stack, GraphMatchingInput& gmInput)
      {
         static_assert(FmcConstruction(FMC{}) == PairwiseConstruction::Left,"");
         std::cout << "construct mp version\n";
         auto& assignment = pd.template GetProblemConstructor<0>();
         auto& mrfLeft = pd.template GetProblemConstructor<1>();
         auto& mrfRight = pd.template GetProblemConstructor<2>();

         std::vector<std::vector<INDEX>> graph = BuildGraph(gmInput);
         std::vector<std::vector<REAL>> unaryCost(graph.size());
         const INDEX noConstraints = gmInput.constraint_.size();
         // get costs of unaries to give to assignment problem
         for(INDEX i=0; i<gmInput.numberOfCliques_; ++i) {
            if(gmInput.cliqueScope_[i].size() == 1) {
               const INDEX var = gmInput.cliqueScope_[i][0];
               assert(gmInput.functionTable_[i].size() == gmInput.cardinality_[var]);
               for(INDEX j=0; j<gmInput.functionTable_[i].size(); ++j) {
                  unaryCost[var].push_back(gmInput.functionTable_[i][j]);
               }
            } else if(gmInput.cliqueScope_[i].size() > 2) {
               throw std::runtime_error("only pairwise models are accepted currently");
            }
         }
         // now add assignment problem
         for(INDEX i=0; i<graph.size(); ++i) {
            for(INDEX j=0; j<graph[i].size(); ++j) {
               const REAL cost = unaryCost[i][j];
               if(graph[i][j] < noConstraints) {
                  assignment.AddAssignment(i,graph[i][j],cost);
               } else if(graph[i][j] == noConstraints) {
                  //std::cout << "Add left slack for node " << i << " and label " << j << "\n";
                  assignment.AddLeftSlack(i,cost);
               } else {
                  assert(false);
                  throw std::runtime_error("assignment problem not well-defined");
               }
            }
         }
         // currently we have <= in constraints. Add slack because of that.
         for(INDEX i=0; i<noConstraints; ++i) {
            assignment.AddRightSlack(i,0.0);
         }
         assignment.Construct(pd);

         // now add pairwise factors and get the unary ones from the assignment.
         /// first register the unaries, as pairwise potentials need them to be able to link to them
         for(INDEX i=0; i<gmInput.numberOfCliques_; ++i) {
            if(gmInput.cliqueScope_[i].size() == 1) {
               const INDEX var = gmInput.cliqueScope_[i][0];
               assert(gmInput.functionTable_[i].size() == gmInput.cardinality_[var]);
               mrfLeft.RegisterUnaryFactor(var,assignment.GetLeftFactor(var));
            } else if(gmInput.cliqueScope_[i].size() > 2) {
               throw std::runtime_error("only pairwise models are accepted currently");
            }
         }
         for(INDEX i=0; i<assignment.GetNumberOfRightFactors(); ++i) {
               mrfRight.RegisterUnaryFactor(i,assignment.GetRightFactor(i));
         }
         // now the pairwise potentials. 
         for(INDEX i=0; i<gmInput.numberOfCliques_; ++i) {
            if(gmInput.cliqueScope_[i].size() == 2) {
               const INDEX var1 = gmInput.cliqueScope_[i][0];
               const INDEX var2 = gmInput.cliqueScope_[i][1];
               assert(var1<var2);
               assert(gmInput.functionTable_[i].size() == gmInput.cardinality_[var1]*gmInput.cardinality_[var2]);
               mrfLeft.AddPairwiseFactor(var1,var2,gmInput.functionTable_[i]);
            }
         }
         std::cout << "Constructed gm with " << mrfLeft.GetNumberOfVariables() << " unary factors and " << mrfLeft.GetNumberOfPairwiseFactors() << " pairwise factors\n";

         // connect left and right mp to mcf labeling factor
         // construct labeling factor
         //auto& mcf = pd.template GetProblemConstructor<3>();
         //mcf.ConstructLinearAssignmentGraph(graph,false);
         //mcf.LinkLeftUnaries(mrfLeft);
         //mcf.LinkRightUnaries(mrfRight,graph,false);
      }
   };

   template<typename FMC> struct action<FMC, typename std::enable_if<FmcTypeCheck<FMC_MCF>(FMC{}) || FmcTypeCheck<FMC_MCF_T>(FMC{}),pegtl::eof>::type> {
      static void apply(const pegtl::input & in, Solver<FMC>& pd, std::stack<SIGNED_INDEX>& integer_stack, std::stack<REAL>& real_stack, GraphMatchingInput& gmInput)
      {
         static_assert(FmcConstruction(FMC{}) == PairwiseConstruction::Left,"");

         std::cout << "problem has " << gmInput.numberOfVariables_ << " variables and " << gmInput.numberOfCliques_ << " cliques\n";
         std::cout << "construct mcf version\n";
         assert(gmInput.numberOfVariables_ == gmInput.cardinality_.size());
         assert(gmInput.numberOfCliques_ == gmInput.functionTable_.size());

         auto& mrf = pd.template GetProblemConstructor<1>();

         BuildMrf(mrf,gmInput);

         // now check if underlying constraints lead to assignment problem. This will be the case, when no variable label pair appears in two constraints and constraints sum to one (this is ensured by the grammar already)
         // build a bipartite network such each variable is a node on the left side, each constraint is a node on the right side, and variable/label pair is an edge which goes from the variable to the constraint.
         const INDEX noVariables = gmInput.numberOfVariables_;
         const INDEX noConstraints = gmInput.constraint_.size();
         std::cout << "no constraints = " << noConstraints << "\n";
         auto& assignment = pd.template GetProblemConstructor<0>();
         std::vector<std::vector<INDEX>> graph = BuildGraph(gmInput);

         // check if pairwise potentials are infinity, whenever two constraints are present
         for(INDEX c=0; c<gmInput.constraint_.size(); ++c) {
            for(INDEX c1=0; c1<gmInput.constraint_[c].size(); ++c1) {
               for(INDEX c2=0; c2<c1; ++c2) {
                  INDEX var1 = gmInput.constraint_[c][c1].first;
                  INDEX label1 = gmInput.constraint_[c][c1].second;
                  INDEX var2 = gmInput.constraint_[c][c2].first;
                  INDEX label2 = gmInput.constraint_[c][c2].second;
                  assert(var1 != var2);
                  if(var1 > var2) {
                     std::swap(var1,var2);
                     std::swap(label1,label2);
                  }

                  if(mrf.HasPairwiseFactor(var1,var2)) {
                     const INDEX factorId = mrf.GetPairwiseFactorId(var1,var2);
                     const REAL val = mrf.GetPairwiseValue(factorId,label1,label2);
                     //std::cout << "diagonal value = " << val << "\n";
                     assert(val > 1000000);
                  }
               }
            }
         }
         // take care of non-assignment
         for(INDEX i=0; i<graph.size(); ++i) {
            for(INDEX j=0; j<graph[i].size(); ++j) {
               if(graph[i][j] < noConstraints) {
                  assignment.AddAssignment(i,graph[i][j],0.0);
               } else if(graph[i][j] == noConstraints) {
                  //std::cout << "Add left slack for node " << i << " and label " << j << "\n";
                  assignment.AddLeftSlack(i,0.0);
               } else {
                  assert(false);
                  throw std::runtime_error("assignment problem not well-defined");
               }
            }
         }
         // currently we have <= in constraints. Add slack because of that.
         for(INDEX i=0; i<noConstraints; ++i) {
            assignment.AddRightSlack(i,0.0);
         }
         assignment.Construct(pd);
         // now link unaries with assignment problem.
         // do zrobienia: assignment problem has non-matching edge as well. Handle this uniformly with TorresaniEtAlInput
         for(INDEX i=0; i<noVariables; ++i) {
            auto& edgeList = assignment.GetEdgeList(i);
            assert(edgeList.size() == gmInput.cardinality_[i]);
            auto unary = mrf.GetUnaryFactor(i);

            UnaryToAssignmentMessage<1> msg(edgeList); // covering factor is 1 as we only have potentials on left side
            typename FMC::UnaryToAssignmentMessageContainer* msgContainer = new typename FMC::UnaryToAssignmentMessageContainer(msg, unary, assignment.GetMinCostFlowFactorContainer(), edgeList.size());
            pd.GetLP().AddMessage(msgContainer);
         }
      }
   };
   // the action class templates for the grammar must only depend on one parameter, hence we take out FMC
   template<typename FMC>
      struct actionSpecialization {
         template<typename RULE> struct type : public action<FMC,RULE> {};
      };

   template<typename FMC>
   bool ParseProblem(const std::string filename, Solver<FMC>& pd)
   {
      std::stack<SIGNED_INDEX> integer_stack;
      std::stack<REAL> real_stack;
      GraphMatchingInput gmInput;

      pegtl::file_parser problem(filename);
      std::cout << "parsing " << filename << "\n";

      return problem.parse< grammar, actionSpecialization<FMC>::template type >(pd, integer_stack, real_stack, gmInput);
   }
}



#endif // LP_MP_GRAPH_MATCHING_H

