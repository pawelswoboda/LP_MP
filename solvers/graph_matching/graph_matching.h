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
// do zrobienia: min cost flow constructors are not used
#include "../cosegmentation/assignment_via_min_cost_flow_constructor.hxx" // move file to problem_constructors
#include "../cosegmentation/assignment_via_message_passing_problem_constructor.hxx" // move file to problem_constructors

#include "problem_constructors/cycle_inequalities.hxx"

#include "parse_rules.h"

#include <vector>
#include <fstream>

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
// UaiGraphMatchingInput contains the grammar in uai MRF format plus constraints section.

// do zrobienia: remove mcf constructors from FMCs and from include

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

   typedef MessageContainer<EqualityMessage, 0, 0, variableMessageNumber, variableMessageNumber, 1, FMC_MP_PARAM, 0 > AssignmentConstraintMessage;
   typedef MessageContainer<LeftMargMessage, 0, 1, variableMessageNumber, 1, variableMessageSize, FMC_MP_PARAM, 1 > UnaryPairwiseMessageLeft;
   typedef MessageContainer<RightMargMessage, 0, 1, variableMessageNumber, 1, variableMessageSize, FMC_MP_PARAM, 2 > UnaryPairwiseMessageRight;

   typedef FactorContainer<Simplex, ExplicitRepamStorage, FMC_MP_PARAM, 2 > EmptyTripletFactor;
   typedef MessageContainer<PairwiseTriplet12Message, 1, 2, variableMessageNumber, 1, variableMessageSize, FMC_MP_PARAM, 4> PairwiseTriplet12MessageContainer;
   typedef MessageContainer<PairwiseTriplet13Message, 1, 2, variableMessageNumber, 1, variableMessageSize, FMC_MP_PARAM, 5> PairwiseTriplet13MessageContainer;
   typedef MessageContainer<PairwiseTriplet23Message, 1, 2, variableMessageNumber, 1, variableMessageSize, FMC_MP_PARAM, 6> PairwiseTriplet23MessageContainer;

   using FactorList = meta::list< UnaryFactor, PairwiseFactor, EmptyTripletFactor>;//McfLabelingFactor, EmptyTripletFactor >;
   using MessageList = meta::list< 
      AssignmentConstraintMessage,
      UnaryPairwiseMessageLeft,
      UnaryPairwiseMessageRight,
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
   typedef MessageContainer<UnaryToAssignmentMessageCS2<McfCoveringFactor>, 1, 0, 1, variableMessageNumber, variableMessageSize, FMC_MCF_PARAM, 2> UnaryToAssignmentMessageContainer;

   using FactorList = meta::list<MinCostFlowAssignmentFactor, UnaryFactor, PairwiseFactor>;
   using MessageList = meta::list<
       UnaryPairwiseMessageLeft,  
       UnaryPairwiseMessageRight, 
       UnaryToAssignmentMessageContainer
      >;

   //using assignment = AssignmentViaMinCostFlowConstructor<FMC_MCF_PARAM,0>;
   //using mcf = AssignmentConstructor<MinCostFlowConstructorCS2<FMC_MCF_PARAM,0>>;
   using mrf = StandardMrfConstructor<FMC_MCF_PARAM,1,2,0,1>;
   using mrfLeft = mrf;
   using mrfRight = mrf;
   //using ProblemDecompositionList = meta::list<assignment,mrfLeft,mrfRight>;
   using ProblemDecompositionList = meta::list<mrf,mrfLeft,mrfRight>;
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
   typedef MessageContainer<UnaryToAssignmentMessageCS2<McfCoveringFactor>, 1, 0, 1, variableMessageNumber, variableMessageSize, FMC_MCF_PARAM, 2> UnaryToAssignmentMessageContainer;

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

   //using assignment = AssignmentViaMinCostFlowConstructor<FMC_MCF_PARAM,0>;
   using mrf = StandardMrfConstructor<FMC_MCF_PARAM,1,2,0,1>;
   using tighteningMrf = TighteningMRFProblemConstructor<mrf,3,3,4,5>;
   using mrfLeft = tighteningMrf;
   using mrfRight = tighteningMrf;
   using ProblemDecompositionList = meta::list<mrf,mrfLeft,mrfRight>; // put away mrf
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

   typedef MessageContainer<LeftMargMessage, 0, 1, variableMessageNumber, 1, variableMessageSize, FMC_GM_PARAM, 0 > UnaryPairwiseMessageLeft;
   typedef MessageContainer<RightMargMessage, 0, 1, variableMessageNumber, 1, variableMessageSize, FMC_GM_PARAM, 1 > UnaryPairwiseMessageRight;

   using FactorList = meta::list< UnaryFactor, PairwiseFactor>;//, McfLabelingFactor >;
   using MessageList = meta::list< UnaryPairwiseMessageLeft, UnaryPairwiseMessageRight>;//, UnaryMcfLabelingMessage >;

   using mrf = StandardMrfConstructor<FMC_GM_PARAM,0,1,0,1>;
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

   typedef MessageContainer<LeftMargMessage, 0, 1, variableMessageNumber, 1, variableMessageSize, FMC_GM_PARAM, 0 > UnaryPairwiseMessageLeft;
   typedef MessageContainer<RightMargMessage, 0, 1, variableMessageNumber, 1, variableMessageSize, FMC_GM_PARAM, 1 > UnaryPairwiseMessageRight;

   // tightening
   typedef FactorContainer<Simplex, ExplicitRepamStorage, FMC_GM_PARAM, 2 > EmptyTripletFactor;
   typedef MessageContainer<PairwiseTriplet12Message, 1, 3, variableMessageNumber, 1, variableMessageSize, FMC_GM_PARAM, 2> PairwiseTriplet12MessageContainer;
   typedef MessageContainer<PairwiseTriplet13Message, 1, 3, variableMessageNumber, 1, variableMessageSize, FMC_GM_PARAM, 3> PairwiseTriplet13MessageContainer;
   typedef MessageContainer<PairwiseTriplet23Message, 1, 3, variableMessageNumber, 1, variableMessageSize, FMC_GM_PARAM, 4> PairwiseTriplet23MessageContainer;

   using FactorList = meta::list< UnaryFactor, PairwiseFactor, EmptyTripletFactor >;
   using MessageList = meta::list<
      UnaryPairwiseMessageLeft,
      UnaryPairwiseMessageRight,
      PairwiseTriplet12MessageContainer, 
      PairwiseTriplet13MessageContainer, 
      PairwiseTriplet23MessageContainer 
         >;

   using mrf = StandardMrfConstructor<FMC_GM_PARAM,0,1,0,1>;
   using tighteningMrf = TighteningMRFProblemConstructor<mrf,2,2,3,4>;
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
   /* file format
   // Angular parentheses mean that it should be replaced with an integer number,
   // curly parentheses mean a floating point number.
   // Point and assignment id's are integers starting from 0.

   c comment line
   p <N0> <N1> <A> <E>     // # points in the left image, # points in the right image, # assignments, # edges
   a <a> <i0> <i1> {cost}  // specify assignment
   e <a> <b> {cost}        // specify edge

   i0 <id> {xi} {yi}       // optional - specify coordinate of a point in the left image
   i1 <id> {xi} {yi}       // optional - specify coordinate of a point in the left image
   n0 <i> <j>              // optional - specify that points <i> and <j> in the left image are neighbors
   n1 <i> <j>              // optional - specify that points <i> and <j> in the right image are neighbors
   */ 

   // here we collect information about graph matching problems. This structure is then given to the respective solver

   using Parsing::mand_whitespace;
   using Parsing::opt_whitespace;
   using Parsing::positive_integer;
   using Parsing::real_number;

   struct GraphMatchingInput {
      struct Assignment { INDEX left_node_, right_node_; REAL cost_; };
      std::vector<Assignment> assignment_;
      std::vector<std::vector<INDEX> >  leftGraph_, rightGraph_;
      // meaning of tuple: (leftNode1, leftNode2) matched to (rightNode1, rightNode2) with given cost
      std::vector<std::tuple<INDEX,INDEX,INDEX,INDEX,REAL>> pairwise_potentials;
   };

   // first two integers are number of left nodes, number of right nodes, then comes number of assignments, and then number of quadratic terms
   struct no_left_nodes : pegtl::seq< positive_integer > {};
   struct no_right_nodes : pegtl::seq< positive_integer > {};
   struct init_line : pegtl::seq< opt_whitespace, pegtl::string<'p'>, mand_whitespace, no_left_nodes, mand_whitespace, no_right_nodes, mand_whitespace, positive_integer, mand_whitespace, positive_integer, opt_whitespace > {};
   // numbers mean: assignment number (consecutive), then comes left node number, right node number, cost
   struct assignment : pegtl::seq < positive_integer, mand_whitespace, positive_integer, mand_whitespace, positive_integer, mand_whitespace, real_number > {};
   struct assignment_line : pegtl::seq< opt_whitespace, pegtl::string<'a'>, mand_whitespace, assignment, opt_whitespace> {};
   // numbers mean: number of left assignment, number of right assignment, cost
   struct quadratic_pot : pegtl::seq< positive_integer, mand_whitespace, positive_integer, mand_whitespace, real_number > {};
   struct quadratic_pot_line : pegtl::seq<opt_whitespace, pegtl::string<'e'>, mand_whitespace, quadratic_pot, opt_whitespace > {};

   struct comment_line : pegtl::seq< opt_whitespace, pegtl::string<'c'>, pegtl::until< pegtl::eol >> {};

   // artifacts from dual decomposition file format. We do not make use of them.
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

   template< typename Rule >
      struct action
      : pegtl::nothing< Rule > {};

   template<> struct action< no_left_nodes > {
      static void apply(const pegtl::action_input& in, GraphMatchingInput& gmInput)
      {
         gmInput.leftGraph_.resize(std::stoul(in.string()));
      }
   };
    
   template<> struct action< no_right_nodes > {
      static void apply(const pegtl::action_input& in, GraphMatchingInput& gmInput)
      {
         gmInput.rightGraph_.resize(std::stoul(in.string()));
      }
   };
    
   template<> struct action< assignment > {
      static void apply(const pegtl::action_input& in, GraphMatchingInput& gmInput)
      {
         std::istringstream iss(in.string());
         INDEX assignment_no; iss >> assignment_no;
         INDEX left_node; iss >> left_node;
         INDEX right_node; iss >> right_node;
         REAL cost; iss >> cost;

         assert(assignment_no == gmInput.assignment_.size());
         gmInput.assignment_.push_back({left_node,right_node,cost});

         assert(left_node < gmInput.leftGraph_.size());
         gmInput.leftGraph_[left_node].push_back(right_node);

         assert(right_node < gmInput.rightGraph_.size());
         gmInput.rightGraph_[right_node].push_back(left_node);
      }
   };
   template<> struct action< quadratic_pot > {
      static void apply(const pegtl::action_input & in, GraphMatchingInput& gmInput)
      {
         std::istringstream iss(in.string());
         INDEX assignment1; iss >> assignment1;
         INDEX assignment2; iss >> assignment2;
         REAL cost; iss >> cost;

         const INDEX leftNode1 = gmInput.assignment_[assignment1].left_node_;
         const INDEX leftNode2 = gmInput.assignment_[assignment2].left_node_;
         const INDEX rightNode1 = gmInput.assignment_[assignment1].right_node_;
         const INDEX rightNode2 = gmInput.assignment_[assignment2].right_node_;

         gmInput.pairwise_potentials.push_back( std::make_tuple(leftNode1,leftNode2,rightNode1,rightNode2,cost) );
      }
   };

   std::vector<std::vector<REAL>> build_left_unaries(const GraphMatchingInput& gm, const REAL weight = 0.5)
   {
      const INDEX no_left_nodes = std::accumulate(gm.assignment_.begin(), gm.assignment_.end(), 0, [](INDEX no, auto a) { return std::max(no, a.left_node_); }) + 1;
      std::vector<std::vector<REAL>> unaries(no_left_nodes);
      for(auto a : gm.assignment_) {
         unaries[a.left_node_].push_back(weight*a.cost_);
      }
      for(auto& u : unaries) { // label for non-assignment
         u.push_back(0.0);
      }
      return std::move(unaries);
   }

   std::vector<std::vector<REAL>> build_right_unaries(const GraphMatchingInput& gm, const REAL weight = 0.5)
   {
      const INDEX no_right_nodes = std::accumulate(gm.assignment_.begin(), gm.assignment_.end(), 0, [](INDEX no, auto a) { return std::max(no, a.right_node_); }) + 1;
      std::vector<std::vector<REAL>> unaries(no_right_nodes);
      for(auto a : gm.assignment_) {
         unaries[a.right_node_].push_back(weight*a.cost_);
      }
      for(auto& u : unaries) { // label for non-assignment
         u.push_back(0.0);
      }
      return std::move(unaries);
   }

   template<typename MAP>
   void AddQuadraticPotential(MAP& q, INDEX node1, INDEX node2, INDEX oppositeNode1, INDEX oppositeNode2, REAL cost, const std::vector<std::vector<INDEX>>& graph)
   {
      INDEX index1 = oppositeNode1; //std::find(graph[node1].begin(), graph[node1].end(), oppositeNode1) - graph[node1].begin();
      INDEX index2 = oppositeNode2; //std::find(graph[node2].begin(), graph[node2].end(), oppositeNode2) - graph[node2].begin();
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
      for(auto i : gmInput.pairwise_potentials) {
         const INDEX leftNode1 = std::get<0>(i);
         const INDEX leftNode2 = std::get<1>(i);
         const INDEX rightNode1 = std::get<2>(i);
         const INDEX rightNode2 = std::get<3>(i);
         const REAL cost = std::get<4>(i);

         const INDEX rightIndex1 = std::find(gmInput.leftGraph_[leftNode1].begin(), gmInput.leftGraph_[leftNode1].end(), rightNode1) - gmInput.leftGraph_[leftNode1].begin();
         const INDEX rightIndex2 = std::find(gmInput.leftGraph_[leftNode2].begin(), gmInput.leftGraph_[leftNode2].end(), rightNode2) - gmInput.leftGraph_[leftNode2].begin();

         AddQuadraticPotential(q, leftNode1,leftNode2, rightIndex1,rightIndex2, weight*cost, gmInput.leftGraph_);
      }
      return q;
   }
   std::map<std::pair<INDEX,INDEX>, std::vector<REAL>> BuildRightPairwisePotentials(const GraphMatchingInput& gmInput, const REAL weight = 1.0)
   {
      std::map<std::pair<INDEX,INDEX>, std::vector<REAL>> q;
      for(auto i : gmInput.pairwise_potentials) {
         const INDEX leftNode1 = std::get<0>(i);
         const INDEX leftNode2 = std::get<1>(i);
         const INDEX rightNode1 = std::get<2>(i);
         const INDEX rightNode2 = std::get<3>(i);
         const REAL cost = std::get<4>(i);

         const INDEX leftIndex1 = std::find(gmInput.rightGraph_[rightNode1].begin(), gmInput.rightGraph_[rightNode1].end(), leftNode1) - gmInput.rightGraph_[rightNode1].begin();
         const INDEX leftIndex2 = std::find(gmInput.rightGraph_[rightNode2].begin(), gmInput.rightGraph_[rightNode2].end(), leftNode2) - gmInput.rightGraph_[rightNode2].begin();

         AddQuadraticPotential(q, rightNode1,rightNode2, leftIndex1,leftIndex2, weight*cost, gmInput.rightGraph_);
      }
      return q;
   }

   template<typename LEFT_MRF_CONSTRUCTOR>
   void construct_left_mrf(GraphMatchingInput& gmInput, LEFT_MRF_CONSTRUCTOR& left_mrf, const REAL unary_weight, const REAL pairwise_weight) 
   {
      // construct unaries of mrfs
      auto left_unaries = build_left_unaries(gmInput,unary_weight);
      for(auto& unary : left_unaries) {
         auto p = left_mrf.AddUnaryFactor(unary);
      }

      // construct pairwise potentials of mrfs
      if(pairwise_weight > 0) {
         std::map<std::pair<INDEX,INDEX>, std::vector<REAL>> leftQuadraticPot = BuildLeftPairwisePotentials(gmInput,pairwise_weight);
         for(auto& q : leftQuadraticPot) {
            auto p = left_mrf.AddPairwiseFactor(q.first.first, q.first.second, q.second);
         }
      }
   }
   template<typename RIGHT_MRF_CONSTRUCTOR>
   void construct_right_mrf(GraphMatchingInput& gmInput, RIGHT_MRF_CONSTRUCTOR& right_mrf, const REAL unary_weight, const REAL pairwise_weight) 
   {
      // construct unaries of mrfs
      auto right_unaries = build_right_unaries(gmInput,unary_weight);
      for(auto& unary : right_unaries) {
         auto p = right_mrf.AddUnaryFactor(unary);
      }

      // construct pairwise potentials of mrfs
      if(pairwise_weight > 0) {
         std::map<std::pair<INDEX,INDEX>, std::vector<REAL>> rightQuadraticPot = BuildRightPairwisePotentials(gmInput,pairwise_weight);
         for(auto& q : rightQuadraticPot) {
            auto p = right_mrf.AddPairwiseFactor(q.first.first, q.first.second, q.second);
         }
      }
   }

   template<typename FMC>
   void construct_mp(Solver<FMC>& s, GraphMatchingInput& gm_input)
   {
      auto& mrf_left = s.template GetProblemConstructor<1>();
      auto& mrf_right = s.template GetProblemConstructor<2>();
      constexpr PairwiseConstruction pc = FmcConstruction(FMC{});
      if(pc == PairwiseConstruction::Left) {
         construct_left_mrf(gm_input, mrf_left, 0.5, 1.0);
         construct_right_mrf(gm_input, mrf_right, 0.5, 0.0);
      }
      if(pc == PairwiseConstruction::BothSides) {
         construct_left_mrf(gm_input, mrf_left, 0.5, 0.5);
         construct_right_mrf(gm_input, mrf_right, 0.5, 0.5);
      }
      if(pc == PairwiseConstruction::Right) {
         construct_left_mrf(gm_input, mrf_left, 0.5, 0.0);
         construct_right_mrf(gm_input, mrf_right, 0.5, 1.0);
      }

      std::vector<INDEX> left_label_count(mrf_left.GetNumberOfVariables(),0);
      std::vector<INDEX> right_label_count(mrf_right.GetNumberOfVariables(),0);
      for(auto& a : gm_input.assignment_) {
         // get left and right unaries and connect them with a message
         auto *l = mrf_left.GetUnaryFactor(a.left_node_);
         auto *r = mrf_right.GetUnaryFactor(a.right_node_);
         const INDEX left_label = left_label_count[a.left_node_];
         const INDEX right_label = right_label_count[a.right_node_];
         auto* m = new typename FMC::AssignmentConstraintMessage( EqualityMessage(left_label, right_label), l, r, 1);
         s.GetLP().AddMessage(m);
         ++left_label_count[a.left_node_];
         ++right_label_count[a.right_node_];
      }
      mrf_left.Construct(s);
      mrf_right.Construct(s);
   }

   template<typename FMC>
   void construct_mcf(Solver<FMC>& s, GraphMatchingInput& gm_input)
   {
      auto& mrf_left = s.template GetProblemConstructor<1>();
      auto& mrf_right = s.template GetProblemConstructor<2>();
      constexpr PairwiseConstruction pc = FmcConstruction(FMC{});
      if(pc == PairwiseConstruction::Left) {
         construct_left_mrf(gm_input, mrf_left, 1.0, 1.0);
      }
      if(pc == PairwiseConstruction::BothSides) {
         construct_left_mrf(gm_input, mrf_left, 0.5, 0.5);
         construct_right_mrf(gm_input, mrf_right, 0.5, 0.5);
      }
      if(pc == PairwiseConstruction::Right) {
         construct_right_mrf(gm_input, mrf_right, 1.0, 1.0);
      }

      // build assignment problem
      const INDEX no_left_nodes = std::accumulate(gm_input.assignment_.begin(), gm_input.assignment_.end(), 0, [](INDEX no, auto a) { return std::max(no, a.left_node_); }) + 1;
      const INDEX no_right_nodes = std::accumulate(gm_input.assignment_.begin(), gm_input.assignment_.end(), 0, [](INDEX no, auto a) { return std::max(no, a.right_node_); }) + 1;
      const INDEX no_edges = gm_input.assignment_.size() + no_left_nodes + no_right_nodes + 1;

      std::vector<typename MinCostFlowFactorCS2::Edge> edges;
      edges.reserve(no_edges);
      for(auto a : gm_input.assignment_) {
         edges.push_back({a.left_node_, no_left_nodes + a.right_node_, 0, 1, 0.0});
      }
      std::vector<SIGNED_INDEX> demands;
      demands.reserve(no_left_nodes + no_right_nodes + 2);
      for(INDEX i=0; i<no_left_nodes; ++i) {
         edges.push_back({i,no_left_nodes + no_right_nodes + 1, 0, 1, 0.0}); // for non-assignment
         demands.push_back(1);
      }
      for(INDEX i=0; i<no_right_nodes; ++i) {
         edges.push_back({no_left_nodes + no_right_nodes, i, 0, 1, 0.0}); // for non-assignment
         demands.push_back(-1);
      }
      edges.push_back({no_left_nodes + no_right_nodes, no_left_nodes + no_right_nodes + 1, 0, std::max(no_left_nodes, no_right_nodes), 0.0});
      demands.push_back(no_right_nodes);
      demands.push_back(-no_left_nodes);

      auto* f = new typename FMC::MinCostFlowAssignmentFactor( MinCostFlowFactorCS2(edges, demands) );
      auto* mcf = f->GetFactor()->GetMinCostFlowSolver();

      // connect assignment factor with unaries
      if(pc == PairwiseConstruction::Left || pc == PairwiseConstruction::BothSides) {
         for(INDEX i=0; i<no_left_nodes; ++i) {
            auto *u = mrf_left.GetUnaryFactor(i);
            auto *m = new typename FMC::UnaryToAssignmentMessageContainer( UnaryToAssignmentMessageCS2<FMC::McfCoveringFactor>(mcf->StartingArc(i), mcf->NoArcs(i)), u, f, mrf_left.GetNumberOfLabels(i));
            s.GetLP().AddMessage(m);
         }
      }
      if(pc == PairwiseConstruction::Right || pc == PairwiseConstruction::BothSides) {
         for(INDEX i=0; i<no_right_nodes; ++i) {
            auto *u = mrf_right.GetUnaryFactor(i);
            auto *m = new typename FMC::UnaryToAssignmentMessageContainer( UnaryToAssignmentMessageCS2<FMC::McfCoveringFactor>(mcf->StartingArc(i + no_left_nodes), mcf->NoArcs(i + no_left_nodes)), u, f, mrf_right.GetNumberOfLabels(i));
            s.GetLP().AddMessage(m);
         }
      }
   }

   template<typename FMC>
   void construct_gm(Solver<FMC>& s, GraphMatchingInput& gm_input)
   {
      // note: possibly more pairwise potentials have to be added for some problems to ensure uniqueness constraint in matching. For house and hotel datasets this is not needed
      auto& mrf = s.template GetProblemConstructor<0>();
      constexpr PairwiseConstruction pc = FmcConstruction(FMC{});
      if(pc == PairwiseConstruction::Left) {
         construct_left_mrf(gm_input, mrf, 1.0, 1.0);
      }
      if(pc == PairwiseConstruction::Right) {
         construct_right_mrf(gm_input, mrf, 1.0, 1.0);
      }
   }

   
   GraphMatchingInput ParseFile(const std::string& filename)
   {
      GraphMatchingInput gmInput;
      pegtl::file_parser problem(filename);
      std::cout << "parsing " << filename << "\n";

      const bool ret = problem.parse< grammar, action >( gmInput );
      if(!ret) {
         throw std::runtime_error("could not read file " + filename);
      }
      return std::move(gmInput);
   }

   template<typename FMC>
   bool ParseProblemGM(const std::string& filename, Solver<FMC>& s)
   {
      auto input = ParseFile(filename);
      construct_gm( s, input );
      return true;
   }

   template<typename FMC>
   bool ParseProblemMP(const std::string& filename, Solver<FMC>& s)
   {
      auto input = ParseFile(filename);
      construct_mp( s, input );
      return true;
   }

   template<typename FMC>
   bool ParseProblemMCF(const std::string& filename, Solver<FMC>& s)
   {
      auto input = ParseFile(filename);
      construct_mcf( s, input );
      return true;
   }


} // end TorresaniEtAlInput

// the UAI mrf format followed by custom constraints describing an underlying assignment problem.
namespace UaiGraphMatchingInput {

   // first use uai input of mrf constructor, afterwards continue with special constraints section
   using Parsing::opt_whitespace;
   using Parsing::mand_whitespace;
   using Parsing::opt_invisible;
   using Parsing::mand_invisible;
   using Parsing::positive_integer;

   struct matching_init_line : pegtl::seq< opt_whitespace, pegtl::string<'m','a','t','c','h','i','n','g'>, opt_whitespace, pegtl::eol> {};
   struct variable : positive_integer {};
   struct label : pegtl::sor< positive_integer, pegtl::string<'s','l','a','c','k'> > {};
   struct matching_line : pegtl::seq< variable, mand_whitespace, pegtl::until< pegtl::eolf, pegtl::seq< label, opt_whitespace> > > {};


   struct grammar : pegtl::seq<
                    pegtl::until<matching_init_line>, 
                    pegtl::plus<matching_line> > {};

   using matching = std::vector<std::vector<INDEX> >;
   using input = std::tuple<UaiMrfInput::MrfInput, matching>;

   template< typename Rule >
      struct action
      : pegtl::nothing< Rule > {};

   template<> struct action< variable > {
      static void apply(const pegtl::action_input & in, matching & m)
      { 
         const INDEX var = std::stoul(in.string());
         assert(m.size() == var);
         m.push_back(std::vector<INDEX>(0));
      }
   };

   template<> struct action< label > {
      static void apply(const pegtl::action_input & in, matching & m)
      { 
         const std::string slack = "slack";
         if(slack == in.string()) {
            m.back().push_back(std::numeric_limits<INDEX>::max());
         } else {
            const INDEX label = std::stoul(in.string());
            assert(label < 100); // do zrobienia: remove this
            m.back().push_back(label);
         }
      }
   };

   input ParseFile(const std::string& filename)
   {
      UaiMrfInput::MrfInput gm_input;

      pegtl::file_parser problem(filename);
      bool ret = problem.parse< UaiMrfInput::grammar, UaiMrfInput::action >( gm_input );
      if(!ret) {
         throw std::runtime_error("could not read mrf input"); 
      }

      matching m;
      ret = problem.parse< grammar, action >( m );
      if(!ret) {
         throw std::runtime_error("could not read matching input"); 
      }
      return std::move(std::make_tuple(std::move(gm_input), std::move(m)));
   }

   // build datastructure indexed by labels and whose elements are the nodes in the gm
   using constraints = std::vector<std::vector<std::tuple<INDEX,INDEX>>>; // variable label pairs which match to specific label
   constraints invert_matching(const matching& m)
   {
      constraints m_inv;
      for(INDEX i=0; i<m.size(); ++i) {
         for(INDEX j=0; j<m[i].size(); ++j) { // first number is the variable
            if(m[i][j] != std::numeric_limits<INDEX>::max() && m[i][j] >= m_inv.size()) {
               m_inv.resize(m[i][j]+1);
            }
            if(m[i][j] != std::numeric_limits<INDEX>::max()) {
               m_inv[m[i][j]].push_back(std::make_tuple(i,j));
            }
         }
      }
      return std::move(m_inv);
   }

   template<typename FMC>
   bool ParseProblemGM(const std::string& filename, Solver<FMC>& s)
   {
      static_assert(FmcConstruction(FMC{}) == PairwiseConstruction::Left, "in uai format only left construction makes sense"); 
      auto i = ParseFile(filename);
      auto& mrf = s.template GetProblemConstructor<0>();
      auto& mrf_input = std::get<0>(i);
      UaiMrfInput::build_mrf(mrf, std::get<0>(i));

      // add additional empty pairwise potentials with infty on diagonals for each constraint, if not already there.

      auto& m = std::get<1>(i);
      auto constraints = invert_matching(m);
      std::map<std::tuple<INDEX,INDEX>, std::vector<REAL>> pairwisePot;
      for(INDEX c=0; c<constraints.size(); ++c) {
         for(INDEX c1=0; c1<constraints[c].size(); ++c1) {
            for(INDEX c2=0; c2<c1; ++c2) {
               INDEX var1 = std::get<0>(constraints[c][c1]);
               INDEX label1 = std::get<1>(constraints[c][c1]);
               INDEX var2 = std::get<0>(constraints[c][c2]);
               INDEX label2 = std::get<1>(constraints[c][c2]);
               assert(var1 != var2);
               if(var1 > var2) {
                  std::swap(var1,var2);
                  std::swap(label1,label2);
               }
               assert(var2 < mrf_input.number_of_variables_);
               assert(label1 < mrf_input.cardinality_[var1]);
               assert(label2 < mrf_input.cardinality_[var2]);

               if(!mrf.HasPairwiseFactor(var1,var2)) {
                  const INDEX potentialSize = mrf_input.cardinality_[var1] * mrf_input.cardinality_[var2];
                  if(pairwisePot.find(std::make_tuple(var1,var2)) == pairwisePot.end()) {
                     pairwisePot.insert( std::make_pair(std::make_pair(var1,var2), std::vector<REAL>(potentialSize, 0.0)) );
                  } 
                  auto it = pairwisePot.find(std::make_pair(var1,var2));
                  assert(it != pairwisePot.end());
                  assert(it->second.size() == potentialSize);
                  assert(label1 + label2*mrf_input.cardinality_[var1] < it->second.size());
                  it->second.operator[](label1 + label2*mrf_input.cardinality_[var1]) = 3e11;
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
      return true;
   }

   template<typename FMC>
   void construct_mp(Solver<FMC>& s, const input& i)
   {
      auto& mrf_left = s.template GetProblemConstructor<1>();
      auto& mrf_input = std::get<0>(i);
      UaiMrfInput::build_mrf(mrf_left, mrf_input);

      // now build unaries for mrf_right. There will be as many unaries (=labels) on the right as constraints
      auto& mrf_right = s.template GetProblemConstructor<2>();
      auto& matching = std::get<1>(i);
      auto constraints = invert_matching(matching);
      for(auto& c : constraints) {
         const INDEX unary_no = mrf_right.AddUnaryFactor(std::vector<REAL>(c.size()+1,0.0)); // extra label is for non-assignment of label
         auto* u_r = mrf_right.GetUnaryFactor(unary_no);
         for(INDEX var_label=0; var_label<c.size(); ++var_label) {
            const INDEX var = std::get<0>(c[var_label]);
            const INDEX label = std::get<1>(c[var_label]);
            auto* u_l = mrf_left.GetUnaryFactor(var);
            auto* m = new typename FMC::AssignmentConstraintMessage( EqualityMessage(label, var_label), u_l, u_r, 1);
            s.GetLP().AddMessage(m);
         }
      }
      
      std::cout << "Constructed gm with " << mrf_left.GetNumberOfVariables() << " unary factors and " << mrf_left.GetNumberOfPairwiseFactors() << " pairwise factors\n";
   }

   template<typename FMC>
   void construct_mcf(Solver<FMC>& s, const input& i)
   {
      auto& mrf_left = s.template GetProblemConstructor<1>();
      const auto& mrf_input = std::get<0>(i);
      UaiMrfInput::build_mrf(mrf_left, mrf_input);

      // build assignment problem
      // We have two types of nodes: matching nodes and slack nodes, the latter taking care whenever a label says slack. Because CS2 orders edges, we must insert slack nodes after the matching nodes, when the need arises. Hence we may have to shift matchign node numbers after constructing slack nodes and inserting them between the matching nodes.
      const auto& matching = std::get<1>(i);
      for(const auto& m : matching) {
         // do zrobienia: this should not be done when assert is not called
         std::vector<INDEX> m_filtered;
         std::copy_if(m.begin(), m.end(), std::back_inserter(m_filtered), [](auto i) { return i != std::numeric_limits<INDEX>::max(); });
         assert(std::is_sorted( m_filtered.begin(), m_filtered.end()));
      }
      const INDEX no_left_nodes = mrf_left.GetNumberOfVariables();
      INDEX no_right_nodes = 0;
      for(auto& m : matching) {
         for(auto c : m) {
            if(c < std::numeric_limits<INDEX>::max()) { // is not a slack node
               no_right_nodes = std::max(no_right_nodes, c+1);
            }
         }
      }
      std::vector<INDEX> slack_node_count(no_right_nodes+1,0); // denotes number of slack nodes that come before normal matching node
      for(INDEX i=0; i<matching.size(); ++i) {
         SIGNED_INDEX last_matching_node = -1;
         for(INDEX j=0; j<matching[i].size(); ++j) {
            if(matching[i][j] == std::numeric_limits<INDEX>::max()) { // is a slack node
               // increase number of slack nodes coming right after right node matching[i][j-1]
               slack_node_count[last_matching_node+1]++;
            } else {
               last_matching_node = matching[i][j];
            }
         }
      }
         
      std::vector<INDEX> matching_node_idx(no_right_nodes+1);
      std::vector<INDEX> cum_slack_node_count(no_right_nodes+1,0);
      std::partial_sum(slack_node_count.begin(), slack_node_count.end(), cum_slack_node_count.begin()); 
      for(INDEX i=0; i<matching_node_idx.size(); ++i) {
         matching_node_idx[i] = i + cum_slack_node_count[i];
      }
      const INDEX total_no_right_nodes = no_right_nodes + std::accumulate(slack_node_count.begin(), slack_node_count.end(), 0);
      std::vector<INDEX> slack_nodes_used(no_right_nodes+1,0);
      std::vector<typename MinCostFlowFactorCS2::Edge> edges;
      for(INDEX i=0; i<matching.size(); ++i) {
         SIGNED_INDEX last_matching_node = -1;
         for(INDEX j=0; j<matching[i].size(); ++j) {
            if(matching[i][j] == std::numeric_limits<INDEX>::max()) { // is a slack node
               const INDEX right_node_no = matching_node_idx[ last_matching_node+1 ] + slack_nodes_used[ last_matching_node+1 ] - slack_node_count[ last_matching_node+1 ];
               ++(slack_nodes_used[ last_matching_node+1 ]);
               edges.push_back({i, no_left_nodes + right_node_no, 0, 1, 0.0});
            } else {
               edges.push_back({i, no_left_nodes + matching_node_idx[ matching[i][j] ], 0, 1, 0.0});
               last_matching_node = matching[i][j];
            }
         }
      }
      for(INDEX i=0; i<total_no_right_nodes; ++i) {
         edges.push_back({no_left_nodes + i, no_left_nodes + total_no_right_nodes, 0, 1, 0.0});
      }

      std::vector<SIGNED_INDEX> demands(no_left_nodes + total_no_right_nodes + 1);
      std::fill(demands.begin(), demands.begin() + no_left_nodes, 1);
      std::fill(demands.begin() + no_left_nodes, demands.end(), 0);
      demands.back() = -no_left_nodes;

      auto* f = new typename FMC::MinCostFlowAssignmentFactor( MinCostFlowFactorCS2(edges, demands) );
      s.GetLP().AddFactor(f);
      auto* mcf = f->GetFactor()->GetMinCostFlowSolver();

      for(INDEX i=0; i<no_left_nodes; ++i) {
         auto *u = mrf_left.GetUnaryFactor(i);
         auto *m = new typename FMC::UnaryToAssignmentMessageContainer( UnaryToAssignmentMessageCS2<FMC::McfCoveringFactor>(mcf->StartingArc(i), mcf->NoArcs(i)), u, f, mrf_left.GetNumberOfLabels(i));
         s.GetLP().AddMessage(m);
      }

      /*

      for(INDEX right_node=0; right_node<constraints.size(); ++right_node) {
         for(INDEX var_label=0; var_label<constraints[right_node].size(); ++var_label) {
            const INDEX var = std::get<0>(constraints[right_node][var_label]);
            const INDEX label = std::get<1>(constraints[right_node][var_label]);
            edges.push_back({var,right_node,0,1,0.0});
         }
      }
      std::vector<SIGNED_INDEX> demands;
      demands.reserve(no_left_nodes + no_right_nodes + 2);
      for(INDEX i=0; i<no_left_nodes; ++i) {
         edges.push_back({i,no_left_nodes + no_right_nodes + 1, 0, 1, 0.0}); // for non-assignment
         demands.push_back(1);
      }
      for(INDEX i=0; i<no_right_nodes; ++i) {
         edges.push_back({no_left_nodes + no_right_nodes, i, 0, 1, 0.0}); // for non-assignment
         demands.push_back(-1);
      }
      edges.push_back({no_left_nodes + no_right_nodes, no_left_nodes + no_right_nodes + 1, 0, std::max(no_left_nodes, no_right_nodes), 0.0});
      demands.push_back(no_right_nodes);
      demands.push_back(-no_left_nodes);

      auto* f = new typename FMC::MinCostFlowAssignmentFactor( MinCostFlowFactorCS2(edges, demands) );
      auto* mcf = f->GetFactor()->GetMinCostFlowSolver();

      // do zrobienia: is this correct? this means that all labels have been used in constraints
      // connect assignment factor with unaries
      for(INDEX i=0; i<no_left_nodes; ++i) {
         auto *u = mrf_left.GetUnaryFactor(i);
         auto *m = new typename FMC::UnaryToAssignmentMessageContainer( UnaryToAssignmentMessageCS2<FMC::McfCoveringFactor>(mcf->StartingArc(i), mcf->NoArcs(i)), u, f, mrf_left.GetNumberOfLabels(i));
         s.GetLP().AddMessage(m);
      }
      */
   }



   
   template<typename FMC>
   bool ParseProblemMP(const std::string& filename, Solver<FMC>& s)
   {
      // do zrobienia: FMC must be left type -> static_assert this
      static_assert(FmcConstruction(FMC{}) == PairwiseConstruction::Left, "in uai format only left construction makes sense"); 
      const auto input = ParseFile(filename);
      construct_mp(s, input);
      return true;
   }

   template<typename FMC>
   bool ParseProblemMCF(const std::string& filename, Solver<FMC>& s)
   {
      static_assert(FmcConstruction(FMC{}) == PairwiseConstruction::Left, "in uai format only left construction makes sense"); 
      const auto input = ParseFile(filename);
      construct_mcf(s, input);
      return true;
   }


   namespace old_format {

      struct constraints_init_line : pegtl::seq< opt_whitespace, pegtl::string<'c','o','n','s','t','r','a','i','n','t','s'>, opt_whitespace, pegtl::eol> {};
      struct variable : positive_integer {};
      struct label : positive_integer {};
      struct variable_label_pair : pegtl::seq< pegtl::string<'('>, pegtl::opt<pegtl::string<'v'>>, variable, pegtl::string<','>, pegtl::opt<pegtl::string<'l'>>, label, pegtl::string<')'>> {};
      struct sum : positive_integer {};
      struct constraint_begin : opt_whitespace {};
      struct constraint_line : pegtl::seq< constraint_begin, variable_label_pair, pegtl::star<opt_whitespace, pegtl::string<'+'>, opt_whitespace, variable_label_pair>, opt_whitespace, pegtl::string<'<','='>, opt_whitespace, sum > {};

      struct grammar : pegtl::seq<
                       pegtl::until<constraints_init_line>, 
                       pegtl::star<constraint_line,pegtl::eol>,
                       pegtl::opt<constraint_line>,
                       pegtl::star<pegtl::sor<mand_whitespace,pegtl::eol>>,
                       pegtl::eof> {};


      using constraints = std::vector<std::vector<std::pair<INDEX,INDEX>>>;
      using input = std::tuple<UaiMrfInput::MrfInput, constraints>;

      template< typename Rule >
         struct action
         : pegtl::nothing< Rule > {};

      template<> struct action< variable > {
         static void apply(const pegtl::action_input & in, constraints & c)
         { 
            //std::cout << "variable \n";
            c.back().push_back(std::make_pair(std::stoul(in.string()),0));
         }
      };

      template<> struct action< label > {
         static void apply(const pegtl::action_input & in, constraints & c)
         { 
            std::get<1>(c.back().back()) = std::stoul(in.string());
         }
      };

      template<> struct action< constraint_begin > {
         static void apply(const pegtl::action_input & in, constraints& c)
         {
            //std::cout << "new constraint\n";
            c.push_back(std::vector<std::pair<INDEX,INDEX>>(0));
         }
      };

      template<> struct action< sum > {
         static void apply(const pegtl::action_input & in, constraints & c)
         { 
            assert(std::stoul(in.string()) == 1);
         }
      };


      input ParseFile(const std::string& filename)
      {
         UaiMrfInput::MrfInput gm_input;

         pegtl::file_parser problem(filename);
         bool ret = problem.parse< UaiMrfInput::grammar, UaiMrfInput::action >( gm_input );
         if(!ret) {
            throw std::runtime_error("could not read mrf input"); 
         }

         constraints c;
         ret = problem.parse< grammar, action >( c );
         if(!ret) {
            throw std::runtime_error("could not read constraints input"); 
         }

         std::ofstream m_file;
         m_file.open (filename, std::ios::out | std::ios::app);

         m_file << "matching\n";
         std::vector<std::vector<INDEX> > matching(gm_input.number_of_variables_);
         for(INDEX i=0; i<gm_input.number_of_variables_; ++i) {
            matching[i].resize(gm_input.cardinality_[i], std::numeric_limits<INDEX>::max());
         }
         for(INDEX idx=0; idx<c.size(); ++idx) {
            for(auto elem : c[idx]) {
               const INDEX var = std::get<0>(elem);
               const INDEX label = std::get<1>(elem);
               matching[var][label] = idx;
            }
         }
         for(INDEX i=0; i<matching.size(); ++i) {
            m_file << i << " ";
            for(INDEX l=0; l<matching[i].size(); ++l) {
               if(matching[i][l] == std::numeric_limits<INDEX>::max()) {
                  m_file << "slack" << " ";
               } else {
                  m_file << matching[i][l] << " ";
               }
            }
            m_file << "\n";
         }
         m_file.close();

         return std::move(std::make_tuple(std::move(gm_input), std::move(c)));
      }

   template<typename FMC>
      void construct_gm(Solver<FMC>& s, const input& i)
      {
         auto& mrf = s.template GetProblemConstructor<0>();
         auto& mrf_input = std::get<0>(i);
         UaiMrfInput::build_mrf(mrf, std::get<0>(i));

         // add additional empty pairwise potentials with infty on diagonals for each constraint, if not already there.
         auto& constraints = std::get<1>(i);
         std::map<std::tuple<INDEX,INDEX>, std::vector<REAL>> pairwisePot;
         for(INDEX c=0; c<constraints.size(); ++c) {
            for(INDEX c1=0; c1<constraints[c].size(); ++c1) {
               for(INDEX c2=0; c2<c1; ++c2) {
                  INDEX var1 = constraints[c][c1].first;
                  INDEX label1 = constraints[c][c1].second;
                  INDEX var2 = constraints[c][c2].first;
                  INDEX label2 = constraints[c][c2].second;
                  assert(var1 != var2);
                  if(var1 > var2) {
                     std::swap(var1,var2);
                     std::swap(label1,label2);
                  }
                  assert(var2 < mrf_input.number_of_variables_);
                  assert(label1 < mrf_input.cardinality_[var1]);
                  assert(label2 < mrf_input.cardinality_[var2]);

                  if(!mrf.HasPairwiseFactor(var1,var2)) {
                     const INDEX potentialSize = mrf_input.cardinality_[var1] * mrf_input.cardinality_[var2];
                     if(pairwisePot.find(std::make_tuple(var1,var2)) == pairwisePot.end()) {
                        pairwisePot.insert( std::make_pair(std::make_pair(var1,var2), std::vector<REAL>(potentialSize, 0.0)) );
                     } 
                  auto it = pairwisePot.find(std::make_pair(var1,var2));
                  assert(it != pairwisePot.end());
                  assert(it->second.size() == potentialSize);
                  assert(label1 + label2*mrf_input.cardinality_[var1] < it->second.size());
                  it->second.operator[](label1 + label2*mrf_input.cardinality_[var1]) = 3e11;
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
   }

   template<typename FMC>
   void construct_mp(Solver<FMC>& s, const input& i)
   {
      auto& mrf_left = s.template GetProblemConstructor<1>();
      auto& mrf_input = std::get<0>(i);
      UaiMrfInput::build_mrf(mrf_left, mrf_input);

      // now build unaries for mrf_right. There will be as many unaries (=labels) on the right as constraints
      auto& mrf_right = s.template GetProblemConstructor<2>();
      auto& constraints = std::get<1>(i);
      for(auto& c : constraints) {
         const INDEX unary_no = mrf_right.AddUnaryFactor(std::vector<REAL>(c.size()+1,0.0)); // extra label is for non-assignment of label
         auto* u_r = mrf_right.GetUnaryFactor(unary_no);
         for(INDEX var_label=0; var_label<c.size(); ++var_label) {
            const INDEX var = std::get<0>(c[var_label]);
            const INDEX label = std::get<1>(c[var_label]);
            auto* u_l = mrf_left.GetUnaryFactor(var);
            auto* m = new typename FMC::AssignmentConstraintMessage( EqualityMessage(label, var_label), u_l, u_r, 1);
            s.GetLP().AddMessage(m);
         }
      }
      
      std::cout << "Constructed gm with " << mrf_left.GetNumberOfVariables() << " unary factors and " << mrf_left.GetNumberOfPairwiseFactors() << " pairwise factors\n";
   }

   template<typename FMC>
   void construct_mcf(Solver<FMC>& s, const input& i)
   {
      auto& mrf_left = s.template GetProblemConstructor<1>();
      auto& mrf_input = std::get<0>(i);
      UaiMrfInput::build_mrf(mrf_left, mrf_input);

      // build assignment problem
      auto& constraints = std::get<1>(i);
      const INDEX no_left_nodes = mrf_left.GetNumberOfVariables();
      const INDEX no_right_nodes = constraints.size();

      std::vector<typename MinCostFlowFactorCS2::Edge> edges;
      for(INDEX right_node=0; right_node<constraints.size(); ++right_node) {
         for(INDEX var_label=0; var_label<constraints[right_node].size(); ++var_label) {
            const INDEX var = std::get<0>(constraints[right_node][var_label]);
            const INDEX label = std::get<1>(constraints[right_node][var_label]);
            edges.push_back({var,right_node,0,1,0.0});
         }
      }
      std::vector<SIGNED_INDEX> demands;
      demands.reserve(no_left_nodes + no_right_nodes + 2);
      for(INDEX i=0; i<no_left_nodes; ++i) {
         edges.push_back({i,no_left_nodes + no_right_nodes + 1, 0, 1, 0.0}); // for non-assignment
         demands.push_back(1);
      }
      for(INDEX i=0; i<no_right_nodes; ++i) {
         edges.push_back({no_left_nodes + no_right_nodes, i, 0, 1, 0.0}); // for non-assignment
         demands.push_back(-1);
      }
      edges.push_back({no_left_nodes + no_right_nodes, no_left_nodes + no_right_nodes + 1, 0, std::max(no_left_nodes, no_right_nodes), 0.0});
      demands.push_back(no_right_nodes);
      demands.push_back(-no_left_nodes);

      auto* f = new typename FMC::MinCostFlowAssignmentFactor( MinCostFlowFactorCS2(edges, demands) );
      auto* mcf = f->GetFactor()->GetMinCostFlowSolver();

      // do zrobienia: is this correct? this means that all labels have been used in constraints
      // connect assignment factor with unaries
      for(INDEX i=0; i<no_left_nodes; ++i) {
         auto *u = mrf_left.GetUnaryFactor(i);
         auto *m = new typename FMC::UnaryToAssignmentMessageContainer( UnaryToAssignmentMessageCS2<FMC::McfCoveringFactor>(mcf->StartingArc(i), mcf->NoArcs(i)), u, f, mrf_left.GetNumberOfLabels(i));
         s.GetLP().AddMessage(m);
      }
   }

   template<typename FMC>
   bool ParseProblem(const std::string filename, Solver<FMC>& pd)
   {
      ParseFile(filename);
      return true;
   }


   } // end old_format



   /*
      template<typename FMC> struct action< FMC, positive_integer > {
      static void apply(const pegtl::action_input & in, Solver<FMC>&, std::stack<SIGNED_INDEX>& integer_stack, std::stack<REAL>&, GraphMatchingInput &)
      { 
      integer_stack.push(std::stoul(in.string())); 
      }
      };
      template<typename FMC> struct action< FMC, real_number > {
      static void apply(const pegtl::action_input & in, Solver<FMC>&, std::stack<SIGNED_INDEX>&, std::stack<REAL>& real_stack, GraphMatchingInput&)
      { 
      real_stack.push(std::stod(in.string())); 
      }
      };
      template<typename FMC> struct action< FMC, numberOfVariables_line > {
      static void apply(const pegtl::action_input & in, Solver<FMC>&, std::stack<SIGNED_INDEX>& integer_stack, std::stack<REAL>& real_stack, GraphMatchingInput& gmInput)
      {
      gmInput.numberOfVariables_ = integer_stack.top();
      integer_stack.pop();
      }
      };
      template<typename FMC> struct action< FMC, numberOfCliques > {
      static void apply(const pegtl::action_input & in, Solver<FMC>&, std::stack<SIGNED_INDEX>& integer_stack, std::stack<REAL>& real_stack, GraphMatchingInput& gmInput)
      {
      gmInput.numberOfCliques_ = integer_stack.top();
      integer_stack.pop();
      }
      };
      template<typename FMC> struct action< FMC, cardinality_line > {
      static void apply(const pegtl::action_input & in, Solver<FMC>&, std::stack<SIGNED_INDEX>& integer_stack, std::stack<REAL>& real_stack, GraphMatchingInput& gmInput)
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
      static void apply(const pegtl::action_input & in, Solver<FMC>&, std::stack<SIGNED_INDEX>& integer_stack, std::stack<REAL>& real_stack, GraphMatchingInput& gmInput)
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
      static void apply(const pegtl::action_input & in, Solver<FMC>&, std::stack<SIGNED_INDEX>& integer_stack, std::stack<REAL>& real_stack, GraphMatchingInput& gmInput)
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
      static void apply(const pegtl::action_input & in, Solver<FMC>&, std::stack<SIGNED_INDEX>& integer_stack, std::stack<REAL>& real_stack, GraphMatchingInput& gmInput)
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
      static void apply(const pegtl::action_input & in, Solver<FMC>& pd, std::stack<SIGNED_INDEX>& integer_stack, std::stack<REAL>& real_stack, GraphMatchingInput& gmInput)
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
      static void apply(const pegtl::action_input & in, Solver<FMC>& pd, std::stack<SIGNED_INDEX>& integer_stack, std::stack<REAL>& real_stack, GraphMatchingInput& gmInput)
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
      static void apply(const pegtl::action_input & in, Solver<FMC>& pd, std::stack<SIGNED_INDEX>& integer_stack, std::stack<REAL>& real_stack, GraphMatchingInput& gmInput)
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
            assert(false);
            //auto& edgeList = assignment.GetEdgeList(i);
            //assert(edgeList.size() == gmInput.cardinality_[i]);
            //auto unary = mrf.GetUnaryFactor(i);

            //UnaryToAssignmentMessage<1> msg(edgeList); // covering factor is 1 as we only have potentials on left side
            //typename FMC::UnaryToAssignmentMessageContainer* msgContainer = new typename FMC::UnaryToAssignmentMessageContainer(msg, unary, assignment.GetMinCostFlowFactorContainer(), edgeList.size());
            //pd.GetLP().AddMessage(msgContainer);
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
   */
}



#endif // LP_MP_GRAPH_MATCHING_H

