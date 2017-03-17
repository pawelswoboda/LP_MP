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

// disable write primal in constructor, for right side mrf constructor
template<typename BASE>
class disable_write_constructor : public BASE
{
   public:
   using BASE::BASE;
   template<typename STREAM>
   void WritePrimal(STREAM& s) const 
   {} 
};

// graph matching with assignment via message passing
template<PairwiseConstruction PAIRWISE_CONSTRUCTION = PairwiseConstruction::Left>
struct FMC_MP {
   using FMC_MP_PARAM = FMC_MP<PAIRWISE_CONSTRUCTION>;
   constexpr static const char* name =
      PAIRWISE_CONSTRUCTION == PairwiseConstruction::Left ? "AMP-O"
      : (PAIRWISE_CONSTRUCTION == PairwiseConstruction::Right ? "AMP-I"
      : (PAIRWISE_CONSTRUCTION == PairwiseConstruction::BothSides ? "AMP-B"
      : "unknown variant"));
      
   typedef FactorContainer<UnarySimplexFactor, FMC_MP_PARAM, 0, true > UnaryFactor; // set to true if labeling by unaries is desired
   typedef FactorContainer<PairwiseSimplexFactor, FMC_MP_PARAM, 1, false > PairwiseFactor;

   constexpr static const Chirality primal_propagation_direction = (PAIRWISE_CONSTRUCTION == PairwiseConstruction::Left || PAIRWISE_CONSTRUCTION == PairwiseConstruction::BothSides) ? Chirality::right : Chirality::left;
   using EqualityMessageType = EqualityMessage<primal_propagation_direction>;
   typedef MessageContainer<EqualityMessageType, 0, 0, variableMessageNumber, variableMessageNumber, FMC_MP_PARAM, 0 > AssignmentConstraintMessage;
   typedef MessageContainer<UnaryPairwiseMessageLeft<MessageSendingType::SRMP>, 0, 1, variableMessageNumber, 1, FMC_MP_PARAM, 1 > UnaryPairwiseMessageLeftContainer;
   typedef MessageContainer<UnaryPairwiseMessageRight<MessageSendingType::SRMP>, 0, 1, variableMessageNumber, 1, FMC_MP_PARAM, 2 > UnaryPairwiseMessageRightContainer;

   using FactorList = meta::list< UnaryFactor, PairwiseFactor>;
   using MessageList = meta::list< AssignmentConstraintMessage, UnaryPairwiseMessageLeftContainer, UnaryPairwiseMessageRightContainer>;//, UnaryMcfLabelingMessage >;

   using mrf = AssignmentGmConstructor<StandardMrfConstructor<FMC_MP_PARAM,0,1,1,2>>;
   using mrfLeft = mrf;
   using mrfRight = disable_write_constructor<mrf>;
   using ProblemDecompositionList = meta::list<mrfLeft, mrfRight>;
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
      
   typedef FactorContainer<UnarySimplexFactor, FMC_MP_PARAM, 0, true > UnaryFactor; // set to true if labeling by unaries is desired
   typedef FactorContainer<PairwiseSimplexFactor, FMC_MP_PARAM, 1, false > PairwiseFactor;

   constexpr static const Chirality primal_propagation_direction = (PAIRWISE_CONSTRUCTION == PairwiseConstruction::Left || PAIRWISE_CONSTRUCTION == PairwiseConstruction::BothSides) ? Chirality::right : Chirality::left;
   using EqualityMessageType = EqualityMessage<primal_propagation_direction>;
   typedef MessageContainer<EqualityMessageType, 0, 0, variableMessageNumber, variableMessageNumber, FMC_MP_PARAM, 0 > AssignmentConstraintMessage;
   typedef MessageContainer<UnaryPairwiseMessageLeft<MessageSendingType::SRMP>, 0, 1, variableMessageNumber, 1, FMC_MP_PARAM, 1 > UnaryPairwiseMessageLeftContainer;
   typedef MessageContainer<UnaryPairwiseMessageRight<MessageSendingType::SRMP>, 0, 1, variableMessageNumber, 1, FMC_MP_PARAM, 2 > UnaryPairwiseMessageRightContainer;

   typedef FactorContainer<SimpleTighteningTernarySimplexFactor, FMC_MP_PARAM, 2 > EmptyTripletFactor;
   typedef MessageContainer<PairwiseTripletMessage12<MessageSendingType::SRMP>, 1, 2, variableMessageNumber, 1, FMC_MP_PARAM, 3> PairwiseTriplet12MessageContainer;
   typedef MessageContainer<PairwiseTripletMessage13<MessageSendingType::SRMP>, 1, 2, variableMessageNumber, 1, FMC_MP_PARAM, 4> PairwiseTriplet13MessageContainer;
   typedef MessageContainer<PairwiseTripletMessage23<MessageSendingType::SRMP>, 1, 2, variableMessageNumber, 1, FMC_MP_PARAM, 5> PairwiseTriplet23MessageContainer;

   using FactorList = meta::list< UnaryFactor, PairwiseFactor, EmptyTripletFactor>;
   using MessageList = meta::list< 
      AssignmentConstraintMessage,
      UnaryPairwiseMessageLeftContainer,
      UnaryPairwiseMessageRightContainer,
      PairwiseTriplet12MessageContainer, 
      PairwiseTriplet13MessageContainer, 
      PairwiseTriplet23MessageContainer 
         >;

   using mrf = AssignmentGmConstructor<StandardMrfConstructor<FMC_MP_PARAM,0,1,1,2>>;
   using tighteningMrf = TighteningMRFProblemConstructor<mrf,2,3,4,5>;
   using mrfLeft = tighteningMrf;
   using mrfRight = disable_write_constructor<tighteningMrf>;
   using ProblemDecompositionList = meta::list<mrfLeft, mrfRight>;
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

   typedef FactorContainer<MinCostFlowFactorCS2, FMC_MCF_PARAM, 0, false > MinCostFlowAssignmentFactor;
   typedef FactorContainer<UnarySimplexFactor, FMC_MCF_PARAM, 1, true > UnaryFactor;
   typedef FactorContainer<PairwiseSimplexFactor, FMC_MCF_PARAM, 2 > PairwiseFactor;

   typedef MessageContainer<UnaryPairwiseMessageLeft<MessageSendingType::SRMP>, 1, 2, variableMessageNumber, 1, FMC_MCF_PARAM, 0 > UnaryPairwiseMessageLeftContainer;
   typedef MessageContainer<UnaryPairwiseMessageRight<MessageSendingType::SRMP>, 1, 2, variableMessageNumber, 1, FMC_MCF_PARAM, 1 > UnaryPairwiseMessageRightContainer;
   //typedef MessageContainer<LeftMargMessage, 1, 2, variableMessageNumber, 1, FMC_MCF_PARAM, 0 > UnaryPairwiseMessageLeft;
   //typedef MessageContainer<RightMargMessage, 1, 2, variableMessageNumber, 1, FMC_MCF_PARAM, 1 > UnaryPairwiseMessageRight;
   using UnaryToAssignmentMessageType = UnaryToAssignmentMessageCS2<McfCoveringFactor,MessagePassingType::SRMP>;
   typedef MessageContainer<UnaryToAssignmentMessageType, 1, 0, 1, variableMessageNumber, FMC_MCF_PARAM, 2> UnaryToAssignmentMessageContainer;

   constexpr static const Chirality primal_propagation_direction = (PAIRWISE_CONSTRUCTION == PairwiseConstruction::Left || PAIRWISE_CONSTRUCTION == PairwiseConstruction::BothSides) ? Chirality::left : Chirality::right;
   using EqualityMessageType = EqualityMessage<primal_propagation_direction,false>;
   typedef MessageContainer<EqualityMessageType, 1, 1, variableMessageNumber, variableMessageNumber, FMC_MCF_PARAM, 3 > AssignmentConstraintMessage;


   using FactorList = meta::list<MinCostFlowAssignmentFactor, UnaryFactor, PairwiseFactor>;
   using MessageList = meta::list<
       UnaryPairwiseMessageLeftContainer,  
       UnaryPairwiseMessageRightContainer, 
       UnaryToAssignmentMessageContainer,
       AssignmentConstraintMessage
      >;

   //using assignment = AssignmentViaMinCostFlowConstructor<FMC_MCF_PARAM,0>;
   //using mcf = AssignmentConstructor<MinCostFlowConstructorCS2<FMC_MCF_PARAM,0>>;
   //using mrf = StandardMrfConstructor<FMC_MCF_PARAM,1,2,0,1>;
   using mrf = AssignmentGmConstructor<StandardMrfConstructor<FMC_MCF_PARAM,1,2,0,1>>;
   using mrfLeft = mrf;
   using mrfRight = disable_write_constructor<mrf>;
   //using ProblemDecompositionList = meta::list<assignment,mrfLeft,mrfRight>;
   using ProblemDecompositionList = meta::list<mrfLeft,mrfRight>;
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

   typedef FactorContainer<MinCostFlowFactorCS2, FMC_MCF_PARAM, 0, false > MinCostFlowAssignmentFactor;
   typedef FactorContainer<UnarySimplexFactor, FMC_MCF_PARAM, 1, true > UnaryFactor;
   typedef FactorContainer<PairwiseSimplexFactor, FMC_MCF_PARAM, 2 > PairwiseFactor;
   typedef FactorContainer<SimpleTighteningTernarySimplexFactor, FMC_MCF_PARAM, 3 > EmptyTripletFactor;

   typedef MessageContainer<UnaryPairwiseMessageLeft<MessageSendingType::SRMP>, 1, 2, variableMessageNumber, 1, FMC_MCF_PARAM, 0 > UnaryPairwiseMessageLeftContainer;
   typedef MessageContainer<UnaryPairwiseMessageRight<MessageSendingType::SRMP>, 1, 2, variableMessageNumber, 1, FMC_MCF_PARAM, 1 > UnaryPairwiseMessageRightContainer;
   //typedef MessageContainer<LeftMargMessage, 1, 2, variableMessageNumber, 1, variableMessageSize, FMC_MCF_PARAM, 0 > UnaryPairwiseMessageLeft;
   //typedef MessageContainer<RightMargMessage, 1, 2, variableMessageNumber, 1, variableMessageSize, FMC_MCF_PARAM, 1 > UnaryPairwiseMessageRight;
   using UnaryToAssignmentMessageType = UnaryToAssignmentMessageCS2<McfCoveringFactor,MessagePassingType::SRMP>;
   typedef MessageContainer<UnaryToAssignmentMessageType, 1, 0, 1, variableMessageNumber, FMC_MCF_PARAM, 2> UnaryToAssignmentMessageContainer;

   typedef MessageContainer<PairwiseTripletMessage12<MessageSendingType::SRMP>, 2, 3, variableMessageNumber, 1, FMC_MCF_PARAM, 3> PairwiseTriplet12MessageContainer;
   typedef MessageContainer<PairwiseTripletMessage13<MessageSendingType::SRMP>, 2, 3, variableMessageNumber, 1, FMC_MCF_PARAM, 4> PairwiseTriplet13MessageContainer;
   typedef MessageContainer<PairwiseTripletMessage23<MessageSendingType::SRMP>, 2, 3, variableMessageNumber, 1, FMC_MCF_PARAM, 5> PairwiseTriplet23MessageContainer;
   //typedef MessageContainer<PairwiseTriplet12Message, 2, 3, variableMessageNumber, 1, variableMessageSize, FMC_MCF_PARAM, 3> PairwiseTriplet12MessageContainer;
   //typedef MessageContainer<PairwiseTriplet13Message, 2, 3, variableMessageNumber, 1, variableMessageSize, FMC_MCF_PARAM, 4> PairwiseTriplet13MessageContainer;
   //typedef MessageContainer<PairwiseTriplet23Message, 2, 3, variableMessageNumber, 1, variableMessageSize, FMC_MCF_PARAM, 5> PairwiseTriplet23MessageContainer;

   constexpr static const Chirality primal_propagation_direction = (PAIRWISE_CONSTRUCTION == PairwiseConstruction::Left || PAIRWISE_CONSTRUCTION == PairwiseConstruction::BothSides) ? Chirality::left : Chirality::right;
   using EqualityMessageType = EqualityMessage<primal_propagation_direction, false>;
   typedef MessageContainer<EqualityMessageType, 1, 1, variableMessageNumber, variableMessageNumber, FMC_MCF_PARAM, 6 > AssignmentConstraintMessage;

   using FactorList = meta::list<MinCostFlowAssignmentFactor, UnaryFactor, PairwiseFactor, EmptyTripletFactor>;
   using MessageList = meta::list<
       UnaryPairwiseMessageLeftContainer,  
       UnaryPairwiseMessageRightContainer, 
       UnaryToAssignmentMessageContainer, 
       PairwiseTriplet12MessageContainer, 
       PairwiseTriplet13MessageContainer, 
       PairwiseTriplet23MessageContainer,
       AssignmentConstraintMessage
      >;

   using mrf = AssignmentGmConstructor<StandardMrfConstructor<FMC_MCF_PARAM,1,2,0,1>>;
   using tighteningMrf = TighteningMRFProblemConstructor<mrf,3,3,4,5>;
   using mrfLeft = tighteningMrf;
   using mrfRight = disable_write_constructor<tighteningMrf>;
   using ProblemDecompositionList = meta::list<mrfLeft,mrfRight>; 
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
      
   typedef FactorContainer<UnarySimplexFactor, FMC_GM_PARAM, 0, true > UnaryFactor; // make true, if primal rounding similar to TRW-S is required
   typedef FactorContainer<PairwiseSimplexFactor, FMC_GM_PARAM, 1, false > PairwiseFactor;

   typedef MessageContainer<UnaryPairwiseMessageLeft<MessageSendingType::SRMP>, 0, 1, variableMessageNumber, 1, FMC_GM_PARAM, 0 > UnaryPairwiseMessageLeftContainer;
   typedef MessageContainer<UnaryPairwiseMessageRight<MessageSendingType::SRMP>, 0, 1, variableMessageNumber, 1, FMC_GM_PARAM, 1 > UnaryPairwiseMessageRightContainer;
   //typedef MessageContainer<LeftMargMessage, 0, 1, variableMessageNumber, 1, variableMessageSize, FMC_GM_PARAM, 0 > UnaryPairwiseMessageLeft;
   //typedef MessageContainer<RightMargMessage, 0, 1, variableMessageNumber, 1, variableMessageSize, FMC_GM_PARAM, 1 > UnaryPairwiseMessageRight;

   using FactorList = meta::list< UnaryFactor, PairwiseFactor>;//, McfLabelingFactor >;
   using MessageList = meta::list< UnaryPairwiseMessageLeftContainer, UnaryPairwiseMessageRightContainer>;//, UnaryMcfLabelingMessage >;

   using mrf = AssignmentGmConstructor<StandardMrfConstructor<FMC_GM_PARAM,0,1,0,1>>;
   using mrf_left = mrf;
   using mrf_right = mrf;
   using ProblemDecompositionList = meta::list<mrf_left, mrf_right>;//,mcfLabeling>;
};

// + tightening triplets
template<PairwiseConstruction PAIRWISE_CONSTRUCTION = PairwiseConstruction::Left> // note: both sides makes no sense here
struct FMC_GM_T {
   using FMC_GM_PARAM = FMC_GM_T<PAIRWISE_CONSTRUCTION>;
   constexpr static const char* name =
      PAIRWISE_CONSTRUCTION == PairwiseConstruction::Left ? "GM-O-T"
      : (PAIRWISE_CONSTRUCTION == PairwiseConstruction::Right ? "GM-I-T"
      : "unknown variant");
      
   typedef FactorContainer<UnarySimplexFactor, FMC_GM_PARAM, 0, true > UnaryFactor; // make true, if primal rounding similar to TRW-S is required
   typedef FactorContainer<PairwiseSimplexFactor, FMC_GM_PARAM, 1, false > PairwiseFactor;

   typedef MessageContainer<UnaryPairwiseMessageLeft<MessageSendingType::SRMP>, 0, 1, variableMessageNumber, 1, FMC_GM_PARAM, 0 > UnaryPairwiseMessageLeftContainer;
   typedef MessageContainer<UnaryPairwiseMessageRight<MessageSendingType::SRMP>, 0, 1, variableMessageNumber, 1, FMC_GM_PARAM, 1 > UnaryPairwiseMessageRightContainer;
   //typedef MessageContainer<LeftMargMessage, 0, 1, variableMessageNumber, 1, variableMessageSize, FMC_GM_PARAM, 0 > UnaryPairwiseMessageLeft;
   //typedef MessageContainer<RightMargMessage, 0, 1, variableMessageNumber, 1, variableMessageSize, FMC_GM_PARAM, 1 > UnaryPairwiseMessageRight;

   // tightening
   typedef FactorContainer<SimpleTighteningTernarySimplexFactor, FMC_GM_PARAM, 2 > EmptyTripletFactor;
   typedef MessageContainer<PairwiseTripletMessage12<MessageSendingType::SRMP>, 1, 2, variableMessageNumber, 1, FMC_GM_PARAM, 2> PairwiseTriplet12MessageContainer;
   typedef MessageContainer<PairwiseTripletMessage13<MessageSendingType::SRMP>, 1, 2, variableMessageNumber, 1, FMC_GM_PARAM, 3> PairwiseTriplet13MessageContainer;
   typedef MessageContainer<PairwiseTripletMessage23<MessageSendingType::SRMP>, 1, 2, variableMessageNumber, 1, FMC_GM_PARAM, 4> PairwiseTriplet23MessageContainer;
   //typedef MessageContainer<PairwiseTriplet12Message, 1, 2, variableMessageNumber, 1, variableMessageSize, FMC_GM_PARAM, 2> PairwiseTriplet12MessageContainer;
   //typedef MessageContainer<PairwiseTriplet13Message, 1, 2, variableMessageNumber, 1, variableMessageSize, FMC_GM_PARAM, 3> PairwiseTriplet13MessageContainer;
   //typedef MessageContainer<PairwiseTriplet23Message, 1, 2, variableMessageNumber, 1, variableMessageSize, FMC_GM_PARAM, 4> PairwiseTriplet23MessageContainer;

   using FactorList = meta::list< UnaryFactor, PairwiseFactor, EmptyTripletFactor >;
   using MessageList = meta::list<
      UnaryPairwiseMessageLeftContainer,
      UnaryPairwiseMessageRightContainer,
      PairwiseTriplet12MessageContainer, 
      PairwiseTriplet13MessageContainer, 
      PairwiseTriplet23MessageContainer 
         >;

   using mrf = AssignmentGmConstructor<StandardMrfConstructor<FMC_GM_PARAM,0,1,0,1>>;
   using tighteningMrf = TighteningMRFProblemConstructor<mrf,2,2,3,4>;
   using mrf_left = tighteningMrf;
   using mrf_right = tighteningMrf;
   using ProblemDecompositionList = meta::list<mrf_left, mrf_right>;//,mcfLabeling>;
};

template<PairwiseConstruction PAIRWISE_CONSTRUCTION = PairwiseConstruction::Left> // note: both sides makes no sense here
struct FMC_HUNGARIAN_BP {
   using FMC_PARAM = FMC_HUNGARIAN_BP<PAIRWISE_CONSTRUCTION>;

   constexpr static const char* name =
      PAIRWISE_CONSTRUCTION == PairwiseConstruction::Left ? "HBP-O"
      : (PAIRWISE_CONSTRUCTION == PairwiseConstruction::Right ? "HBP-I"
      : (PAIRWISE_CONSTRUCTION == PairwiseConstruction::BothSides ? "HBP-B"
      : "unknown variant"));
      
   constexpr static INDEX McfCoveringFactor = PAIRWISE_CONSTRUCTION == PairwiseConstruction::BothSides ? 2 : 1;

   // rounding is done by the mcf factor
   typedef FactorContainer<MinCostFlowFactorCS2, FMC_PARAM, 0, true > MinCostFlowAssignmentFactor;
   typedef FactorContainer<UnarySimplexFactor, FMC_PARAM, 1, false > UnaryFactor;
   typedef FactorContainer<PairwiseSimplexFactor, FMC_PARAM, 2 > PairwiseFactor;

   typedef MessageContainer<UnaryPairwiseMessageLeft<MessageSendingType::MPLP>, 1, 2, variableMessageNumber, 1, FMC_PARAM, 0 > UnaryPairwiseMessageLeftContainer;
   typedef MessageContainer<UnaryPairwiseMessageRight<MessageSendingType::MPLP>, 1, 2, variableMessageNumber, 1, FMC_PARAM, 1 > UnaryPairwiseMessageRightContainer;
   //typedef MessageContainer<LeftMargMessageMPLP, 1, 2, variableMessageNumber, 1, variableMessageSize, FMC_PARAM, 0 > UnaryPairwiseMessageLeft;
   //typedef MessageContainer<RightMargMessageMPLP, 1, 2, variableMessageNumber, 1, variableMessageSize, FMC_PARAM, 1 > UnaryPairwiseMessageRight;
   using UnaryToAssignmentMessageType = UnaryToAssignmentMessageCS2<McfCoveringFactor,MessagePassingType::HUNGARIAN>;
   typedef MessageContainer<UnaryToAssignmentMessageType, 1, 0, 1, variableMessageNumber, FMC_PARAM, 2> UnaryToAssignmentMessageContainer;

   using FactorList = meta::list<MinCostFlowAssignmentFactor, UnaryFactor, PairwiseFactor>;
   using MessageList = meta::list<
       UnaryPairwiseMessageLeftContainer,  
       UnaryPairwiseMessageRightContainer, 
       UnaryToAssignmentMessageContainer
      >;

   using mrf = AssignmentGmConstructor<StandardMrfConstructor<FMC_HUNGARIAN_BP,1,2,0,1>>;
   using mrfLeft = mrf;
   using mrfRight = disable_write_constructor<mrf>;
   using ProblemDecompositionList = meta::list<mrfLeft,mrfRight>;
};

template<PairwiseConstruction PAIRWISE_CONSTRUCTION = PairwiseConstruction::Left> // note: both sides makes no sense here
struct FMC_HUNGARIAN_BP_T {
   using FMC_PARAM = FMC_HUNGARIAN_BP_T<PAIRWISE_CONSTRUCTION>;

   constexpr static const char* name =
      PAIRWISE_CONSTRUCTION == PairwiseConstruction::Left ? "HBP-O-T"
      : (PAIRWISE_CONSTRUCTION == PairwiseConstruction::Right ? "HBP-I-T"
      : (PAIRWISE_CONSTRUCTION == PairwiseConstruction::BothSides ? "HBP-B-T"
      : "unknown variant"));
      
   constexpr static INDEX McfCoveringFactor = PAIRWISE_CONSTRUCTION == PairwiseConstruction::BothSides ? 2 : 1;

   // rounding is done by the mcf factor
   typedef FactorContainer<MinCostFlowFactorCS2, FMC_PARAM, 0, true > MinCostFlowAssignmentFactor;
   typedef FactorContainer<UnarySimplexFactor, FMC_PARAM, 1, false > UnaryFactor;
   typedef FactorContainer<PairwiseSimplexFactor, FMC_PARAM, 2 > PairwiseFactor;

   typedef MessageContainer<UnaryPairwiseMessageLeft<MessageSendingType::MPLP>, 1, 2, variableMessageNumber, 1, FMC_PARAM, 0 > UnaryPairwiseMessageLeftContainer;
   typedef MessageContainer<UnaryPairwiseMessageRight<MessageSendingType::MPLP>, 1, 2, variableMessageNumber, 1, FMC_PARAM, 1 > UnaryPairwiseMessageRightContainer;
   //typedef MessageContainer<LeftMargMessageMPLP, 1, 2, variableMessageNumber, 1, variableMessageSize, FMC_PARAM, 0 > UnaryPairwiseMessageLeft;
   //typedef MessageContainer<RightMargMessageMPLP, 1, 2, variableMessageNumber, 1, variableMessageSize, FMC_PARAM, 1 > UnaryPairwiseMessageRight;
   using UnaryToAssignmentMessageType = UnaryToAssignmentMessageCS2<McfCoveringFactor,MessagePassingType::HUNGARIAN>;
   typedef MessageContainer<UnaryToAssignmentMessageType, 1, 0, 1, variableMessageNumber, FMC_PARAM, 2> UnaryToAssignmentMessageContainer;

   // tightening
   typedef FactorContainer<SimpleTighteningTernarySimplexFactor, FMC_PARAM, 3 > EmptyTripletFactor;
   typedef MessageContainer<PairwiseTripletMessage12<MessageSendingType::SRMP>, 2, 3, variableMessageNumber, 1, FMC_PARAM, 3> PairwiseTriplet12MessageContainer;
   typedef MessageContainer<PairwiseTripletMessage13<MessageSendingType::SRMP>, 2, 3, variableMessageNumber, 1, FMC_PARAM, 4> PairwiseTriplet13MessageContainer;
   typedef MessageContainer<PairwiseTripletMessage23<MessageSendingType::SRMP>, 2, 3, variableMessageNumber, 1, FMC_PARAM, 5> PairwiseTriplet23MessageContainer;
   //typedef MessageContainer<PairwiseTriplet12MessageMPLP, 2, 3, variableMessageNumber, 1, variableMessageSize, FMC_PARAM, 3> PairwiseTriplet12MessageContainer;
   //typedef MessageContainer<PairwiseTriplet13MessageMPLP, 2, 3, variableMessageNumber, 1, variableMessageSize, FMC_PARAM, 4> PairwiseTriplet13MessageContainer;
   //typedef MessageContainer<PairwiseTriplet23MessageMPLP, 2, 3, variableMessageNumber, 1, variableMessageSize, FMC_PARAM, 5> PairwiseTriplet23MessageContainer;

   using FactorList = meta::list<MinCostFlowAssignmentFactor, UnaryFactor, PairwiseFactor, EmptyTripletFactor>;
   using MessageList = meta::list<
       UnaryPairwiseMessageLeftContainer,  
       UnaryPairwiseMessageRightContainer, 
       UnaryToAssignmentMessageContainer,
       PairwiseTriplet12MessageContainer, 
       PairwiseTriplet13MessageContainer, 
       PairwiseTriplet23MessageContainer 
      >;

   using mrf = AssignmentGmConstructor<StandardMrfConstructor<FMC_HUNGARIAN_BP_T,1,2,0,1>>;
   using tighteningMrf = TighteningMRFProblemConstructor<mrf,3,3,4,5>;
   using mrfLeft = tighteningMrf;
   using mrfRight = disable_write_constructor<tighteningMrf>;
   using ProblemDecompositionList = meta::list<mrfLeft,mrfRight>;
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
      template<typename INPUT>
      static void apply(const INPUT& in, GraphMatchingInput& gmInput)
      {
         gmInput.leftGraph_.resize(std::stoul(in.string()));
      }
   };
    
   template<> struct action< no_right_nodes > {
      template<typename INPUT>
      static void apply(const INPUT& in, GraphMatchingInput& gmInput)
      {
         gmInput.rightGraph_.resize(std::stoul(in.string()));
      }
   };
    
   template<> struct action< assignment > {
      template<typename INPUT>
      static void apply(const INPUT& in, GraphMatchingInput& gmInput)
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
      template<typename INPUT>
      static void apply(const INPUT & in, GraphMatchingInput& gmInput)
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

   std::vector<std::vector<REAL>> build_left_unaries(const GraphMatchingInput& gm, const REAL weight)
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

   std::vector<std::vector<REAL>> build_right_unaries(const GraphMatchingInput& gm, const REAL weight)
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
         std::cout << "This value is not useful, due to matching constraint: " << cost << "\n";
         cost = 1e10;
         return;
      }
      assert(graph[node1][oppositeNode1] != graph[node2][oppositeNode2]);
      //if(oppositeNode1 == oppositeNode2) {
      //   std::cout << "This value is not useful, due to matching constraint: " << cost << "\n";
      //   cost = 1e100;
      //   return;
      //}

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
                  cost[i1 + i2*leftDim] = std::numeric_limits<REAL>::infinity();
                  // transposed version
                  //assert(i1*rightDim + i2 < cost.size());
                  //cost[i1*rightDim + i2] = std::numeric_limits<REAL>::infinity();
               }
            }
         }
         q.insert(std::make_pair(std::make_pair(node1,node2), std::move(cost)));
      }

      // now add specific cost
      assert(q.find(std::make_pair(node1,node2)) != q.end());
      std::vector<REAL>& costVec = (*q.find(std::make_pair(node1,node2))).second;
      // do zrobienia: correct, or transpose?
      assert(costVec[index1 + index2*leftDim] == 0.0);
      costVec[index1 + index2*leftDim] = cost;
      // transposed version
      //assert(costVec[index1*rightDim + index2] == 0.0);
      //costVec[index1*rightDim + index2] = cost;
   }

   std::map<std::pair<INDEX,INDEX>, std::vector<REAL>> BuildLeftPairwisePotentials(const GraphMatchingInput& gmInput, const REAL weight)
   {
      std::map<std::pair<INDEX,INDEX>, std::vector<REAL>> q;
      for(auto i : gmInput.pairwise_potentials) {
         const INDEX leftNode1 = std::get<0>(i);
         const INDEX leftNode2 = std::get<1>(i);
         const INDEX rightNode1 = std::get<2>(i);
         const INDEX rightNode2 = std::get<3>(i);
         const REAL cost = std::get<4>(i);

         const INDEX rightIndex1 = std::find(gmInput.leftGraph_[leftNode1].begin(), gmInput.leftGraph_[leftNode1].end(), rightNode1) - gmInput.leftGraph_[leftNode1].begin();
         assert(rightIndex1 < gmInput.leftGraph_[leftNode1].size());
         const INDEX rightIndex2 = std::find(gmInput.leftGraph_[leftNode2].begin(), gmInput.leftGraph_[leftNode2].end(), rightNode2) - gmInput.leftGraph_[leftNode2].begin();
         assert(rightIndex2 < gmInput.leftGraph_[leftNode2].size());

         AddQuadraticPotential(q, leftNode1,leftNode2, rightIndex1,rightIndex2, weight*cost, gmInput.leftGraph_);
      }
      return q;
   }
   std::map<std::pair<INDEX,INDEX>, std::vector<REAL>> BuildRightPairwisePotentials(const GraphMatchingInput& gmInput, const REAL weight)
   {
      std::map<std::pair<INDEX,INDEX>, std::vector<REAL>> q;
      for(auto i : gmInput.pairwise_potentials) {
         const INDEX leftNode1 = std::get<0>(i);
         const INDEX leftNode2 = std::get<1>(i);
         const INDEX rightNode1 = std::get<2>(i);
         const INDEX rightNode2 = std::get<3>(i);
         const REAL cost = std::get<4>(i);

         const INDEX leftIndex1 = std::find(gmInput.rightGraph_[rightNode1].begin(), gmInput.rightGraph_[rightNode1].end(), leftNode1) - gmInput.rightGraph_[rightNode1].begin();
         assert(leftIndex1 < gmInput.rightGraph_[rightNode1].size());
         const INDEX leftIndex2 = std::find(gmInput.rightGraph_[rightNode2].begin(), gmInput.rightGraph_[rightNode2].end(), leftNode2) - gmInput.rightGraph_[rightNode2].begin();
         assert(leftIndex2 < gmInput.rightGraph_[rightNode2].size());

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
         std::map<std::pair<INDEX,INDEX>, std::vector<REAL>> leftQuadraticPot = BuildLeftPairwisePotentials(gmInput, pairwise_weight);
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
         std::map<std::pair<INDEX,INDEX>, std::vector<REAL>> rightQuadraticPot = BuildRightPairwisePotentials(gmInput, pairwise_weight);
         for(auto& q : rightQuadraticPot) {
            auto p = right_mrf.AddPairwiseFactor(q.first.first, q.first.second, q.second);
         }
      }
   }

   template<typename SOLVER>
   void construct_mp(SOLVER& s, GraphMatchingInput& gm_input)
   {
      using FMC = typename SOLVER::FMC;
      constexpr PairwiseConstruction pc = FmcConstruction(FMC{});

      auto& mrf_left = s.template GetProblemConstructor<0>();
      auto& mrf_right = s.template GetProblemConstructor<1>();
      /*

      mrf_left.SetGraph(gm_input.leftGraph_);
      mrf_right.SetGraph(gm_input.rightGraph_);

      if(pc == PairwiseConstruction::Left) {
         construct_left_mrf(gm_input, mrf_left, 1.0, 1.0);
         construct_right_mrf(gm_input, mrf_right, 0.0, 0.0);
      }
      if(pc == PairwiseConstruction::BothSides) {
         construct_left_mrf(gm_input, mrf_left, 0.5, 0.5);
         construct_right_mrf(gm_input, mrf_right, 0.5, 0.5);
      }
      if(pc == PairwiseConstruction::Right) {
         construct_left_mrf(gm_input, mrf_left, 0.0, 0.0);
         construct_right_mrf(gm_input, mrf_right, 1.0, 1.0);
      }
      */

      std::vector<INDEX> left_label_count(mrf_left.GetNumberOfVariables(),0);
      std::vector<INDEX> right_label_count(mrf_right.GetNumberOfVariables(),0);
      for(auto& a : gm_input.assignment_) {
         // get left and right unaries and connect them with a message
         auto *l = mrf_left.GetUnaryFactor(a.left_node_);
         auto *r = mrf_right.GetUnaryFactor(a.right_node_);
         const INDEX left_label = left_label_count[a.left_node_];
         const INDEX right_label = right_label_count[a.right_node_];
         auto* m = new typename FMC::AssignmentConstraintMessage( typename FMC::EqualityMessageType(left_label, right_label), l, r);
         s.GetLP().AddMessage(m);
         ++left_label_count[a.left_node_];
         ++right_label_count[a.right_node_];
      }
      for(INDEX i=0; i<mrf_left.GetNumberOfVariables(); ++i) {
         assert(mrf_left.GetNumberOfLabels(i) == left_label_count[i]+1);
      }
      for(INDEX i=0; i<mrf_right.GetNumberOfVariables(); ++i) {
         assert(mrf_right.GetNumberOfLabels(i) == right_label_count[i]+1);
      }
      
      mrf_left.Construct(s);
      mrf_right.Construct(s);
   }

   // add mcf factor, but assume graphical model has already been built.
   template<typename SOLVER>
   void construct_mcf(SOLVER& s, GraphMatchingInput& gm_input)
   {
      using FMC = typename SOLVER::FMC;
      // build assignment problem
      const INDEX no_left_nodes = std::accumulate(gm_input.assignment_.begin(), gm_input.assignment_.end(), 0, [](INDEX no, auto a) { return std::max(no, a.left_node_); }) + 1;
      const INDEX no_right_nodes = std::accumulate(gm_input.assignment_.begin(), gm_input.assignment_.end(), 0, [](INDEX no, auto a) { return std::max(no, a.right_node_); }) + 1;
      const INDEX no_nodes = no_left_nodes + no_right_nodes;
      const INDEX no_edges = gm_input.assignment_.size() + no_left_nodes + no_right_nodes + 1;

      std::vector<typename MinCostFlowFactorCS2::Edge> edges;
      edges.reserve(no_edges);
      for(auto a : gm_input.assignment_) {
         edges.push_back({a.left_node_, no_left_nodes + a.right_node_, 0, 1, 0.0});
      }
      std::vector<SIGNED_INDEX> demands;
      demands.reserve(no_left_nodes + no_right_nodes + 2);

      // slack edges
      for(INDEX i=0; i<no_left_nodes; ++i) {
         edges.push_back({i,no_nodes + 1, 0, 1, 0.0}); // for non-assignment
         demands.push_back(1);
      }
      for(INDEX i=0; i<no_right_nodes; ++i) {
         edges.push_back({no_nodes, no_left_nodes + i, 0, 1, 0.0}); // for non-assignment
         demands.push_back(-1);
      }
      edges.push_back({no_nodes, no_nodes + 1, 0, std::max(no_left_nodes, no_right_nodes), 0.0});
      demands.push_back(no_right_nodes);
      demands.push_back(-no_left_nodes);

      const INDEX no_covered_edges = edges.size() - 1;
      auto* f = new typename FMC::MinCostFlowAssignmentFactor( MinCostFlowFactorCS2(edges, demands, no_covered_edges) );
      auto* mcf = f->GetFactor()->GetMinCostFlowSolver();
      s.GetLP().AddFactor(f);

      // connect assignment factor with unaries
      // left side
      {
      auto& mrf_left = s.template GetProblemConstructor<0>();
      std::vector<std::vector<INDEX>> edgeId(no_left_nodes);
      INDEX c=0;
      for(auto& a : gm_input.assignment_) {
         // get left and right unaries and connect them with a message
         edgeId[a.left_node_].push_back(c);
         ++c;
      }
      // add slack edges
      for(INDEX i=0; i<no_left_nodes; ++i) {
         edgeId[i].push_back(c);
         ++c;
      }
      for(INDEX i=0; i<no_left_nodes; ++i) {
         //assert(mrf_left.GetNumberOfLabels(i) == mcf->NoArcs(i));
         auto *u = mrf_left.GetUnaryFactor(i);
         using MessageType = typename FMC::UnaryToAssignmentMessageType;
         //auto *m = new typename FMC::UnaryToAssignmentMessageContainer( MessageType(mcf->StartingArc(i), mcf->NoArcs(i)), u, f, mrf_left.GetNumberOfLabels(i));
         auto *m = new typename FMC::UnaryToAssignmentMessageContainer( MessageType(edgeId[i]), u, f);
         s.GetLP().AddMessage(m);
      }
      }
      // right side
      {
      auto& mrf_right = s.template GetProblemConstructor<1>();
      std::vector<std::vector<INDEX>> edgeId(no_right_nodes);
      INDEX c=0;
      for(auto& a : gm_input.assignment_) {
         // get left and right unaries and connect them with a message
         edgeId[a.right_node_].push_back(c);
         ++c;
      }
      // add slack edges
      c += no_left_nodes;
      for(INDEX i=0; i<no_right_nodes; ++i) {
         edgeId[i].push_back(c);
         ++c;
      }
      for(INDEX i=0; i<no_right_nodes; ++i) {
         //assert(mrf_right.GetNumberOfLabels(i) == mcf->NoArcs(no_left_nodes + i));
         auto *u = mrf_right.GetUnaryFactor(i);
         using MessageType = typename FMC::UnaryToAssignmentMessageType;
         //auto *m = new typename FMC::UnaryToAssignmentMessageContainer( MessageType(mcf->StartingArc(i + no_left_nodes), mcf->NoArcs(i + no_left_nodes)), u, f, mrf_right.GetNumberOfLabels(i));
         auto *m = new typename FMC::UnaryToAssignmentMessageContainer( MessageType(edgeId[i]), u, f);
         s.GetLP().AddMessage(m);
      }
      }
   }

   /*
   template<typename FMC>
   void construct_mcf(SOLVER& s, GraphMatchingInput& gm_input)
   {
      constexpr PairwiseConstruction pc = FmcConstruction(FMC{});

      auto& mrf_left = s.template GetProblemConstructor<0>();
      if(pc == PairwiseConstruction::Left || pc == PairwiseConstruction::BothSides) {
         mrf_left.SetGraph(gm_input.leftGraph_);
      }
      auto& mrf_right = s.template GetProblemConstructor<1>();
      if(pc == PairwiseConstruction::Right || pc == PairwiseConstruction::BothSides) {
         mrf_right.SetGraph(gm_input.rightGraph_);
      }
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
      const INDEX no_nodes = no_left_nodes + no_right_nodes;
      const INDEX no_edges = gm_input.assignment_.size() + no_left_nodes + no_right_nodes + 1;

      std::vector<typename MinCostFlowFactorCS2::Edge> edges;
      edges.reserve(no_edges);
      for(auto a : gm_input.assignment_) {
         edges.push_back({a.left_node_, no_left_nodes + a.right_node_, 0, 1, 0.0});
      }
      std::vector<SIGNED_INDEX> demands;
      demands.reserve(no_left_nodes + no_right_nodes + 2);

      // slack edges
      for(INDEX i=0; i<no_left_nodes; ++i) {
         edges.push_back({i,no_nodes + 1, 0, 1, 0.0}); // for non-assignment
         demands.push_back(1);
      }
      for(INDEX i=0; i<no_right_nodes; ++i) {
         edges.push_back({no_nodes, no_left_nodes + i, 0, 1, 0.0}); // for non-assignment
         demands.push_back(-1);
      }
      edges.push_back({no_nodes, no_nodes + 1, 0, std::max(no_left_nodes, no_right_nodes), 0.0});
      demands.push_back(no_right_nodes);
      demands.push_back(-no_left_nodes);

      // check for duplicate edges
      //std::sort(edges.begin(), edges.end(), [](const auto& a, const auto& b) { 
      //      if(a.start_node != b.start_node) {
      //      return a.start_node < b.start_node;
      //      } else {
      //      return a.end_node < b.end_node;
      //      }
      //      });
      //for(INDEX a=1; a<edges.size(); ++a) {
      //   assert(!((edges[a].start_node == edges[a-1].start_node) && (edges[a].end_node == edges[a-1].end_node)));
      //}
      

      const INDEX no_covered_edges = edges.size() - no_right_nodes - 1;
      auto* f = new typename FMC::MinCostFlowAssignmentFactor( MinCostFlowFactorCS2(edges, demands, no_covered_edges) );
      auto* mcf = f->GetFactor()->GetMinCostFlowSolver();
      s.GetLP().AddFactor(f);

      // connect assignment factor with unaries
      if(pc == PairwiseConstruction::Left || pc == PairwiseConstruction::BothSides) {
         std::vector<std::vector<INDEX>> edgeId(no_left_nodes);
         INDEX c=0;
         for(auto& a : gm_input.assignment_) {
            // get left and right unaries and connect them with a message
            edgeId[a.left_node_].push_back(c);
            ++c;
         }
         // add slack edges
         for(INDEX i=0; i<no_left_nodes; ++i) {
            edgeId[i].push_back(c);
            ++c;
         }
         for(INDEX i=0; i<no_left_nodes; ++i) {
            //assert(mrf_left.GetNumberOfLabels(i) == mcf->NoArcs(i));
            auto *u = mrf_left.GetUnaryFactor(i);
            using MessageType = typename FMC::UnaryToAssignmentMessageType;
            //auto *m = new typename FMC::UnaryToAssignmentMessageContainer( MessageType(mcf->StartingArc(i), mcf->NoArcs(i)), u, f, mrf_left.GetNumberOfLabels(i));
            auto *m = new typename FMC::UnaryToAssignmentMessageContainer( MessageType(edgeId[i]), u, f, mrf_left.GetNumberOfLabels(i));
            s.GetLP().AddMessage(m);
         }
      }
      if(pc == PairwiseConstruction::Right || pc == PairwiseConstruction::BothSides) {
         std::vector<std::vector<INDEX>> edgeId(no_right_nodes);
         INDEX c=0;
         for(auto& a : gm_input.assignment_) {
            // get left and right unaries and connect them with a message
            edgeId[a.right_node_].push_back(c);
            ++c;
         }
         // add slack edges
         c += no_left_nodes;
         for(INDEX i=0; i<no_right_nodes; ++i) {
            edgeId[i].push_back(c);
            ++c;
         }
         for(INDEX i=0; i<no_right_nodes; ++i) {
            //assert(mrf_right.GetNumberOfLabels(i) == mcf->NoArcs(no_left_nodes + i));
            auto *u = mrf_right.GetUnaryFactor(i);
            using MessageType = typename FMC::UnaryToAssignmentMessageType;
            //auto *m = new typename FMC::UnaryToAssignmentMessageContainer( MessageType(mcf->StartingArc(i + no_left_nodes), mcf->NoArcs(i + no_left_nodes)), u, f, mrf_right.GetNumberOfLabels(i));
            auto *m = new typename FMC::UnaryToAssignmentMessageContainer( MessageType(edgeId[i]), u, f, mrf_right.GetNumberOfLabels(i));
            s.GetLP().AddMessage(m);
         }
      }
   }
*/

   template<typename SOLVER>
   void construct_gm(SOLVER& s, GraphMatchingInput& gm_input)
   {
      using FMC = typename SOLVER::FMC;
      constexpr PairwiseConstruction pc = FmcConstruction(FMC{});

      auto& mrf_left = s.template GetProblemConstructor<0>();
      auto& mrf_right = s.template GetProblemConstructor<1>();

      mrf_left.SetGraph(gm_input.leftGraph_);
      mrf_right.SetGraph(gm_input.rightGraph_);

      if(pc == PairwiseConstruction::Left) {
         construct_left_mrf(gm_input, mrf_left, 1.0, 1.0);
         construct_right_mrf(gm_input, mrf_right, 0.0, 0.0);
      }
      if(pc == PairwiseConstruction::BothSides) {
         construct_left_mrf(gm_input, mrf_left, 0.5, 0.5);
         construct_right_mrf(gm_input, mrf_right, 0.5, 0.5);
      }
      if(pc == PairwiseConstruction::Right) {
         construct_left_mrf(gm_input, mrf_left, 0.0, 0.0);
         construct_right_mrf(gm_input, mrf_right, 1.0, 1.0);
      }

      // note: possibly more pairwise potentials have to be added for some problems to ensure uniqueness constraint in matching. For house and hotel datasets this is not needed
      /*
      auto& mrf = s.template GetProblemConstructor<0>();
      constexpr PairwiseConstruction pc = FmcConstruction(FMC{});
      if(pc == PairwiseConstruction::Left) {
         mrf.SetGraph(gm_input.leftGraph_);
         construct_left_mrf(gm_input, mrf, 1.0, 1.0);
      }
      if(pc == PairwiseConstruction::Right) {
         mrf.SetGraph(gm_input.rightGraph_);
         construct_right_mrf(gm_input, mrf, 1.0, 1.0);
      }
      */
   }

   
   GraphMatchingInput ParseFile(const std::string& filename)
   {
      GraphMatchingInput gmInput;
      pegtl::file_parser problem(filename);
      std::cout << "parsing " << filename << "\n";

      for(auto& edges : gmInput.leftGraph_) {
              assert(std::is_sorted(edges.begin(), edges.end()));
      }
      for(auto& edges : gmInput.rightGraph_) {
              assert(std::is_sorted(edges.begin(), edges.end()));
      }
      const bool ret = problem.parse< grammar, action >( gmInput );
      if(!ret) {
         throw std::runtime_error("could not read file " + filename);
      }
      return std::move(gmInput);
   }

   template<typename SOLVER>
   bool ParseProblemGM(const std::string& filename, SOLVER& s)
   {
      auto input = ParseFile(filename);
      construct_gm( s, input );
      return true;
   }

   template<typename SOLVER>
   bool ParseProblemMP(const std::string& filename, SOLVER& s)
   {
      auto input = ParseFile(filename);
      construct_gm( s, input );
      construct_mp( s, input );
      return true;
   }

   template<typename SOLVER>
   bool ParseProblemMCF(const std::string& filename, SOLVER& s)
   {
      auto input = ParseFile(filename);
      construct_gm( s, input );
      construct_mp( s, input );
      construct_mcf( s, input );
      return true;
   }

   template<typename SOLVER>
   bool ParseProblemHungarian(const std::string& filename, SOLVER& s)
   {
      auto input = ParseFile(filename);
      construct_gm( s, input );
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
      template<typename INPUT>
      static void apply(const INPUT & in, matching & m)
      { 
         const INDEX var = std::stoul(in.string());
         assert(m.size() == var);
         m.push_back(std::vector<INDEX>(0));
      }
   };

   template<> struct action< label > {
      template<typename INPUT>
      static void apply(const INPUT & in, matching & m)
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

   template<typename SOLVER>
   bool ParseProblemGM(const std::string& filename, SOLVER& s)
   {
      using FMC = typename SOLVER::FMC;
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

   template<typename SOLVER>
   void construct_mp(SOLVER& s, const input& i)
   {
      using FMC = typename SOLVER::FMC;
      auto& mrf_left = s.template GetProblemConstructor<0>();
      auto& mrf_input = std::get<0>(i);
      UaiMrfInput::build_mrf(mrf_left, mrf_input);

      // now build unaries for mrf_right. There will be as many unaries (=labels) on the right as constraints
      auto& mrf_right = s.template GetProblemConstructor<1>();
      auto& matching = std::get<1>(i);
      auto constraints = invert_matching(matching);
      for(auto& c : constraints) {
         auto* u_r = mrf_right.AddUnaryFactor(std::vector<REAL>(c.size()+1,0.0)); // extra label is for non-assignment of label
         for(INDEX var_label=0; var_label<c.size(); ++var_label) {
            const INDEX var = std::get<0>(c[var_label]);
            const INDEX label = std::get<1>(c[var_label]);
            auto* u_l = mrf_left.GetUnaryFactor(var);
            auto* m = new typename FMC::AssignmentConstraintMessage( typename FMC::EqualityMessageType(label, var_label), u_l, u_r);
            s.GetLP().AddMessage(m);
         }
      }
      
      std::cout << "Constructed gm with " << mrf_left.GetNumberOfVariables() << " unary factors and " << mrf_left.GetNumberOfPairwiseFactors() << " pairwise factors\n";
   }

   template<typename SOLVER>
   void construct_mcf(SOLVER& s, const input& i)
   {
      using FMC = typename SOLVER::FMC;
      auto& mrf_left = s.template GetProblemConstructor<0>();
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

      std::vector<std::vector<INDEX>> edgeId(no_left_nodes);
      for(INDEX i=0; i<edges.size(); ++i) {
              const INDEX left_node = edges[i].start_node;
              if(left_node < no_left_nodes) {
                      edgeId[left_node].push_back(i);
              }
      }

      std::vector<SIGNED_INDEX> demands(no_left_nodes + total_no_right_nodes + 1);
      std::fill(demands.begin(), demands.begin() + no_left_nodes, 1);
      std::fill(demands.begin() + no_left_nodes, demands.end(), 0);
      demands.back() = -no_left_nodes;

      auto* f = new typename FMC::MinCostFlowAssignmentFactor( MinCostFlowFactorCS2(edges, demands, edges.size()) );
      s.GetLP().AddFactor(f);
      auto* mcf = f->GetFactor()->GetMinCostFlowSolver();

      for(INDEX i=0; i<no_left_nodes; ++i) {
         auto *u = mrf_left.GetUnaryFactor(i);
         using MessageType = typename FMC::UnaryToAssignmentMessageType;
         auto *m = new typename FMC::UnaryToAssignmentMessageContainer( MessageType(edgeId[i]), u, f, mrf_left.GetNumberOfLabels(i));
         s.GetLP().AddMessage(m);
      }
   }

   template<typename SOLVER>
   bool ParseProblemMP(const std::string& filename, SOLVER& s)
   {
      // do zrobienia: FMC must be left type -> static_assert this
      using FMC = typename SOLVER::FMC;
      static_assert(FmcConstruction(FMC{}) == PairwiseConstruction::Left, "in uai format only left construction makes sense"); 
      const auto input = ParseFile(filename);
      construct_mp(s, input);
      return true;
   }

   template<typename SOLVER>
   bool ParseProblemMCF(const std::string& filename, SOLVER& s)
   {
      using FMC = typename SOLVER::FMC;
      static_assert(FmcConstruction(FMC{}) == PairwiseConstruction::Left, "in uai format only left construction makes sense"); 
      const auto input = ParseFile(filename);
      construct_mcf(s, input);
      return true;
   }
}



#endif // LP_MP_GRAPH_MATCHING_H

