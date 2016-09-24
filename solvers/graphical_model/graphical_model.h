#ifndef LP_MP_GRAPH_MATCHING_H
#define LP_MP_GRAPH_MATCHING_H

//#include "problem_decomposition.hxx"
#include "factors_messages.hxx"
#include "LP_MP.h"
#include "factors/simplex_factor.hxx"
//#include "messages/unary_pairwise_mcf_message.hxx"
#include "const_array_types.h"
#include "messages/simplex_marginalization_message.hxx"
#include "problem_constructors/mrf_problem_construction.hxx"

#include "parse_rules.h"
#include "hdf5_routines.hxx"

// this file contains solvers for graphical models, among them SRMP and MPLP

using namespace LP_MP;

typedef UnaryLoop<> UnaryLoopType;
// do zrobienia: use enum chirality type in Pairwise Loop specification
typedef PairwiseLoop<0> LeftLoopType;
typedef PairwiseLoop<1> RightLoopType;

typedef SimplexFactor<> Simplex;

typedef PairwiseTripletLoop<0,1> PairwiseTripletLoopType12;
typedef PairwiseTripletLoop<0,2> PairwiseTripletLoopType13;
typedef PairwiseTripletLoop<1,2> PairwiseTripletLoopType23;
typedef SimplexMarginalizationMessage<UnaryLoopType,PairwiseTripletLoopType12,true,false,false,true> PairwiseTriplet12Message;
typedef SimplexMarginalizationMessage<UnaryLoopType,PairwiseTripletLoopType13,true,false,false,true> PairwiseTriplet13Message;
typedef SimplexMarginalizationMessage<UnaryLoopType,PairwiseTripletLoopType23,true,false,false,true> PairwiseTriplet23Message;



struct FMC_SRMP { // equivalent to SRMP or TRWS
   constexpr static const char* name = "SRMP for pairwise case = TRWS";

   typedef FactorContainer<Simplex, ExplicitRepamStorage, FMC_SRMP, 0, true, true > UnaryFactor;
   typedef FactorContainer<Simplex, ExplicitRepamStorage, FMC_SRMP, 1, false, false > PairwiseFactor;

   typedef SimplexMarginalizationMessage<UnaryLoopType,LeftLoopType,true,false,false,true> LeftMargMessage;
   typedef SimplexMarginalizationMessage<UnaryLoopType,RightLoopType,true,false,false,true> RightMargMessage;

   typedef MessageContainer<LeftMargMessage, 0, 1, variableMessageNumber, 1, variableMessageSize, FMC_SRMP, 0 > UnaryPairwiseMessageLeft;
   typedef MessageContainer<RightMargMessage, 0, 1, variableMessageNumber, 1, variableMessageSize, FMC_SRMP, 1 > UnaryPairwiseMessageRight;

   using FactorList = meta::list< UnaryFactor, PairwiseFactor >;
   using MessageList = meta::list< UnaryPairwiseMessageLeft, UnaryPairwiseMessageRight >;

   using mrf = StandardMrfConstructor<FMC_SRMP,0,1,0,1>;
   using ProblemDecompositionList = meta::list<mrf>;
};

struct FMC_SRMP_T { // equivalent to SRMP or TRWS
   constexpr static const char* name = "SRMP for pairwise case with tightening triplets";

   typedef FactorContainer<Simplex, ExplicitRepamStorage, FMC_SRMP_T, 0, true, true > UnaryFactor;
   typedef FactorContainer<Simplex, ExplicitRepamStorage, FMC_SRMP_T, 1, false, false > PairwiseFactor;

   typedef SimplexMarginalizationMessage<UnaryLoopType,LeftLoopType,true,false,false,true> LeftMargMessage;
   typedef SimplexMarginalizationMessage<UnaryLoopType,RightLoopType,true,false,false,true> RightMargMessage;

   typedef MessageContainer<LeftMargMessage, 0, 1, variableMessageNumber, 1, variableMessageSize, FMC_SRMP_T, 0 > UnaryPairwiseMessageLeft;
   typedef MessageContainer<RightMargMessage, 0, 1, variableMessageNumber, 1, variableMessageSize, FMC_SRMP_T, 1 > UnaryPairwiseMessageRight;

   // tightening
   typedef FactorContainer<Simplex, ExplicitRepamStorage, FMC_SRMP_T, 2, false, false > EmptyTripletFactor;
   typedef MessageContainer<PairwiseTriplet12Message, 1, 2, variableMessageNumber, 1, variableMessageSize, FMC_SRMP_T, 2> PairwiseTriplet12MessageContainer;
   typedef MessageContainer<PairwiseTriplet13Message, 1, 2, variableMessageNumber, 1, variableMessageSize, FMC_SRMP_T, 3> PairwiseTriplet13MessageContainer;
   typedef MessageContainer<PairwiseTriplet23Message, 1, 2, variableMessageNumber, 1, variableMessageSize, FMC_SRMP_T, 4> PairwiseTriplet23MessageContainer;

   using FactorList = meta::list< UnaryFactor, PairwiseFactor, EmptyTripletFactor >;
   using MessageList = meta::list< UnaryPairwiseMessageLeft, UnaryPairwiseMessageRight,
         PairwiseTriplet12MessageContainer, PairwiseTriplet13MessageContainer, PairwiseTriplet23MessageContainer
         >;

   using mrf = StandardMrfConstructor<FMC_SRMP_T,0,1,0,1>;
   using tighteningMrf = TighteningMRFProblemConstructor<mrf,2,2,3,4>;
   using ProblemDecompositionList = meta::list<tighteningMrf>;
};



struct FMC_MPLP {
   constexpr static const char* name = "MPLP for pairwise case";

   typedef FactorContainer<Simplex, ExplicitRepamStorage, FMC_MPLP, 0, true, true > UnaryFactor;
   typedef FactorContainer<Simplex, ExplicitRepamStorage, FMC_MPLP, 1, false, false > PairwiseFactor;

   typedef SimplexMarginalizationMessage<UnaryLoopType,LeftLoopType,false,true,false,true> LeftMargMessage;
   typedef SimplexMarginalizationMessage<UnaryLoopType,RightLoopType,false,true,false,true> RightMargMessage;

   typedef MessageContainer<LeftMargMessage, 0, 1, variableMessageNumber, 1, variableMessageSize, FMC_MPLP, 0 > UnaryPairwiseMessageLeft;
   typedef MessageContainer<RightMargMessage, 0, 1, variableMessageNumber, 1, variableMessageSize, FMC_MPLP, 1 > UnaryPairwiseMessageRight;

   using FactorList = meta::list< UnaryFactor, PairwiseFactor >;
   using MessageList = meta::list< UnaryPairwiseMessageLeft, UnaryPairwiseMessageRight >;

   using mrf = StandardMrfConstructor<FMC_MPLP,0,1,0,1>;
   using ProblemDecompositionList = meta::list<mrf>;
};

/*
struct FMC_MPLP_MCF {
   constexpr static const char* name = "MPLP with minimum cost flow message updates";

   typedef FactorContainer<Simplex, ExplicitRepamStorage, FMC_MPLP_MCF, 0, true, true > UnaryFactor;
   typedef FactorContainer<PairwiseMcfFactor, ExplicitRepamStorage, FMC_MPLP_MCF, 1, false, false > PairwiseFactor;

   typedef MessageContainer<UnaryPairwiseMcfMessage, 0, 1, variableMessageNumber, 2, variableMessageSize, FMC_MPLP_MCF, 0 > UnaryPairwiseMessage;

   using FactorList = meta::list< UnaryFactor, PairwiseFactor >;
   using MessageList = meta::list<UnaryPairwiseMessage>;

   struct MrfConstruction : public MRFProblemConstructor<FMC_MPLP_MCF,0,1,0,0> {
      using MRFProblemConstructor<FMC_MPLP_MCF,0,1,0,0>::MRFProblemConstructor;
      UnaryFactorType ConstructUnaryFactor(const std::vector<REAL>& cost) 
      { return UnaryFactorType(cost); }
      PairwiseFactorType ConstructPairwiseFactor(const std::vector<REAL>& cost, const INDEX leftDim, const INDEX rightDim) 
      { return PairwiseFactorType(leftDim, rightDim); }
      RightMessageType ConstructRightUnaryPairwiseMessage(UnaryFactorContainer* const right, PairwiseFactorContainer* const p)
      { return RightMessageType(right->size(), Chirality::right); }
      LeftMessageType ConstructLeftUnaryPairwiseMessage(UnaryFactorContainer* const left, PairwiseFactorContainer* const p)
      { return LeftMessageType(left->size(), Chirality::left); }
   };
   using mrf = MrfConstruction;
   using ProblemDecompositionList = meta::list<mrf>;
};
*/



#endif // LP_MP_GRAPH_MATCHING_H

