#ifndef LP_MP_COSEGMENTATION_BINS_HXX
#define LP_MP_COSEGMENTATION_BINS_HXX

#include "factors_messages.hxx"
#include "LP_MP.h"
#include "factors/multiplex_factor.hxx"
#include "const_array_types.h"
#include "messages/multiplex_marg_message.hxx"
#include "messages/equality_message.hxx"
#include "marginal_summation_message.hxx"
#include "message_replicator.hxx"
#include "message_replicator_factor.hxx"

#include "assignment_problem_construction.hxx"
#include "potts_problem_construction.hxx"
#include "marginal_summation_replicator.hxx"

namespace LP_MP {

// support Wasserstein distance with histogram bins into which more than one pixel belongs -> big images segmented on pixel basis
class CosegmentationBins {

public:
   typedef MultiplexFactor<std::vector<REAL>, const_ones_array, const_one> Simplex;

   typedef UnaryLoop UnaryLoopType;
   typedef PairwiseLoop<0> LeftLoopType;
   typedef PairwiseLoop<1> RightLoopType;

   typedef MultiplexMargMessage<UnaryLoopType,LeftLoopType,true,false> LeftMargMessage;
   typedef MultiplexMargMessage<UnaryLoopType,RightLoopType,true,false> RightMargMessage;

   // factor-message-network
   struct FMC; // forward declaration
   // factors
   typedef FactorContainer<Simplex, ExplicitRepamStorage, FMC, 0 > UnaryFactor; // 2 entries: foreground, background. do zrobienia: specialize factor
   typedef FactorContainer<Simplex, ExplicitRepamStorage, FMC, 1 > PairwiseFactor; // 4 entries for pairwise potential. do zrobienia: specialize factor
   typedef FactorContainer<MultiplexFactor<std::vector<REAL>, constant_array, INDEX>, ExplicitRepamStorage, FMC, 2 > AssignmentFactor;
   typedef FactorContainer<MessageReplicatorFactor, ImplicitRepamStorage, FMC, 3> SummationReplicatorFactor; // for bin size > 1

   //messages
   typedef MessageContainer<EqualityMessage, FixedMessageStorage<1>, FMC, 0 > AssignmentConstraintMessage;
   typedef MessageContainer<LeftMargMessage, MessageStorageSIMD, FMC, 1 > UnaryPairwiseMessageLeft; // unary has left_dim entries
   typedef MessageContainer<RightMargMessage, MessageStorageSIMD, FMC, 2 > UnaryPairwiseMessageRight; // unary has right_dim entries
   typedef MessageContainer<MessageReplicator<Chirality::left, MarginalSummationMessage>, FixedMessageStorage<2>, FMC, 3> MarginalSummationMessageContainerLeft; // for bin size > 1
   typedef MessageContainer<MessageReplicator<Chirality::right, MarginalSummationMessage>, FixedMessageStorage<2>, FMC, 4> MarginalSummationMessageContainerRight; // for bin size > 1

   struct FMC {
      using factor_list = meta::list< UnaryFactor, PairwiseFactor, AssignmentFactor, SummationReplicatorFactor>;
      using msg_list = meta::list< 
                                   meta::list< AssignmentConstraintMessage, meta::size_t<2>, meta::size_t<2>, std::vector<AssignmentConstraintMessage*>, std::vector<AssignmentConstraintMessage*> >,
                                   meta::list< UnaryPairwiseMessageLeft,  meta::size_t<0>, meta::size_t<1>, std::vector<UnaryPairwiseMessageLeft*>, FixedSizeContainer<UnaryPairwiseMessageLeft*,1> >,
                                   meta::list< UnaryPairwiseMessageRight, meta::size_t<0>, meta::size_t<1>, std::vector<UnaryPairwiseMessageRight*>, FixedSizeContainer<UnaryPairwiseMessageRight*,1> >,
                                   meta::list< MarginalSummationMessageContainerLeft, meta::size_t<0>, meta::size_t<3>, std::vector<MarginalSummationMessageContainerLeft*>, std::vector<MarginalSummationMessageContainerLeft*> >,
                                   meta::list< MarginalSummationMessageContainerRight, meta::size_t<3>, meta::size_t<2>, std::vector<MarginalSummationMessageContainerRight*>, std::vector<MarginalSummationMessageContainerRight*> >
                                      >;

      using lap = AssignmentProblemConstructor<FMC,2,0>;
      using potts1 = PottsProblemConstructor<FMC,0,1,1,2>;
      using potts2 = PottsProblemConstructor<FMC,0,1,1,2>;
      using marg_consistency1 = MarginalSummationReplicator<FMC, 3, 3, 4, 0, 1>;
      using marg_consistency2 = MarginalSummationReplicator<FMC, 3, 3, 4, 0, 1>;
      using problem_decomposition = meta::list< lap, potts1, potts2, marg_consistency1, marg_consistency2 >;
   };
};

} // end namespace LP_MP

#endif // LP_MP_COSEGMENTATION_BINS_HXX

