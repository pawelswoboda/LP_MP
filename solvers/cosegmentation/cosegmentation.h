#ifndef LP_MP_COSEGMENTATION_HXX
#define LP_MP_COSEGMENTATION_HXX

#include "factors_messages.hxx"
#include "LP_MP.h"
#include "factors/multiplex_factor.hxx"
#include "factors/simplex_factor_simd.hxx"
#include "const_array_types.h"
#include "messages/multiplex_marg_message.hxx"
#include "messages/equality_message.hxx"
#include "marginal_summation_message.hxx"
#include "problem_construction_helper.hxx"

#include "message_replicator.hxx"
#include "message_replicator_factor.hxx"

#include "MinCost/MinCost.h"

#include <iostream>
#include <fstream>
#include <string>
#include <tuple>
#include <exception>
#include <algorithm>
#include <memory>

namespace LP_MP {

// similar to quadratic assignment problem, except that we use modified Potts factors and no other quadratic potentials
// many functions are duplicates from graph_matching
// do zrobienia: implement cosegmentation with bins

class Cosegmentation {
public:
   Cosegmentation() : lp_(nullptr) {}
   ~Cosegmentation() { delete lp_; }
   
   Cosegmentation& SetNumberLeftNodes(const INDEX nodeNumber);
   Cosegmentation& SetNumberRightNodes(const INDEX nodeNumber);
   Cosegmentation& AddAssignmentCost(const INDEX leftNode, const INDEX rightNode, const REAL cost);
   Cosegmentation& AddLeftPottsTerm(const INDEX leftNode1, const INDEX leftNode2, const REAL cost);
   Cosegmentation& AddRightPottsTerm(const INDEX rightNode1, const INDEX rightNode2, const REAL cost);
   Cosegmentation& Solve(const INDEX nIter);
   std::vector<int> GetLeftSegmentation();
   std::vector<int> GetRightSegmentation();
   REAL dualBound() {return lp_->LowerBound();};
   REAL evalLeft(const std::vector<int>& seg);
   REAL evalRight(const std::vector<int>& seg);
   REAL evalAssignment(const std::vector<int>& assignment);
   REAL primalBound() { return 0.0; };

private:
   //specification for underlying factor/message-network
   // do zrobienia: hardcode 2 labels for more efficient message passing
   typedef UnaryLoop UnaryLoopType;
   typedef PairwiseLoop<0> LeftLoopType;
   typedef PairwiseLoop<1> RightLoopType;

   typedef MultiplexMargMessage<UnaryLoopType,LeftLoopType,true,false> LeftMargMessage;
   typedef MultiplexMargMessage<UnaryLoopType,RightLoopType,true,false> RightMargMessage;

   struct FMC; // forward declaration

   // factors
   typedef MultiplexFactor<std::vector<REAL>, const_ones_array, const_one> Simplex;
   typedef FactorContainer<Simplex, ExplicitRepamStorage, FMC, 0 > UnaryFactor; // 2 entries: foreground, background. do zrobienia: specialize factor
   typedef FactorContainer<Simplex, ExplicitRepamStorage, FMC, 1 > PairwiseFactor; // 4 entries for pairwise potential. do zrobienia: specialize factor
   typedef FactorContainer<Simplex, ExplicitRepamStorage, FMC, 2 > AssignmentFactor; // for each pixel: number of possible assignments
   typedef FactorContainer<MessageReplicatorFactor, ImplicitRepamStorage, FMC, 3> SummationReplicatorFactor; // do zrobienia: not supported yet

   // messages
   typedef MessageContainer<EqualityMessage, FixedMessageStorage<1>, FMC, 0 > AssignmentConstraintMessage;
   typedef MessageContainer<LeftMargMessage, MessageStorageSIMD, FMC, 1 > UnaryPairwiseMessageLeft; // unary has left_dim entries
   typedef MessageContainer<RightMargMessage, MessageStorageSIMD, FMC, 2 > UnaryPairwiseMessageRight; // unary has right_dim entries
   typedef MessageContainer<MarginalSummationMessage, FixedMessageStorage<2>, FMC, 3> MarginalSummationMessageContainer;
   typedef MessageContainer<MessageReplicator<Chirality::left, MarginalSummationMessage>, FixedMessageStorage<2>, FMC, 4> MarginalSummationMessageContainerLeft; // do zrobienia: not supported yet
   typedef MessageContainer<MessageReplicator<Chirality::right, MarginalSummationMessage>, FixedMessageStorage<2>, FMC, 5> MarginalSummationMessageContainerRight; // do zrobienia: not supported yet

   struct FMC {
      using factor_list = meta::list< UnaryFactor, PairwiseFactor, AssignmentFactor, SummationReplicatorFactor>;
      using msg_list = meta::list< 
                                   meta::list< AssignmentConstraintMessage, meta::size_t<2>, meta::size_t<2>, std::vector<AssignmentConstraintMessage*>, std::vector<AssignmentConstraintMessage*> >,
                                   meta::list< UnaryPairwiseMessageLeft,  meta::size_t<0>, meta::size_t<1>, std::vector<UnaryPairwiseMessageLeft*>, FixedSizeContainer<UnaryPairwiseMessageLeft*,1> >,
                                   meta::list< UnaryPairwiseMessageRight, meta::size_t<0>, meta::size_t<1>, std::vector<UnaryPairwiseMessageRight*>, FixedSizeContainer<UnaryPairwiseMessageRight*,1> >,
                                   meta::list< MarginalSummationMessageContainer, meta::size_t<0>, meta::size_t<2>, FixedSizeContainer<MarginalSummationMessageContainer*,1>, FixedSizeContainer<MarginalSummationMessageContainer*,1> >,
                                   meta::list< MarginalSummationMessageContainerLeft, meta::size_t<0>, meta::size_t<3>, std::vector<MarginalSummationMessageContainerLeft*>, std::vector<MarginalSummationMessageContainerLeft*> >,
                                   meta::list< MarginalSummationMessageContainerRight, meta::size_t<3>, meta::size_t<2>, std::vector<MarginalSummationMessageContainerRight*>, std::vector<MarginalSummationMessageContainerRight*> >
                                      >;
   };

   Cosegmentation& SetNumberNodes(const INDEX nodeNumber, std::vector<std::vector<INDEX> >& graph);
   void ConstructGraph(
         const std::map<std::pair<INDEX,INDEX>, REAL>& PottsCost,
         std::vector<UnaryFactor*>& unaryFactor, 
         std::vector<PairwiseFactor*>& pairwiseFactor);
   void ConstructAssignmentFactors(
         const std::vector<std::vector<REAL> >& cost,
         std::vector<AssignmentFactor*>& assignmentFactor);
   void ConstructSummationConstraints(
         const std::vector<UnaryFactor*>& unaryFactor,
         const std::vector<AssignmentFactor*>& assignmentFactor);
   void ConstructAssignmentConstraints(
         const std::vector<std::vector<INDEX> >& graph1,
         const std::vector<AssignmentFactor*>& assignmentFactor1,
         const std::vector<std::vector<INDEX> >& graph2,
         const std::vector<AssignmentFactor*>& assignmentFactor2);

   REAL eval(
         const std::vector<int>& assignment,
         const std::vector<std::vector<INDEX> >& graph, 
      const std::vector<std::vector<REAL> >& unaryCost,
      const std::map<std::pair<INDEX,INDEX>, REAL>& pottsCost);

   LP* lp_;
   // data structures for the underlying assignment problem
   // {left|right}-graph specify possible assignments from the left and right image. 
   // in unaryCost{Left_|Right_} there is one more entry than in {left|right}graph_, which is the cost for non-assignment
   std::vector<std::vector<INDEX> > leftGraph_, rightGraph_;

   std::vector<UnaryFactor*> leftUnaryFactor_;
   std::vector<UnaryFactor*> rightUnaryFactor_;

   std::vector<PairwiseFactor*> leftPairwiseFactor_;
   std::vector<PairwiseFactor*> rightPairwiseFactor_;

   std::vector<AssignmentFactor*> leftAssignmentFactor_;
   std::vector<AssignmentFactor*> rightAssignmentFactor_;

   std::vector<std::vector<REAL> > unaryCostLeft_;
   std::vector<std::vector<REAL> > unaryCostRight_;

   std::map<std::pair<INDEX,INDEX>, REAL> leftPottsCost_;
   std::map<std::pair<INDEX,INDEX>, REAL> rightPottsCost_;
};

} // end namespace LP_MP

#endif // LP_MP_COSEGMENTATION_HXX
