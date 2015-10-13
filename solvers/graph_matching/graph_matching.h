#ifndef LP_MP_GRAPH_MATCHING_SOLVER
#define LP_MP_GRAPH_MATCHING_SOLVER

#include "factors_messages.hxx"
#include "LP_MP.h"
#include "factors/multiplex_factor.hxx"
#include "factors/simplex_factor_simd.hxx"
#include "factors/pairwise_simplex_factor_simd.hxx"
#include "const_array_types.h"
#include "messages/multiplex_marg_message.hxx"
#include "messages/multiplex_message_simd.hxx"
#include "messages/equality_message.hxx"
#include "problem_construction_helper.hxx"

#include "MinCost/MinCost.h"

#include <iostream>
#include <fstream>
#include <string>
#include <tuple>
#include <exception>
#include <algorithm>

namespace LP_MP {

class GraphMatching
{
public:
   GraphMatching();
   ~GraphMatching();
   void AddLeftNode(const INDEX nodeNumber, const INDEX min_capacity, const INDEX max_capacity);
   void AddRightNode(const INDEX nodeNumber, const INDEX min_capacity, const INDEX max_capacity);
   void AddAssignmentCost(const INDEX leftNode, const INDEX rightNode, const REAL cost);
   void AddPairwiseCost(INDEX leftNode1, INDEX rightNode1, INDEX leftNode2, INDEX rightNode2, const REAL cost);
   void AddLeftPairwiseCost(INDEX leftNode1, INDEX rightNode1, INDEX leftNode2, INDEX rightNode2, const REAL cost);
   void AddRightPairwiseCost(INDEX rightNode1, INDEX leftNode1, INDEX rightNode2, INDEX leftNode2, const REAL cost);
   void AddTernaryCost(const INDEX leftNode1, const INDEX rightNode1, const INDEX leftNode2, const INDEX rightNode2, const INDEX leftNode3, const INDEX rightNode3, const REAL cost);
   void AddLeftTernaryCost(const INDEX leftNode1, const INDEX rightNode1, const INDEX leftNode2, const INDEX rightNode2, const INDEX leftNode3, const INDEX rightNode3, const REAL cost);
   void AddRightTernaryCost(const INDEX rightNode1, const INDEX leftNode1, const INDEX rightNode2, const INDEX leftNode2, const INDEX rightNode3, const INDEX leftNode3, const REAL cost);
   // do zrobienia: think about potts terms for multiplex unaries
   void AddLeftPottsTerm(const INDEX leftNode1, const INDEX leftNode2, const REAL cost);
   void AddRightPottsTerm(const INDEX rightNode1, const INDEX rightNode2, const REAL cost);
   void Solve(const INDEX nIter);

   std::vector<int> GetLeftAssignment();
   std::vector<int> GetRightAssignment();
   REAL dualBound() {return lp_->LowerBound();};
   REAL evalLeft(const std::vector<int>& assignment);
   REAL evalRight(const std::vector<int>& assignment);
   REAL primalBound() { return std::min(evalLeft(GetLeftAssignment()), evalRight(GetRightAssignment())); };
   std::vector<INDEX> arg() {return std::vector<INDEX>(0);}

private:
   //specification for underlying factor/message-network
   typedef UnaryLoop UnaryLoopType;
   typedef PairwiseLoop<0> LeftLoopType;  // 0 is right from test with non-homogenous data
   typedef PairwiseLoop<1> RightLoopType; // 1 is right from test with non-homogenous data

   // normal
   typedef MultiplexMargMessage<UnaryLoopType,LeftLoopType,true,false> LeftMargMessage;
   typedef MultiplexMargMessage<UnaryLoopType,RightLoopType,true,false> RightMargMessage;
   // SIMD
   //typedef MultiplexMessageSIMD<UnaryLoopType,LeftLoopType,true,false,true,false> LeftMargMessage; // for SIMD
   //typedef MultiplexMessageSIMD<UnaryLoopType,RightLoopType,true,false,false,true> RightMargMessage; // for SIMD


   struct FMC; // forward declaration
   
   typedef MultiplexFactor<std::vector<REAL>, const_ones_array, const_one> Simplex;
   typedef FactorContainer<Simplex, ExplicitRepamStorage, FMC, 0 > UnaryFactor;
   //typedef FactorContainer<Simplex, ImplicitRepamStorage, FMC, 1 > PairwiseFactor; // implicit repam storage
   typedef FactorContainer<Simplex, ExplicitRepamStorage, FMC, 1 > PairwiseFactor; // normal
   //typedef FactorContainer<SimplexFactorSIMD<false>, ExplicitRepamStorageSIMD, FMC, 1 > PairwiseFactor; // normal
   //typedef FactorContainer<PairwiseSimplexFactorSIMD, PairwiseStorageSIMD, FMC, 1 > PairwiseFactor; // for SIMD
   typedef MessageContainer<EqualityMessage, FixedMessageStorage<1>, FMC, 0 > AssignmentConstraintMessage;
   typedef MessageContainer<LeftMargMessage, MessageStorageSIMD, FMC, 1 > UnaryPairwiseMessageLeft; // unary has left_dim entries
   typedef MessageContainer<RightMargMessage, MessageStorageSIMD, FMC, 2 > UnaryPairwiseMessageRight; // unary has right_dim entries

   struct FMC {
      using factor_list = meta::list< UnaryFactor, PairwiseFactor >;
      // fourth and fifth entry are message holding containers for the left and right factor type: ugly: the types the container are holding could be deduced automatically, but they are not. meta::list cannot easily hold class templates
      using msg_list = meta::list< 
                                   meta::list< AssignmentConstraintMessage, meta::size_t<0>, meta::size_t<0>, std::vector<AssignmentConstraintMessage*>, std::vector<AssignmentConstraintMessage*> >,
                                   meta::list< UnaryPairwiseMessageLeft,  meta::size_t<0>, meta::size_t<1>, std::vector<UnaryPairwiseMessageLeft*>, FixedSizeContainer<UnaryPairwiseMessageLeft*,1> >,
                                   meta::list< UnaryPairwiseMessageRight, meta::size_t<0>, meta::size_t<1>, std::vector<UnaryPairwiseMessageRight*>, FixedSizeContainer<UnaryPairwiseMessageRight*,1> >
                                      >;
   };



   static std::vector<INDEX> ConstructMultiplexCapacityVector(const std::vector<INDEX>& graph, const std::pair<INDEX,INDEX> nodeCapacity);
   void AddNode(const INDEX nodeNumber, const INDEX min_capacity, const INDEX max_capacity, std::vector<std::vector<INDEX> >& graph, std::vector<std::pair<INDEX,INDEX> >& nodeCapacity);
   void AddToPairwisePotential(
         INDEX leftNode1, INDEX leftNode2, 
         INDEX rightNode1, INDEX rightNode2, 
         const REAL cost, 
         std::map<std::pair<INDEX,INDEX>, std::map<std::pair<INDEX,INDEX>, REAL> >& pairwisePotentials);
   void AddToTernaryPotential(
         INDEX leftNode1, INDEX leftNode2, INDEX leftNode3,
         INDEX rightNode1, INDEX rightNode2, INDEX rightNode3, 
         const REAL cost, 
         std::map<std::tuple<INDEX,INDEX,INDEX>, std::map<std::tuple<INDEX,INDEX,INDEX>, REAL> >& ternaryPotentials);
   void ConstructUnaryCost(
         std::vector<UnaryFactor*>& unaryMultiplex,
         std::vector<std::vector<REAL> >& unaryCost,
         const std::vector<std::pair<INDEX,INDEX> >& nodeCapacity,
         const std::vector<std::vector<INDEX> >& graph);
   void ConstructPairwiseCost(
         std::map<std::pair<INDEX,INDEX>, std::map<std::pair<INDEX,INDEX>, REAL> >& pairwisePotentials,
         std::map<std::pair<INDEX,INDEX>, REAL>& PottsCost,
         const std::vector<std::vector<INDEX> >& graph,
         const std::vector<std::pair<INDEX,INDEX> >& nodeCapacity,
         const std::vector<std::vector<INDEX> >& inverseGraph,
         const std::vector<UnaryFactor*>& unaryMultiplex,
         std::map<std::pair<INDEX,INDEX>, PairwiseFactor*>& pairwiseMultiplex);
   void ConstructTernaryCost(
         std::map<std::tuple<INDEX,INDEX,INDEX>, std::map<std::tuple<INDEX,INDEX,INDEX>, REAL> > ternaryPotentials,
         const std::vector<std::vector<INDEX> >& graph,
         const std::vector<std::pair<INDEX,INDEX> >& nodeCapacity,
         const std::vector<std::vector<INDEX> >& inverseGraph,
         const std::map<std::pair<INDEX,INDEX>, PairwiseFactor*>& pairwiseMultiplexFactor);
   void ConstructAssignmentConstraints();

   std::vector<int> GetAssignment(
         const std::vector<std::vector<REAL> >& cost, 
         const std::vector<std::vector<INDEX> >& graph,
         const std::vector<std::pair<INDEX,INDEX> >& nodeCapacity,
         const std::vector<std::pair<INDEX,INDEX> >& inverseNodeCapacity);
   REAL eval(
      const std::vector<int>& assignment,
      const std::vector<std::vector<INDEX> >& graph, 
      const std::vector<std::vector<REAL> >& unaryCost,
      const std::map<std::pair<INDEX,INDEX>, std::map<std::pair<INDEX,INDEX>, REAL> >& pairwisePotentials);


   LP* lp_; // do zrobienia: consider const
   std::vector<std::vector<INDEX> >  leftGraph_, rightGraph_;
   std::vector<std::pair<INDEX,INDEX> > leftNodeCapacity_, rightNodeCapacity_;

   // convert to second order MRF with assignment constraints
   std::vector<std::vector<REAL> > unaryCostLeft_;
   std::vector<std::vector<REAL> > unaryCostRight_;
   std::vector<UnaryFactor*> leftGraphUnaryMultiplex_;
   std::vector<UnaryFactor*> rightGraphUnaryMultiplex_;
   std::map<std::pair<INDEX,INDEX>, std::map<std::pair<INDEX,INDEX>, REAL> > leftPairwisePotentials_;
   std::map<std::pair<INDEX,INDEX>, std::map<std::pair<INDEX,INDEX>, REAL> > rightPairwisePotentials_;
   std::map<std::pair<INDEX,INDEX>, REAL> leftPottsCost_;
   std::map<std::pair<INDEX,INDEX>, REAL> rightPottsCost_;
   std::map<std::pair<INDEX,INDEX>, PairwiseFactor*> leftGraphPairwiseMultiplex_;
   std::map<std::pair<INDEX,INDEX>, PairwiseFactor*> rightGraphPairwiseMultiplex_;
   std::map<std::tuple<INDEX,INDEX,INDEX>, std::map<std::tuple<INDEX,INDEX,INDEX>, REAL> > leftTernaryPotentials_;
   std::map<std::tuple<INDEX,INDEX,INDEX>, std::map<std::tuple<INDEX,INDEX,INDEX>, REAL> > rightTernaryPotentials_;
};

} // end namespace LP_MP

#endif // LP_MP_GRAPH_MATCHING_SOLVER
