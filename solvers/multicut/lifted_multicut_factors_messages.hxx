#ifndef LP_MP_MULTICUT_LIFTED_HXX
#define LP_MP_MULTICUT_LIFTED_HXX

#include "LP_MP.h"

namespace LP_MP {

class LiftedMulticutCutFactor {
public:
   LiftedMulticutCutFactor(const INDEX noCutEdges) : noCutEdges_(noCutEdges), noLiftedEdges_(0) {}
   LiftedMulticutCutFactor(const INDEX noCutEdges, const INDEX noLiftedEdges) : noCutEdges_(noCutEdges), noLiftedEdges_(noLiftedEdges) {}

   // the feasible set is: When the ordinary edges are all one (cut), then the lifted edges must be one as well
   // if at least one ordinary edge is zero (not cut), then the lifted edges may be arbitrary
   template<typename REPAM_ARRAY>
   REAL LowerBound(const REPAM_ARRAY& repam) const 
   {
      assert(noCutEdges_ > 0);
      assert(noLiftedEdges_ > 0);

      std::cout << "lifted factor repam: ";
      for(INDEX i=0; i<noCutEdges_; ++i) {
         std::cout << repam[i] << ", ";
      }
      std::cout << ";;; ";
      for(INDEX i=0; i<noLiftedEdges_; ++i) {
         std::cout << repam[i+noCutEdges_] << ", ";
      }
      std::cout << "\n";
      // do zrobienia: recompute all statistics here from scratch, once they are introduced
      assert(repam.size() == size());

      // do zrobienia: start at index 1, initialize with zero value
      REAL maxCutEdgeVal = -std::numeric_limits<REAL>::max();
      REAL cutEdgeContrib = 0.0;
      for(INDEX i=0; i<noCutEdges_; ++i) {
         maxCutEdgeVal = std::max(maxCutEdgeVal, repam[i]);
         cutEdgeContrib += std::min(0.0,repam[i]);
      }

      // do zrobienia: start at index 1, initialize with zero value
      REAL liftedEdgeContrib = 0.0; // the value which the lifted edge can contribute to the lower bound
      REAL liftedEdgeForcedContrib = 0.0;
      for(INDEX i=0; i<noLiftedEdges_; ++i) {
         liftedEdgeContrib += std::min(0.0,repam[i + noCutEdges_]);
         liftedEdgeForcedContrib += std::max(0.0,repam[i + noCutEdges_]);
      }

      // all one <=> maxCutEdgeVal <= 0
      // do zrobienia: devise single formula without case study

      return cutEdgeContrib + liftedEdgeContrib + std::max(0.0,std::min(liftedEdgeForcedContrib,-maxCutEdgeVal));

      /*
      if(maxCutEdgeVal >= 0.0) { // there is a cut edge with value 0
         return cutEdgeContrib + liftedEdgeContrib;
      } else { // all cut edges want to be one
         // do zrobienia: correct?
         if(-maxCutEdgeVal <= liftedEdgeForcedContrib) { // we force a cut edge to be zero, lifted edges are free
            return cutEdgeContrib - maxCutEdgeVal + liftedEdgeContrib;
         } else {
            return cutEdgeContrib + liftedEdgeContrib + liftedEdgeForcedContrib; // all lifted edges are one
         }
      }
      */
   }

   INDEX size() const { return noCutEdges_ + noLiftedEdges_; }

   template<typename REPAM_ARRAY>
   REAL EvaluatePrimal(const REPAM_ARRAY& repam, const PrimalSolutionStorage::Element primal) const
   {
      assert(repam.size() == size());
      return 1000000.0;
   }

   // do zrobienia: use own reparametrization storage and make it bigger here as well.
   void IncreaseLifted()
   {
      ++noLiftedEdges_;
   }

   // do zrobienia: make this more efficient by keeping track of needed informations
   // compute by how much we can change the cut edge's cost
   template<typename REPAM>
   REAL CutEdgeBreakpoint(const REPAM& repam, const INDEX c) const 
   { 
      assert(c < noCutEdges_);

      // do zrobienia: start at index 1, initialize with zero value
      REAL maxCutEdgeVal = -std::numeric_limits<REAL>::max();
      REAL cutEdgeContrib = 0.0;
      for(INDEX i=0; i<noCutEdges_; ++i) {
         if(i != c) {
            maxCutEdgeVal = std::max(maxCutEdgeVal, repam[i]);
            cutEdgeContrib += std::min(0.0,repam[i]);
         }
      }

      // do zrobienia: start at index 1, initialize with 0 value
      REAL liftedEdgeContrib = 0.0; // the value which the lifted edge can contribute to the lower bound
      REAL liftedEdgeForcedContrib = 0.0;
      for(INDEX i=0; i<noLiftedEdges_; ++i) {
         liftedEdgeContrib += std::min(0.0,repam[i + noCutEdges_]);
         liftedEdgeForcedContrib += std::max(0.0,repam[i + noCutEdges_]);
      }


      // approach: compute the cost for assignment =1 and assignment -0 of current variable c and let change be difference of those two values
      const INDEX edgeIndex = c;
      const REAL zeroAssignment = cutEdgeContrib + liftedEdgeContrib; // no constraints need to be considered, as all are automatically satisfied
      const REAL oneAssignment = cutEdgeContrib + repam[edgeIndex] + liftedEdgeContrib + std::max(0.0,std::min(liftedEdgeForcedContrib,-maxCutEdgeVal));
      return oneAssignment - zeroAssignment;

      
      assert(false);
   }

   template<typename REPAM>
   REAL LiftedEdgeBreakpoint(const REPAM& repam, const INDEX c) const
   {
      assert(c < noCutEdges_+noLiftedEdges_);
      assert(c >= noCutEdges_);

      REAL maxCutEdgeVal = -std::numeric_limits<REAL>::max();
      REAL cutEdgeContrib = 0.0;
      for(INDEX i=0; i<noCutEdges_; ++i) {
         maxCutEdgeVal = std::max(maxCutEdgeVal, repam[i]);
         cutEdgeContrib += std::min(0.0,repam[i]);
      }

      REAL liftedEdgeContrib = 0.0; // the value which the lifted edge can contribute to the lower bound
      REAL liftedEdgeForcedContrib = 0.0;
      for(INDEX i=0; i<noLiftedEdges_; ++i) {
         if(i + noCutEdges_ != c) {
            liftedEdgeContrib += std::min(0.0,repam[i + noCutEdges_]);
            liftedEdgeForcedContrib += std::max(0.0,repam[i + noCutEdges_]);
         }
      }

      const INDEX edgeIndex = c;

      const REAL zeroAssignment = cutEdgeContrib + liftedEdgeContrib - std::min(0.0,maxCutEdgeVal); // this assignment allows at most noCutEdges_-1 cut edges to be active
      //const REAL liftedEdgeForcedContribOne = liftedEdgeForcedContrib + std::max(0.0,repam[edgeIndex]);
      const REAL oneAssignment = cutEdgeContrib + liftedEdgeContrib + repam[edgeIndex] + std::max(0.0,std::min(liftedEdgeForcedContrib,-maxCutEdgeVal));
      return oneAssignment - zeroAssignment;


      if(maxCutEdgeVal >= 0.0) { // there is a cut edge with value 0, lifted edges not constrained
         return repam[edgeIndex];
      } else { // all cut edges want to be one
         // how much does it take to make cut edge with largest (negative) cost to be zero?
         // make liftedEdgeForcedContrib = -maxCutEdgeVal
         return std::min(repam[edgeIndex],0.0) + liftedEdgeForcedContrib + maxCutEdgeVal;
         //return repam[edgeIndex] + std::min(0.0,maxCutEdgeVal + liftedEdgeForcedContrib);
      }

      assert(false);
   }


private:
   // edges are arranged as follows: first come the lifted edges, then the original ones.
   // do zrobienia: make const
   INDEX noLiftedEdges_; // number of lifted edges that have endpoints in different components
   INDEX noCutEdges_; // number of cut edges in the original graph
};

// message between edge of the base graph and lifted multicut factor
class CutEdgeLiftedMulticutFactorMessage {
public:
   CutEdgeLiftedMulticutFactorMessage(const INDEX i) : i_(i) {}

   constexpr static INDEX size() { return 1; }

   template<typename RIGHT_FACTOR, typename G1, typename G2>
   void
   ReceiveMessageFromRight(RIGHT_FACTOR* const r, const G1& rightPot, G2& msg) 
   {
      assert(msg.size() == 1);
      msg[0] -= r->CutEdgeBreakpoint(rightPot,i_);
   }

   template<typename LEFT_FACTOR, typename G1, typename G3>
   void
   SendMessageToRight(LEFT_FACTOR* const l, const G1& leftPot, G3& msg, const REAL omega)
   {
      assert(msg.size() == 1);
      assert(leftPot.size() == 1);
      msg[0] -= omega*leftPot[0];
   }

   template<typename G>
   void RepamLeft(G& repamPot, const REAL msg, const INDEX msg_dim)
   {
      assert(msg_dim == 0);
      assert(repamPot.size() == 1);
      repamPot[0] += msg;
   }
   template<typename G>
   void RepamRight(G& repamPot, const REAL msg, const INDEX msg_dim)
   {
      assert(msg_dim == 0);
      repamPot[i_] += msg; 
   }

private:
   INDEX i_; // index of cut edge in the multicut factor
};

// message between lifted edge and lifted multicut factor
class LiftedEdgeLiftedMulticutFactorMessage {
public:
   LiftedEdgeLiftedMulticutFactorMessage(const INDEX i) : i_(i) {}

   constexpr static INDEX size() { return 1; }

   template<typename RIGHT_FACTOR, typename G1, typename G2>
   void
   ReceiveMessageFromRight(RIGHT_FACTOR* const r, const G1& rightPot, G2& msg) 
   {
      assert(msg.size() == 1);
      msg[0] -= r->LiftedEdgeBreakpoint(rightPot,i_);
   }

   template<typename LEFT_FACTOR, typename G1, typename G3>
   void
   SendMessageToRight(LEFT_FACTOR* const l, const G1& leftPot, G3& msg, const REAL omega)
   {
      assert(msg.size() == 1);
      assert(leftPot.size() == 1);
      msg[0] -= omega*leftPot[0];
   }

   template<typename G>
   void RepamLeft(G& repamPot, const REAL msg, const INDEX msg_dim)
   {
      assert(msg_dim == 0);
      assert(repamPot.size() == 1);
      repamPot[0] += msg;
   }
   template<typename G>
   void RepamRight(G& repamPot, const REAL msg, const INDEX msg_dim)
   {
      assert(msg_dim == 0);
      repamPot[i_] += msg; 
   }


private:
   INDEX i_; // index of lifted edge in the multicut factor. Must be bigger than number of cut edges, smaller than size of cut factor
};

} // end namespace LP_MP

#endif // LP_MP_MULTICUT_LIFTED_HXX

