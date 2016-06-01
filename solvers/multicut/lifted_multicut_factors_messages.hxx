#ifndef LP_MP_MULTICUT_LIFTED_HXX
#define LP_MP_MULTICUT_LIFTED_HXX

#include "LP_MP.h"

namespace LP_MP {

class LiftedMulticutCutFactor {
public:
   LiftedMulticutCutFactor(const INDEX noCutEdges) 
      : noCutEdges_(noCutEdges),
      noLiftedEdges_(0),
      maxCutEdgeVal_(0.0),
      cutEdgeContrib_(0.0),
      liftedEdgeContrib_(0.0),
      liftedEdgeForcedContrib_(0.0) 
   {}
   //LiftedMulticutCutFactor(const INDEX noCutEdges, const INDEX noLiftedEdges) : noCutEdges_(noCutEdges), noLiftedEdges_(noLiftedEdges) {}

   // the feasible set is: When the ordinary edges are all one (cut), then the lifted edges must be one as well
   // if at least one ordinary edge is zero (not cut), then the lifted edges may be arbitrary
   template<typename REPAM_ARRAY>
   REAL LowerBound(const REPAM_ARRAY& repam) const 
   {
      assert(noCutEdges_ > 0);
      assert(noLiftedEdges_ > 0);
      return cutEdgeContrib_ + liftedEdgeContrib_ + std::max(0.0,std::min(liftedEdgeForcedContrib_,-maxCutEdgeVal_));


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

      assert(std::abs(liftedEdgeContrib - liftedEdgeContrib_) < eps);
      assert(std::abs(liftedEdgeForcedContrib - liftedEdgeForcedContrib_) < eps);
      assert(std::abs(maxCutEdgeVal - maxCutEdgeVal_) < eps);
      assert(std::abs(cutEdgeContrib - cutEdgeContrib_) < eps);

      // all one <=> maxCutEdgeVal <= 0
      return cutEdgeContrib + liftedEdgeContrib + std::max(0.0,std::min(liftedEdgeForcedContrib,-maxCutEdgeVal));

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

      const INDEX edgeIndex = c;
      if(repam[edgeIndex] < maxCutEdgeVal_ || liftedEdgeForcedContrib_ <= -maxCutEdgeVal_) { // do zrobienia: take into account std::max(0.0,...)
         return repam[edgeIndex] + std::max(0.0,std::min(liftedEdgeForcedContrib_,-maxCutEdgeVal_));
      } else { // we must recompute maxCutEdgeVal_ without the active repam value
         REAL maxCutEdgeValExcl = -std::numeric_limits<REAL>::max();
         for(INDEX i=0; i<noCutEdges_; ++i) {
            if(i != c) {
               maxCutEdgeValExcl = std::max(maxCutEdgeValExcl, repam[i]);
            }
         }
         return repam[edgeIndex] + std::max(0.0,std::min(liftedEdgeForcedContrib_,-maxCutEdgeValExcl));
      }
      assert(false);

      //const INDEX edgeIndex = c;
      //const REAL zeroAssignment = cutEdgeContrib - std::min(0.0,repam[edgeIndex]) + liftedEdgeContrib; // no constraints need to be considered, as all are automatically satisfied
      //const REAL oneAssignment = cutEdgeContrib - std::min(0.0,repam[edgeIndex]) + repam[edgeIndex] + liftedEdgeContrib + std::max(0.0,std::min(liftedEdgeForcedContrib,-maxCutEdgeValExcl_[edgeIndex]));
      //return oneAssignment - zeroAssignment;
 

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
      const REAL zeroAssignment = cutEdgeContrib + liftedEdgeContrib; // no constraints need to be considered, as all are automatically satisfied
      const REAL oneAssignment = cutEdgeContrib + repam[edgeIndex] + liftedEdgeContrib + std::max(0.0,std::min(liftedEdgeForcedContrib,-maxCutEdgeVal));
      const REAL diff = repam[edgeIndex] + std::max(0.0,std::min(liftedEdgeForcedContrib,-maxCutEdgeVal));
      assert(std::abs(oneAssignment - zeroAssignment - diff) < eps);
      return diff;

      
      assert(false);
   }

   template<typename REPAM>
   REAL LiftedEdgeBreakpoint(const REPAM& repam, const INDEX c) const
   {
      assert(c < noCutEdges_+noLiftedEdges_);
      assert(c >= noCutEdges_);
      const INDEX edgeIndex = c;

      return repam[edgeIndex] + std::min(0.0,maxCutEdgeVal_) + std::max(0.0,std::min(liftedEdgeForcedContrib_ - std::max(0.0,repam[edgeIndex]),-maxCutEdgeVal_));

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


      const REAL zeroAssignment = cutEdgeContrib + liftedEdgeContrib - std::min(0.0,maxCutEdgeVal); // this assignment allows at most noCutEdges_-1 cut edges to be active
      const REAL oneAssignment = cutEdgeContrib + liftedEdgeContrib + repam[edgeIndex] + std::max(0.0,std::min(liftedEdgeForcedContrib,-maxCutEdgeVal));
      const REAL diff = repam[edgeIndex] + std::min(0.0,maxCutEdgeVal) + std::max(0.0,std::min(liftedEdgeForcedContrib,-maxCutEdgeVal));
      assert(std::abs(oneAssignment - zeroAssignment - diff) < eps);
      return oneAssignment - zeroAssignment;
   }

   // statistics with which evaluation is fast.
   REAL& CutEdgeContrib() { return cutEdgeContrib_; }
   REAL& LiftedEdgeContrib() { return liftedEdgeContrib_; }
   REAL& LiftedEdgeForcedContrib() { return liftedEdgeForcedContrib_; }
   REAL& MaxCutEdgeVal() { return maxCutEdgeVal_; }

   INDEX NoLiftedEdges() const { return noLiftedEdges_; }
   INDEX NoCutEdges() const { return noCutEdges_; }

private:
   // edges are arranged as follows: first come the lifted edges, then the original ones.
   // do zrobienia: make const
   INDEX noLiftedEdges_; // number of lifted edges that have endpoints in different components
   INDEX noCutEdges_; // number of cut edges in the original graph

   // hold these quantities explicitly to avoid having to recompute them after every update -> constant time operations
   REAL maxCutEdgeVal_;
   REAL cutEdgeContrib_;
   REAL liftedEdgeContrib_;
   REAL liftedEdgeForcedContrib_;
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
      //const REAL cutEdgeContribDiff = -std::min(0.0,repamPot[i_]) + std::min(0.0,msg);
      //repamPot.GetFactor()->CutEdgeContrib() += cutEdgeContribDiff;

      repamPot.GetFactor()->CutEdgeContrib() -= std::min(0.0,repamPot[i_]);
      const REAL prevRepam = repamPot[i_];
      repamPot[i_] += msg; 
      repamPot.GetFactor()->CutEdgeContrib() += std::min(0.0,repamPot[i_]);
 
      // update maxCutEdgeVal_, if it needs to be updated.
      if(repamPot[i_] > repamPot.GetFactor()->MaxCutEdgeVal()) {
         repamPot.GetFactor()->MaxCutEdgeVal() = repamPot[i_];
      } else if(prevRepam == repamPot.GetFactor()->MaxCutEdgeVal() && msg < 0.0) {
         // search through all indices to find possibly new maxCutEdgeVal_
         REAL maxCutEdgeVal = -std::numeric_limits<REAL>::max();
         for(INDEX i=0; i<repamPot.GetFactor()->NoCutEdges(); ++i) {
            maxCutEdgeVal = std::max(repamPot[i],maxCutEdgeVal);
         }
         repamPot.GetFactor()->MaxCutEdgeVal() = maxCutEdgeVal;
      }
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
      //const REAL liftedEdgeContribDiff = -std::min(0.0,repamPot[i_]) + std::min(0.0,msg);
      //repamPot.GetFactor()->LiftedEdgeContrib() += liftedEdgeContribDiff;
      //const REAL liftedEdgeForcedContribDiff = -std::max(0.0,repamPot[i_]) + std::max(0.0,msg);
      //repamPot.GetFactor()->LiftedEdgeForcedContrib() += liftedEdgeForcedContribDiff;
      repamPot.GetFactor()->LiftedEdgeContrib() -= std::min(0.0,repamPot[i_]);
      repamPot.GetFactor()->LiftedEdgeForcedContrib() -= std::max(0.0,repamPot[i_]);
      repamPot[i_] += msg; 
      repamPot.GetFactor()->LiftedEdgeContrib() += std::min(0.0,repamPot[i_]);
      repamPot.GetFactor()->LiftedEdgeForcedContrib() += std::max(0.0,repamPot[i_]);
   }

private:
   INDEX i_; // index of lifted edge in the multicut factor. Must be bigger than number of cut edges, smaller than size of cut factor
};

} // end namespace LP_MP

#endif // LP_MP_MULTICUT_LIFTED_HXX

