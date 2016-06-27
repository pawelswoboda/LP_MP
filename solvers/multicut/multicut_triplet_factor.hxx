#ifndef LP_MP_MULTICUT_TRIPLET_FACTOR_HXX
#define LP_MP_MULTICUT_TRIPLET_FACTOR_HXX

#include "LP_MP.h"

namespace LP_MP {

// encoding with 4 entries corresponding to the four possible labelings in this order:
// 011 101 110 111. Labeling 000 is always assigned cost 0. Labelings 100 010 001 are forbidden.
class MulticutTripletFactor
{
public:
   MulticutTripletFactor() {};

   using TripletEdges = std::array<std::array<INDEX,2>,3>;
   static TripletEdges SortEdges(const INDEX i1, const INDEX i2, const INDEX i3)
   {
      std::array<INDEX,3> ti{i1,i2,i3}; // the node indices of the triplet
      std::sort(ti.begin(), ti.end()); // do zrobienia: use faster sorting
      assert(ti[0] < ti[1] < ti[2]);
      TripletEdges te{{{ti[0],ti[1]},{ti[0],ti[2]},{ti[1],ti[2]}}}; 
      return te;
   }

   constexpr static INDEX size() { return 4; }

   template<typename REPAM_ARRAY>
   static REAL LowerBound(const REPAM_ARRAY& repamPot) 
   {
      assert(repamPot.size() == 4);
      return std::min(std::min(std::min( repamPot[0], repamPot[1]), std::min(repamPot[2], repamPot[3])),0.0); // possibly, if we have a SIMD factor, use a simdized minimum
   }

   template<typename REPAM_ARRAY, typename PRIMAL>
   static REAL EvaluatePrimal(const REPAM_ARRAY& repam, const PRIMAL primal)
   {
      assert(repam.size() == 4);
      for(INDEX i=0; i<4; ++i) {
         assert(primal[i] != unknownState);
      }

      REAL cost = 0.0;
      INDEX noActive = 0;
      for(INDEX i=0; i<4; ++i) {
         cost += primal[i]*repam[i];
         noActive += primal[i];
      }
      if(noActive > 1) {
         assert(false);
         return std::numeric_limits<REAL>::infinity();
      }
      return cost;
   }

};


// message between unary factor and triplet factor
//enum class MessageSending { SRMP, MPLP }; // do zrobienia: place this possibly more global, also applies to pairwise factors in MRFs
template<MessageSendingType MST = MessageSendingType::SRMP>
class MulticutUnaryTripletMessage
{
public:
   // i is the index in the triplet factor
   MulticutUnaryTripletMessage(const INDEX i) : i_(i) 
   { 
      assert(i < 3); 
   }; 
   ~MulticutUnaryTripletMessage()
   {
      static_assert(MST == MessageSendingType::SRMP,"");
   }

   constexpr static INDEX size() { return 1; }

   template<typename RIGHT_FACTOR, typename G1, typename G2, MessageSendingType MST_TMP = MST>
   //typename std::enable_if<MST_TMP == MessageSending::SRMP,void>::type
   void
   ReceiveMessageFromRight(RIGHT_FACTOR* const r, const G1& rightPot, G2& msg) const
   {
      assert(msg.size() == 1);
      msg[0] -= std::min(std::min(rightPot[(i_+1)%3], rightPot[(i_+2)%3]), rightPot[3]) - std::min(rightPot[i_],0.0);
   }

   template<typename RIGHT_FACTOR, typename G1, typename G2, MessageSendingType MST_TMP = MST>
   //typename std::enable_if<MST_TMP == MessageSending::SRMP,void>::type
   void
   ReceiveRestrictedMessageFromRight(RIGHT_FACTOR* const r, const G1& rightPot, G2& msg, PrimalSolutionStorage::Element primal) const
   {
      assert(msg.size() == 1);
      if(primal[(i_+1)%3] == true || primal[(i_+2)%3] == true || primal[3] == true) { // force unary to one
         msg[0] += std::numeric_limits<REAL>::infinity();
         assert(msg.GetLeftFactor()->operator[](0) == -std::numeric_limits<REAL>::infinity() );
      } else if(primal[i_] == true || (primal[0] == false && primal[1] == false && primal[2] == false && primal[3] == false)) { // force unary to zero 
         msg[0] -= std::numeric_limits<REAL>::infinity(); 
         assert(msg.GetLeftFactor()->operator[](0) == std::numeric_limits<REAL>::infinity() );
      } else { // compute message on unknown values only. No entry of primal is true
         assert(4 == r->size());
         std::array<REAL,4> restrictedPot;
         restrictedPot.fill(std::numeric_limits<REAL>::infinity());
         for(INDEX i=0; i<4; ++i) {
            assert(primal[i] != true);
            if(primal[i] == unknownState) {
               restrictedPot[i] = rightPot[i];
            }
         }
         msg[0] -= std::min(std::min(restrictedPot[(i_+1)%3], restrictedPot[(i_+2)%3]), restrictedPot[3]) - std::min(restrictedPot[i_],0.0);
      }
   }

   template<typename LEFT_FACTOR, typename G1, typename G3, MessageSendingType MST_TMP = MST>
   //typename std::enable_if<MST_TMP == MessageSending::SRMP,void>::type
   void
   SendMessageToRight(LEFT_FACTOR* const l, const G1& leftPot, G3& msg, const REAL omega) const
   {
      assert(msg.size() == 1);
      static_assert(MST_TMP == MST,"");
      msg[0] -= omega*leftPot[0];
   }

   template<typename G>
   void RepamLeft(G& repamPot, const REAL msg, const INDEX msg_dim) const
   {
      assert(msg_dim == 0);
      repamPot[0] += msg;
   }
   template<typename G>
   void RepamRight(G& repamPot, const REAL msg, const INDEX msg_dim) const
   {
      assert(msg_dim == 0);
      // the labeling with two 1s in them
      repamPot[(i_+1)%3] += msg;
      repamPot[(i_+2)%3] += msg;
      // for the 111 labeling which is always affected
      repamPot[3] += msg; 
   }

   void ComputeRightFromLeftPrimal(PrimalSolutionStorage::Element leftPrimal, typename PrimalSolutionStorage::Element rightPrimal) const
   {
      if(leftPrimal[0] == true) { 
         // check whether one of 011,101,110 can be inferred
         if(rightPrimal[3] != false && rightPrimal[i_] == false) { // this means that some other message has set current variable to zero, but then we can infer that one of 110,101,011 must be true
            if(rightPrimal[(i_+1)%3] == false && rightPrimal[(i_+2)%3] != false) {
               rightPrimal[(i_+2)%3] = true;
            } else if(rightPrimal[(i_+2)%3] == false && rightPrimal[(i_+1)%3] != false) { 
               rightPrimal[(i_+1)%3] = true;
            }
         }
         rightPrimal[i_] = false;
      } else if( leftPrimal[0] == false) {
         // check, whether some other variable is one. This is true, whenever exactly one other 0-configuration is forbidden. If so, we can infer one of 011,101,110
         if(rightPrimal[3] != false && ((rightPrimal[(i_+1)%3] == false || rightPrimal[(i_+2)%3] != false) || (rightPrimal[(i_+2)%3] == false || rightPrimal[(i_+1)%3] != false)) ) {
            rightPrimal[i_] = true;
         }
         rightPrimal[(i_+1)%3] = false;
         rightPrimal[(i_+2)%3] = false;
         rightPrimal[3] = false;
      } else { // we assume that unary factor has been set in rounding previously.
         assert(false);
      }
      // now check if 111 can be inferred:
      // if 011,101,110 are false, but 111 is not, then make it true. We know in this case that all the unaries are true
      if(rightPrimal[0] == false && rightPrimal[1] == false && rightPrimal[2] == false && rightPrimal[3] != false) {
         rightPrimal[3] = true;
      } 
   }

   bool CheckPrimalConsistency(PrimalSolutionStorage::Element leftPrimal, PrimalSolutionStorage::Element rightPrimal) const
   {
      assert(leftPrimal[0] != unknownState);
      if(leftPrimal[0] == false) {
         if(!(rightPrimal[i_] == true || (rightPrimal[0] == false && rightPrimal[1] == false && rightPrimal[2] == false && rightPrimal[3] == false))) {
            return false;
         }
      } else {
         if(!(rightPrimal[(i_+1)%3] == true || rightPrimal[(i_+2)%3] == true || rightPrimal[3] == true)) {
            return false;
         }
      }
      return true;
   }

private:
   const SHORT_INDEX i_;
   // do zrobienia: measure whether precomputing (i_+1)%3 and (i_+2)%3 would make algorithm faster.
   // all values can be stored in one byte
};

} // end namespace LP_MP

#endif // LP_MP_MULTICUT_TRIPLET_FACTOR

