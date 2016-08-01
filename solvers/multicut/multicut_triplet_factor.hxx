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
   constexpr static INDEX PrimalSize() { return 5; } // primal states: 011 101 110 111 000. Last one receives cost 0 always, however is needed to keep track of primal solution in PropagatePrimal etc.

   template<typename REPAM_ARRAY>
   static REAL LowerBound(const REPAM_ARRAY& repamPot) 
   {
      assert(repamPot.size() == 4);
      return std::min(std::min(std::min( repamPot[0], repamPot[1]), std::min(repamPot[2], repamPot[3])),0.0); // possibly, if we have a SIMD factor, use a simdized minimum
   }

   // if one entry is unknown, set it to true. If one entry is true, set all other to false
   void PropagatePrimal(PrimalSolutionStorage::Element primal) const
   {
      INDEX noTrue = 0;
      INDEX noUnknown = 0;
      INDEX unknownIndex;
      for(INDEX i=0; i<PrimalSize(); ++i) {
         if(primal[i] == true) { 
            ++noTrue;
         }
         if(primal[i] == unknownState) {
            ++noUnknown;
            unknownIndex = i;
         }
      }
      if(noTrue == 0 && noUnknown == 1) {
         primal[unknownIndex] = true;
      } else if(noTrue == 1 && noUnknown > 0) { // should not happen in code currently. Left it for completeness
         for(INDEX i=0; i<PrimalSize(); ++i) {
            if(primal[i] == unknownState) {
               primal[i] = false;
            }
         }
      }
   }

   template<typename REPAM_ARRAY, typename PRIMAL>
   static REAL EvaluatePrimal(const REPAM_ARRAY& repam, const PRIMAL primal)
   {
      assert(repam.size() == 4);
      for(INDEX i=0; i<PrimalSize(); ++i) {
         assert(primal[i] != unknownState);
      }

      REAL cost = 0.0;
      INDEX noActive = 0;
      for(INDEX i=0; i<size(); ++i) {
         cost += primal[i]*repam[i];
         noActive += primal[i];
      }
      noActive += primal[PrimalSize()-1];
      if(noActive != 1) {
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
      //return;
      assert(msg.size() == 1);
      if(primal[(i_+1)%3] == true || primal[(i_+2)%3] == true || primal[3] == true) { // force unary to one
         msg[0] += std::numeric_limits<REAL>::infinity();
         //std::cout << msg.GetLeftFactor()->operator[](0) << "\n";
         //assert(msg.GetLeftFactor()->operator[](0) == -std::numeric_limits<REAL>::infinity() );
      } else if(primal[i_] == true || primal[4] == true) { // force unary to zero 
         //assert(msg.GetLeftFactor()->operator[](0) != std::numeric_limits<REAL>::infinity() );
         msg[0] -= std::numeric_limits<REAL>::infinity(); 
         //std::cout << msg.GetLeftFactor()->operator[](0) << "\n";
         //assert(msg.GetLeftFactor()->operator[](0) == std::numeric_limits<REAL>::infinity() );
      } else { // compute message on unknown values only. No entry of primal is true
         assert(4 == r->size());
         assert(RIGHT_FACTOR::PrimalSize() == 5);
         std::array<REAL,RIGHT_FACTOR::PrimalSize()> restrictedPot;
         restrictedPot.fill(std::numeric_limits<REAL>::infinity());
         for(INDEX i=0; i<size(); ++i) {
            assert(primal[i] != true);
            if(primal[i] == unknownState) {
               restrictedPot[i] = rightPot[i];
            }
         }
         if(primal[4] == unknownState) {
            restrictedPot[4] = 0.0;
         }
         msg[0] -= std::min(std::min(restrictedPot[(i_+1)%3], restrictedPot[(i_+2)%3]), restrictedPot[3]) - std::min(restrictedPot[i_],restrictedPot[4]);
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

   template<typename LEFT_FACTOR, typename RIGHT_FACTOR>
   void ComputeRightFromLeftPrimal(PrimalSolutionStorage::Element leftPrimal, LEFT_FACTOR* l, typename PrimalSolutionStorage::Element rightPrimal, RIGHT_FACTOR* r) const
   {
      if(leftPrimal[0] == true) {
         rightPrimal[i_] = false;
         rightPrimal[4] = false;
      } else {
         assert(leftPrimal[0] == false);
         rightPrimal[(i_+1)%3] = false;
         rightPrimal[(i_+2)%3] = false;
         rightPrimal[3] = false;
      }
      return;
   }

   template<typename LEFT_FACTOR, typename RIGHT_FACTOR>
   bool CheckPrimalConsistency(PrimalSolutionStorage::Element leftPrimal, LEFT_FACTOR* l, typename PrimalSolutionStorage::Element rightPrimal, RIGHT_FACTOR* r) const
   {
      return true;
      assert(leftPrimal[0] != unknownState);
      if(leftPrimal[0] == false) {
         if(!(rightPrimal[i_] == true || rightPrimal[4] == true)) {
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

