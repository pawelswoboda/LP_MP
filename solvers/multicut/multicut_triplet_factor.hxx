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

   template<typename REPAM_ARRAY>
   void MaximizePotential(const REPAM_ARRAY& repam) {};

   constexpr static INDEX size() { return 4; }

   template<typename REPAM_ARRAY>
   REAL LowerBound(const REPAM_ARRAY& repamPot) const 
   {
      assert(repamPot.size() == 4);
      return std::min(std::min(std::min( repamPot[0], repamPot[1]), std::min(repamPot[2], repamPot[3])),0.0); // possibly, if we have a SIMD factor, use a simdized minimum
   }

   template<typename REPAM_ARRAY, typename PRIMAL>
   REAL EvaluatePrimal(const REPAM_ARRAY& repam, const PRIMAL primal) const
   {
      assert(repam.size() == 4);
      return 3e111;
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
   ReceiveMessageFromRight(RIGHT_FACTOR* const r, const G1& rightPot, G2& msg) 
   {
      assert(msg.size() == 1);
      msg[0] -= std::min(std::min(rightPot[(i_+1)%3], rightPot[(i_+2)%3]), rightPot[3]) - std::min(rightPot[i_],0.0);
   }

   template<typename LEFT_FACTOR, typename G1, typename G3, MessageSendingType MST_TMP = MST>
   //typename std::enable_if<MST_TMP == MessageSending::SRMP,void>::type
   void
   SendMessageToRight(LEFT_FACTOR* const l, const G1& leftPot, G3& msg, const REAL omega)
   {
      assert(msg.size() == 1);
      static_assert(MST_TMP == MST,"");
      msg[0] -= omega*leftPot[0];
   }

   template<typename G>
   void RepamLeft(G& repamPot, const REAL msg, const INDEX msg_dim)
   {
      assert(msg_dim == 0);
      repamPot[0] += msg;
   }
   template<typename G>
   void RepamRight(G& repamPot, const REAL msg, const INDEX msg_dim)
   {
      assert(msg_dim == 0);
      // the labeling with two 1s in them
      repamPot[(i_+1)%3] += msg;
      repamPot[(i_+2)%3] += msg;
      // for the 111 labeling which is always affected
      repamPot[3] += msg; 
   }

   // compute primal functions, how to do it? look into ultra-metric rounding
   void ComputeRightFromLeftPrimal(const bool leftPrimal, typename PrimalSolutionStorage::Element rightPrimal)
   {
      // do zrobienia
   }

private:
   const INDEX i_;
};

} // end namespace LP_MP

#endif // LP_MP_MULTICUT_TRIPLET_FACTOR

