#ifndef LP_MP_MULTICUT_MESSAGE_HXX
#define LP_MP_MULTICUT_MESSAGE_HXX

#include "LP_MP.h"
#include "multicut_unary_factor.hxx"
#include "multicut_triplet_factor.hxx"

namespace LP_MP {
 
// left factor must be MulticutUnaryFactor and right factor must be MulticutTripletFactor
// possibly templatize for index i_
// templatize for either SRMP- or MPLP-type message passing, i.e. uanry factors are active or triplet factors are active
enum class MessageSending { SRMP, MPLP }; // do zrobienia: place this possibly more global, also applies to pairwise factors in MRFs
template<MessageSending MST = MessageSending::SRMP>
class MulticutUnaryTripletMessage
{
public:
   MulticutUnaryTripletMessage(const INDEX i) : i_(i) { assert(i < 3); }; // i is the index in the triplet factor

   template<typename RIGHT_FACTOR, typename G1, typename G2, MessageSending MST_TMP = MST>
   typename std::enable_if<MST_TMP == MessageSending::SRMP,void>::type
   ReceiveMessageFromRight(RIGHT_FACTOR* const r, const G1& rightPot, G2& msg) 
   {
      static_assert(MST_TMP == MST,"");
      const auto sortIndices = r->SortIndices(rightPot);
      // now adjust costs such that labeling just stays optimal for modified rightPot
      const REAL x = rightPot[sortIndices[0]] + rightPot[sortIndices[1]];

      // investigate all possiblities: labelings + position in labeling -> 3x3 decisions
      // all three possible states (up to permutation)
      // do zrobienia: devise one formula which takes into account all those decisions. faster?
      if(x > 0.0) { // labeling 000
         if(sortIndices[0] == i_) {
            msg[0] += -x;
         } else if(sortIndices[1] == i_) {
            msg[0] += -x;
         } else {
            msg[0] += -rightPot[i_] - rightPot[sortIndices[0]];
         }
      } else if(rightPot[sortIndices[2]] > 0) { // labeling 110
         if(sortIndices[0] == i_) {
            msg[0] += -rightPot[i_] + std::min(-rightPot[sortIndices[1]], rightPot[sortIndices[2]]);
         } else if(sortIndices[1] == i_) {
            msg[0] += -rightPot[i_] + std::min(-rightPot[sortIndices[0]], rightPot[sortIndices[2]]);
         } else {
            msg[0] += -rightPot[i_] + std::max(rightPot[sortIndices[1]], 0.0);
         }
      } else { // labeling 111 -> all reparametrized values <= 0
         if(sortIndices[0] == i_) {
            msg[0] -= rightPot[i_];
         } else if(sortIndices[1] == i_) {
            msg[0] -= rightPot[i_];
         } else {
            msg[0] -= rightPot[i_];
         }
      }
   }

   template<typename LEFT_FACTOR, typename RIGHT_FACTOR, typename G1, typename G2, typename G3, MessageSending MST_TMP = MST>
   typename std::enable_if<MST_TMP == MessageSending::SRMP,void>::type
   SendMessageToRight(LEFT_FACTOR* const l, RIGHT_FACTOR* const r, const G1& leftPot, const G2& rightPot, G3& msg, const REAL omega)
   {
      static_assert(MST_TMP == MST,"");
      msg[0] += omega*leftPot[0];
   }

   template<typename LEFT_FACTOR, typename G1, typename G2, MessageSending MST_TMP = MST>
   typename std::enable_if<MST_TMP == MessageSending::MPLP,void>::type
   ReceiveMessageFromLeft(LEFT_FACTOR* const l, const G1& leftPot, G2& msg) 
   {
      static_assert(MST_TMP == MST,"");
      msg[0] += leftPot[0];
   }

   template<typename RIGHT_FACTOR, typename MSG_ARRAY, typename RIGHT_REPAM, typename ITERATOR, MessageSending MST_TMP = MST>
   static typename std::enable_if<MST_TMP == MessageSending::MPLP,void>::type
   SendMessagesToLeft(const RIGHT_FACTOR& rightFactor, const RIGHT_REPAM& rightRepam, const MSG_ARRAY& msgs, ITERATOR omegaIt)
   {
      static_assert(MST_TMP == MST,"");
      const REAL omega = std::accumulate(omegaIt, omegaIt+3,0.0);
      //const REAL omega = 1.0; // do zrobienia: for now, in general this will not converge
      assert(omega < 1.0 + eps);
      //std::cout << "(" << *omegaIt << "," << *(omegaIt+1) << "," << *(omegaIt+2) << ")\n";
      const auto sortIndices = rightFactor.SortIndices(rightRepam);
      const REAL x = rightRepam[sortIndices[0]] + rightRepam[sortIndices[1]];

      assert(msgs.size() == 3);
      // do zrobienia: stupid interface. Make references. See also in factors_messages.hxx
      assert(msgs[0]->GetMessageOp().i_ == 0);
      assert(msgs[1]->GetMessageOp().i_ == 1);
      assert(msgs[2]->GetMessageOp().i_ == 2);

      if(x > 0.0) { // labeling 000
         if(rightRepam[sortIndices[0]] < 0.0) {
            //msgs[sortIndices[2]]->operator[](0) += omega*(-rightRepam[sortIndices[2]] + rightRepam[sortIndices[1]] - x);
            //msgs[sortIndices[1]]->operator[](0) -= omega*x;
            const REAL delta = rightRepam[sortIndices[1]] - 0.5*x;
            msgs[sortIndices[0]]->operator[](0) -= 0.5*omega*x;
            msgs[sortIndices[1]]->operator[](0) -= 0.5*omega*x;
            msgs[sortIndices[2]]->operator[](0) += omega*(-rightRepam[sortIndices[2]] + delta);
         } else {
            for(INDEX i=0; i<3; ++i) {
               msgs[i]->operator[](0) -= omega*rightRepam[i];
            }
         }
      } else if(rightRepam[sortIndices[2]] > 0.0) { //labeling 110
         if(rightRepam[sortIndices[1]] > 0.0) {
            //msgs[sortIndices[0]]->operator[](0) -= omega*x;
            //msgs[sortIndices[2]]->operator[](0) += omega*(-rightRepam[sortIndices[2]] + rightRepam[sortIndices[1]]);
            // do zrobiebia: make decrease of rightRepam[sortIndices[2]] and increase of rightRepam[sortIndices[1]] as equal as pssible.
            const REAL b = std::min(-0.5*x,0.5*(rightRepam[sortIndices[2]] - rightRepam[sortIndices[1]]));
            const REAL delta = rightRepam[sortIndices[1]] + b;
            msgs[sortIndices[1]]->operator[](0) += omega*b;
            msgs[sortIndices[2]]->operator[](0) += omega*(-rightRepam[sortIndices[2]] + delta);
            msgs[sortIndices[0]]->operator[](0) -= omega*(0.5*x + b); // do zrobienia: can possibly be heightened more
         } else {
            for(INDEX i=0; i<3; ++i) {
               msgs[i]->operator[](0) -= omega*rightRepam[i];
            }
         }
      } else { // labeling 111
         for(INDEX i=0; i<3; ++i) {
            msgs[i]->operator[](0) -= omega*rightRepam[i];
         }
      }
   }

   template<typename G>
   void RepamLeft(G& repamPot, const REAL msg, const INDEX msg_dim)
   {
      assert(msg_dim == 0);
      repamPot[0] -= msg;
   }
   template<typename G>
   void RepamRight(G& repamPot, const REAL msg, const INDEX msg_dim)
   {
      assert(msg_dim == 0);
      repamPot[i_] += msg;
   }

   // compute primal functions, how to do it? look into ultra-metric rounding
   void ComputeRightFromLeftPrimal(const bool leftPrimal, MulticutTripletFactor::LabelingType& rightPrimal)
   {
      rightPrimal[i_] = leftPrimal;
   }


private:
   const INDEX i_; // index of the affected variable in the cycle factor, do zrobienia: needs only two bits
};

} // end namespace LP_MP

#endif // LP_MP_MULTICUT_MESSAGE_HXX

