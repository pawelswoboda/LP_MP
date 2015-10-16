#ifndef LP_MP_EMPTY_MESSAGE_TEMPLATE_HXX
#define LP_MP_EMPTY_MESSAGE_TEMPLATE_HXX

#include "instances.inc"

namespace LP_MP {

class EmptyMessageTemplate {
public:
   template<typename RIGHT_FACTOR, typename G1, typename G2>
   void ReceiveMessageFromRight(RIGHT_FACTOR* const r, const G1& rightPot, G2& msg)
   {}

   template<typename LEFT_FACTOR, typename G1, typename G2>
   void ReceiveMessageFromLeft(LEFT_FACTOR* l, const G1& leftPot, G2& msg)
   {}

   template<typename LEFT_FACTOR, typename RIGHT_FACTOR, typename G1, typename G2, typename G3>
   void SendMessageToRight(LEFT_FACTOR* const l, RIGHT_FACTOR* const r, const G1& leftPot, const G2& rightPot, G3& msg, const REAL omega)
   {}

   template<typename LEFT_FACTOR, typename RIGHT_FACTOR, typename G1, typename G2, typename G3>
   void SendMessageToLeft(LEFT_FACTOR* const l, RIGHT_FACTOR* const r, const G1& leftPot, const G2& rightPot, G3& msg, const REAL omega)
   {}

   template<typename MSG_ARRAY, typename RIGHT_REPAM, typename ITERATOR>
   static void SendMessagesToLeft(const MSG_ARRAY& msgs, const RIGHT_REPAM& rightRepam, ITERATOR omegaIt)
   {}

   template<typename MSG_ARRAY, typename LEFT_REPAM, typename ITERATOR>
   static void SendMessagesToRight(const MSG_ARRAY& msgs, const LEFT_REPAM& leftRepam, ITERATOR omegaIt)
   {}

   template<typename G>
   void RepamLeft(G& leftRepamPot, const REAL msg, const INDEX dim)
   {}
   template<typename G>
   void RepamRight(G& rightRepamPot, const REAL msg, const INDEX dim)
   {}

   template<typename M>
   const REAL GetLeftMessage(const INDEX i, const M& msg) const
   {}
   template<typename M>
   const REAL GetRightMessage(const INDEX i, const M& msg) const
   {}
};

} // end namespace LP_MP

#endif // LP_MP_EMPTY_MESSAGE_TEMPLATE_HXX

