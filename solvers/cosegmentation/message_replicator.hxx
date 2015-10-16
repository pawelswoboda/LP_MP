#ifndef LP_MP_MESSAGE_REPLICATOR_HXX
#define LP_MP_MESSAGE_REPLICATOR_HXX

#include "instances.inc"
#include <type_traits>

namespace LP_MP {

// MessageReplicator is used to replicate messages, so that in certain circumstances, more than two factors can be joined via a single message. 
// We defer reparametrization updates, so that it only happens, when ReceiveMessageFrom* is called

// SFINAE seems not to work with inheritance, implement explicit SFINAE switches by defining every function twice and invoking inherited function explicitly, as done with RepamLeft

template<Chirality ChiralityType, class BaseMessageClass> // Chirality can be left or right, indicating whether true message is on left or right side
class MessageReplicator : public BaseMessageClass {
public:
   // for this to work, ReceiveMessageFrom* must be called every time in factor_messages.hxx
   template<typename RIGHT_FACTOR, typename G1, typename G2, Chirality CT = ChiralityType>
   typename std::enable_if<CT == Chirality::left>::type
   ReceiveMessageFromRight(RIGHT_FACTOR* const r, const G1& rightPot, G2& msg)
   {
      //static_assert(RIGHT_FACTOR == MessageReplicatorFactor, "");
      //msg = r->message_;
   }

   template<typename LEFT_FACTOR, typename G1, typename G2, Chirality CT = ChiralityType>
   typename std::enable_if<CT == Chirality::right>::type
   ReceiveMessageFromLeft(LEFT_FACTOR* l, const G1& leftPot, G2& msg)
   {
      //msg = l->message_;
   }

   template<typename LEFT_FACTOR, typename RIGHT_FACTOR, typename G1, typename G2, typename G3, Chirality CT = ChiralityType>
   typename std::enable_if<CT == Chirality::right>::type
   SendMessageToRight(LEFT_FACTOR* const l, RIGHT_FACTOR* const r, const G1& leftPot, const G2& rightPot, G3& msg, const REAL omega)
   {
      //msg = r->message_;
   }

   template<typename LEFT_FACTOR, typename RIGHT_FACTOR, typename G1, typename G2, typename G3, Chirality CT = ChiralityType>
   typename std::enable_if<CT == Chirality::right>::type
   SendMessageToLeft(LEFT_FACTOR* const l, RIGHT_FACTOR* const r, const G1& leftPot, const G2& rightPot, G3& msg, const REAL omega)
   {
      //msg = l->message_;
   }

   template<typename G, Chirality CT = ChiralityType>
   typename std::enable_if<CT == Chirality::right>::type
   RepamLeft(G& leftRepamPot, const REAL msg, const INDEX dim)
   {}
   template<typename G, Chirality CT = ChiralityType>
   typename std::enable_if<CT == Chirality::left>::type
   RepamLeft(G& leftRepamPot, const REAL msg, const INDEX dim)
   {
      BaseMessageClass::RepamLeft(leftRepamPot,msg,dim);
   }

      /*
   template<typename G, Chirality CT = ChiralityType>
   typename std::enable_if<CT == Chirality::left>::type
   RepamRight(G& rightRepamPot, const REAL msg, const INDEX dim)
   {}
   */

   template<typename M, Chirality CT = ChiralityType>
   typename std::enable_if<CT == Chirality::right>::type
   GetLeftMessage(const INDEX i, const M& msg) const
   { return 0.0; }

   template<typename M, Chirality CT = ChiralityType>
   typename std::enable_if<CT == Chirality::left>::type
   GetRightMessage(const INDEX i, const M& msg) const
   { return 0.0; }
};

} // end namespace LP_MP

#endif // LP_MP_MESSAGE_REPLICATOR_HXX
