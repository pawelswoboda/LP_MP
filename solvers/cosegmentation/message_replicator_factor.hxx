#ifndef LP_MP_MESSAGE_REPLICATOR_FACTOR_HXX
#define LP_MP_MESSAGE_REPLICATOR_FACTOR_HXX

#include "instances.inc"

namespace LP_MP {

// factor for circumventing many factor message: replicate messages and require them to be equal: 
// Perform this in every ReceiveMessageFrom{Left|Right}-step.

// make MessageReplicator a friend
class MessageReplicatorFactor {
   template<Chirality ChiralityType, class BaseMessageClass>
   friend class MessageReplicator;
public:
   MessageReplicatorFactor(const INDEX MessageSize) : message_(MessageSize,0.0) {}

   void MaximizePotential() {}
   template<typename REPAM_POT>
   REAL LowerBound(const REPAM_POT& repamPot) const {
      return 0.0;
   }

   REAL operator[](const INDEX i) const { return 0.0; }
   INDEX size() const { return message_.size(); }

   REAL GetTrueMessage(const INDEX i) { return message_[i]; }

private:
   std::vector<REAL> message_; // the true value of the message

};

} // end namespace LP_MP

#endif // LP_MP_MESSAGE_REPLICATOR_FACTOR_HXX
