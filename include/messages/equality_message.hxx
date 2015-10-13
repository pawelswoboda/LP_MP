#ifndef LP_MP_EQUALITY_MESSAGE
#define LP_MP_EQUALITY_MESSAGE

namespace LP_MP {

// maximize/minimize to second min/max
// assume FactorType is Simplex. 
// do zrobienia: or multiplex
// do zrobienia: use breakpoINDEX cost for message updates
class EqualityMessage 
{
public:
   EqualityMessage(const INDEX lV, const INDEX rV)
      : leftVar_(lV), rightVar_(rV)
   {}

   template<typename G1, typename G2>
      void MakeRightFactorUniform(const G1& rightPot, G2& msg, const REAL omega = 1.0);
   template<typename G1, typename G2>
      void MakeLeftFactorUniform(const G1& leftPot, G2& msg, const REAL omega = 1.0);

   template<typename OP, typename REPAM_ARRAY, typename MSG>
   void MakeFactorUniform(const OP op, const REPAM_ARRAY& repamPot, MSG& msg, const INDEX var_idx, const REAL omega = 1.0)
   {
      assert(msg.size() == 1);
      assert(var_idx < repamPot.size());

      // possibly do it differently: search for two second smallest entries and then select first or second one depending upon whether it is rightVar_ or not. Faster?
      REAL min_val = std::numeric_limits<REAL>::max();
      for(INDEX i=0; i<repamPot.size(); ++i) {
         if(i!=var_idx) {
            min_val = std::min(min_val, repamPot[i]);
         }
      }
      // could possibly be replaced by, but repamPot need not support iterators
      //min_val = *std::min_element(repamPot.cbegin(), repamPot.cbegin()+var_idx);
      //min_val = std::min(*std::min_element(repamPot.cbegin() + var_idx+1, repamPot.cend()), min_val);
      // faster?

      const REAL new_msg = msg[0] + op(min_val - repamPot[var_idx]);
      const REAL old_msg = msg[0];
      msg[0] = (1.0-omega)*old_msg + omega*new_msg;
   }

   template<typename RIGHT_FACTOR, typename G1, typename G2>
   void ReceiveMessageFromRight(RIGHT_FACTOR* const r, const G1& rightPot, G2& msg)
   { 
      auto op = [](const REAL x) { return +x; };
      MakeFactorUniform<decltype(op)>(op,rightPot, msg, rightVar_);
      //MakeRightFactorUniform(rightPot, msg); 
   }

   template<typename LEFT_FACTOR, typename G1, typename G2>
   void ReceiveMessageFromLeft(LEFT_FACTOR* l, const G1& leftPot, G2& msg)
   { 
      auto op = [](const REAL x) { return -x; };
      MakeFactorUniform<decltype(op)>(op,leftPot, msg, leftVar_);
      //MakeLeftFactorUniform(leftPot, msg); 
   }

   template<typename LEFT_FACTOR, typename RIGHT_FACTOR, typename G1, typename G2, typename G3>
   void SendMessageToRight(LEFT_FACTOR* const l, RIGHT_FACTOR* const r, const G1& leftPot, const G2& rightPot, G3& msg, const REAL omega)
   { //std::cout << "Send message to right equal\n";
      auto op = [](const REAL x) { return -x; };
      MakeFactorUniform<decltype(op)>(op,leftPot, msg, leftVar_, omega);
      //MakeLeftFactorUniform(leftPot,msg,omega); 
   }
   
   template<typename LEFT_FACTOR, typename RIGHT_FACTOR, typename G1, typename G2, typename G3>
   void SendMessageToLeft(LEFT_FACTOR* const l, RIGHT_FACTOR* const r, const G1& leftPot, const G2& rightPot, G3& msg, const REAL omega)
   { 
      auto op = [](const REAL x) { return +x; };
      MakeFactorUniform<decltype(op)>(op,rightPot, msg, rightVar_, omega);
      //MakeRightFactorUniform(rightPot,msg,omega); 
   }

   // send all messages of the same type at once
   // possibly templatize this so that case, where all variables of left factor are accessed exatly once, is processed with greater speed (only one minimum searching required then)
   // template specialization is also possible for every dimension accessed except for last one
   // assume that all attached factors have different endpoINDEXs

   // do zrobienia: not best encapsulation: we have access to MessageTypeCRTP holding EqualityMessage, not to the EqualityMessage
   // idea: construct proxy object that will give access directly to EqualityMessage, but has no need to construct  the array explicitly
   //
   // code duplication: templatize code

   // for sending multiple messages at once: makes factor uniform by sending all messages at once
   template<typename SIGN_OP, typename VAR_ACCESS_OP, typename MSG_ARRAY, typename RIGHT_REPAM, typename ITERATOR>
   static void MakeFactorUniformParallel(const SIGN_OP sign_op, VAR_ACCESS_OP var_access_op, const MSG_ARRAY& msgs, const RIGHT_REPAM& repam, ITERATOR omegaIt)
   {
      assert(msgs.size() >= 2); // otherwise calling this method makes no sense
      assert(msgs.size() <= repam.size());
      const ITERATOR omegaEnd = omegaIt + msgs.size();

      const REAL omega_sum = 0.7 * std::accumulate(omegaIt, omegaEnd, 0.0); // strangely, a smaller factor makes the algortihm faster
      assert(omega_sum <= 1.0 + 1e-6);

      // find minimal value of potential over all indices accessed by messages
      REAL min_val_covered = std::numeric_limits<REAL>::max();
      for(INDEX msg_idx=0; msg_idx<msgs.size(); ++msg_idx) {
         const INDEX var_idx = var_access_op(msgs[msg_idx]->msg_op_);
         //std::cout << "leftVar = " << leftVar << "\n";
         min_val_covered = std::min(min_val_covered, repam[var_idx]);
      }

      REAL min_val = std::numeric_limits<REAL>::max();
      REAL second_min_val = std::numeric_limits<REAL>::max();
      for(INDEX i=0; i<repam.size(); ++i) {
         const REAL cur_val = repam[i];
         //std::cout << "cur_val = " << cur_val << "\n";
         if(min_val >= cur_val) { // if two values are equally small, second_min_val should be the lowest value, too
            second_min_val = min_val;
            min_val = cur_val;
         } else if(second_min_val > cur_val) {
            second_min_val = cur_val;
         }
      }
      assert(std::make_pair(min_val, second_min_val) == SmallestValues<REAL>(repam));
      //assert(second_min_val != std::numeric_limits<REAL>::max());

      REAL new_val;  // this value will be taken by the new reparametrized entries
      if(min_val < min_val_covered) { new_val = min_val; }
      else { new_val = second_min_val; }

      //std::cout << "omega_sum = " << omega_sum << ", no messages = " << msgs.size() << ", potential size = " << leftRepam.size() << "\n";
      //std::cout << "min_val_covered = " << min_val_covered << ", min_val = " << min_val << ", second_min_val = " << second_min_val << ", new_val = " << new_val << "\n";
      //const INDEX last_idx = leftRepam.size() - 1;
      //std::cout << "not covered = " << leftRepam[last_idx] << "\n";

      for(INDEX msg_idx=0; msg_idx<msgs.size(); ++msg_idx, omegaIt++) {
         if(*omegaIt > 0) {
         const INDEX var_idx = var_access_op(msgs[msg_idx]->msg_op_);
         const REAL new_msg = msgs[msg_idx]->operator[](0) + sign_op(new_val - repam[var_idx]);
         const REAL old_msg = msgs[msg_idx]->operator[](0);
         //std::cout << "new message = " << (1.0-omega_sum)*old_msg + omega_sum*new_msg << ", delta = " << -(new_val - leftRepam[leftVar]) << "\n";
         //msgs[msg_idx]->operator[](0) = (1.0-*omegaBegin)*old_msg + *omegaBegin*new_msg;
         msgs[msg_idx]->operator[](0) = (1.0-omega_sum)*old_msg + omega_sum*new_msg;
         }
      }
   }

   template<typename MSG_ARRAY, typename RIGHT_REPAM, typename ITERATOR>
   static void SendMessagesToLeft(const MSG_ARRAY& msgs, const RIGHT_REPAM& rightRepam, ITERATOR omegaIt)
   {
      auto sign_op = [](const REAL x) -> REAL { return +x; };
      auto var_access_op = [](const EqualityMessage& msg) -> INDEX { return msg.rightVar_; };
      MakeFactorUniformParallel(sign_op, var_access_op, msgs, rightRepam, omegaIt);
   }

   template<typename MSG_ARRAY, typename LEFT_REPAM, typename ITERATOR>
   static void SendMessagesToRight(const MSG_ARRAY& msgs, const LEFT_REPAM& leftRepam, ITERATOR omegaIt)
   {
      auto sign_op = [](const REAL x) -> REAL { return -x; };
      auto var_access_op = [](const EqualityMessage& msg) -> INDEX { return msg.leftVar_; };
      MakeFactorUniformParallel(sign_op, var_access_op, msgs, leftRepam, omegaIt);
   }

   template<typename G>
   void RepamLeft(G& leftRepamPot, const REAL msg, const INDEX dim) { 
      assert(dim == 0); 
      leftRepamPot[leftVar_] = leftRepamPot[leftVar_] - msg; 
   }
   template<typename G>
   void RepamRight(G& rightRepamPot, const REAL msg, const INDEX dim) { 
      assert(dim == 0); 
      rightRepamPot[rightVar_] = rightRepamPot[rightVar_] + msg; 
   }

   // for implicit repam storage, not really needed, only test
   template<typename M>
   const REAL GetLeftMessage(const INDEX i, const M& msg) const
   {
      if(i==leftVar_) return -msg[0];
      else return 0.0;
   }
   template<typename M>
   const REAL GetRightMessage(const INDEX i, const M& msg) const
   {
      if(i==rightVar_) return msg[0];
      else return 0.0;
   }
   

private:
   //do zrobienia: possibly SHORT_INDEX
   const INDEX leftVar_, rightVar_; // variables affected 
};

/*
template<typename G1, typename G2>
void EqualityMessage::MakeRightFactorUniform(const G1& rightPot, G2& msg, const REAL omega)
{
   assert(msg.size() == 1);
   assert(rightVar_ < rightPot.size());
   
   // possibly do it differently: search for two second smallest entries and then select first or second one depending upon whether it is rightVar_ or not. Faster?
   REAL min_val = std::numeric_limits<REAL>::max();
   for(INDEX i=0; i<rightPot.size(); ++i) {
      if(i!=rightVar_) {
         min_val = std::min(min_val, rightPot[i]);
      }
   }
   // can be replaced by
   //min_val = *std::min_element(rightPot.cbegin(), rightPot.cbegin()+rightVar_);
   //min_val = std::min(*std::min_element(rightPot.cbegin() + rightVar+1, rightPot.cend()), min_val);
   // faster?

   REAL omega_ = std::min(15.0*omega,1.0);
   REAL new_msg = msg[0] + (min_val - rightPot[rightVar_]);
   REAL old_msg = msg[0];
   msg[0] = (1.0-omega_)*old_msg + omega_*new_msg;
}

template<typename G1, typename G2>
void EqualityMessage::MakeLeftFactorUniform(const G1& leftPot, G2& msg, const REAL omega)
{
   assert(msg.size() == 1);
   assert(leftVar_ < leftPot.size());
   
   REAL min_val = std::numeric_limits<REAL>::max();
   for(INDEX i=0; i<leftPot.size(); ++i) {
      if(i!=leftVar_) {
         min_val = std::min(min_val, leftPot[i]);
      }
   }

   REAL omega_ = std::min(15.0*omega,1.0);
   REAL new_msg = msg[0] - (min_val - leftPot[leftVar_]);
   REAL old_msg = msg[0];
   msg[0] = (1.0-omega_)*old_msg + omega_*new_msg;
}
*/

} // end namespace LP_MP

#endif // LP_MP_EQUALITY_MESSAGE
