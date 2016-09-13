#ifndef LP_MP_EQUALITY_MESSAGE
#define LP_MP_EQUALITY_MESSAGE

#include "config.hxx"

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

      msg[0] -= omega*(repamPot[var_idx] - min_val);
      return;
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

   template<typename RIGHT_FACTOR, typename G1, typename G2>
   void ReceiveRestrictedMessageFromRight(RIGHT_FACTOR* const r, const G1& rightPot, G2& msg, typename PrimalSolutionStorage::Element rightPrimal)
   {
      if(rightPrimal[rightVar_] == false) {
         msg[0] -= std::numeric_limits<REAL>::infinity();
      } else if(rightPrimal[rightVar_] == true) {
         msg[0] += std::numeric_limits<REAL>::infinity();
      }
   }

   template<typename LEFT_FACTOR, typename G1, typename G2>
   void ReceiveRestrictedMessageFromLeft(LEFT_FACTOR* l, const G1& leftPot, G2& msg, typename PrimalSolutionStorage::Element leftPrimal)
   { 
      if(leftPrimal[leftVar_] == false) {
         msg[0] -= std::numeric_limits<REAL>::infinity();
      } else if(leftPrimal[leftVar_] == true) {
         msg[0] += std::numeric_limits<REAL>::infinity();
      }
   }
   /*
   template<typename LEFT_FACTOR, typename G1, typename G3>
   void SendMessageToRight(LEFT_FACTOR* const l, const G1& leftPot, G3& msg, const REAL omega)
   { //std::cout << "Send message to right equal\n";
      auto op = [](const REAL x) { return -x; };
      MakeFactorUniform<decltype(op)>(op,leftPot, msg, leftVar_, omega);
      //MakeLeftFactorUniform(leftPot,msg,omega); 
   }
   
   template<typename RIGHT_FACTOR, typename G2, typename G3>
   void SendMessageToLeft(RIGHT_FACTOR* const r, const G2& rightPot, G3& msg, const REAL omega)
   { 
      auto op = [](const REAL x) { return +x; };
      MakeFactorUniform<decltype(op)>(op,rightPot, msg, rightVar_, omega);
      //MakeRightFactorUniform(rightPot,msg,omega); 
   }
   */

   // send all messages of the same type at once
   // possibly templatize this so that case, where all variables of left factor are accessed exatly once, is processed with greater speed (only one minimum searching required then)
   // template specialization is also possible for every dimension accessed except for last one
   // assume that all attached factors have different endpoINDEXs

   // do zrobienia: not best encapsulation: we have access to MessageTypeCRTP holding EqualityMessage, not to the EqualityMessage
   // idea: construct proxy object that will give access directly to EqualityMessage, but has no need to construct  the array explicitly
   //
   // code duplication: templatize code

   // for sending multiple messages at once: makes factor uniform by sending all messages at once
   template<typename VAR_ACCESS_OP, typename MSG_ARRAY, typename RIGHT_REPAM, typename ITERATOR>
   static void MakeFactorUniformParallel(VAR_ACCESS_OP var_access_op, const MSG_ARRAY& msgs, const RIGHT_REPAM& repam, ITERATOR omegaIt)
   {
      //assert(msgs.size() >= 2); // otherwise calling this method makes no sense, but it can happen for some trivial problems.
      assert(msgs.size() <= repam.size());
      //assert(msgs.size() == repam.size()-1); // special case of one edge is not assignment in QAP or cosegmentation. For now. Only for hotel and house
      const ITERATOR omegaEnd = omegaIt + msgs.size();

      // do zrobienia:
      const REAL omega_sum = 0.5 * std::accumulate(omegaIt, omegaEnd, 0.0); // strangely, a smaller factor makes the algorithm faster
      //const REAL omega_sum = std::accumulate(omegaIt, omegaEnd, 0.0);
      assert(omega_sum <= 1.0 + eps);

      // find minimal value of potential over all indices accessed by messages
      REAL min_val_covered = std::numeric_limits<REAL>::max();
      for(INDEX msg_idx=0; msg_idx<msgs.size(); ++msg_idx) {
         const INDEX var_idx = var_access_op(msgs[msg_idx].GetMessageOp());
         //assert(var_idx != repam.size()-1); // this is only valied for assignment problems from house and hotel
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

      REAL new_val;  // this value will be taken by the new reparametrized entries
      if(min_val < min_val_covered) { new_val = min_val; }
      else { new_val = second_min_val; }

      //std::cout << "omega_sum = " << omega_sum << ", no messages = " << msgs.size() << ", potential size = " << leftRepam.size() << "\n";
      //std::cout << "min_val_covered = " << min_val_covered << ", min_val = " << min_val << ", second_min_val = " << second_min_val << ", new_val = " << new_val << "\n";
      //const INDEX last_idx = leftRepam.size() - 1;
      //std::cout << "not covered = " << leftRepam[last_idx] << "\n";

      for(INDEX msg_idx=0; msg_idx<msgs.size(); ++msg_idx, omegaIt++) {
         if(*omegaIt > 0) {
            const INDEX var_idx = var_access_op(msgs[msg_idx].GetMessageOp());
            msgs[msg_idx].operator[](0) -= omega_sum*(repam[var_idx] - new_val);
         }
      }
   }

   // do zrobienia: enable again
   template<typename RIGHT_FACTOR, typename MSG_ARRAY, typename RIGHT_REPAM, typename ITERATOR>
   static void SendMessagesToLeft(const RIGHT_FACTOR& rightFactor, const RIGHT_REPAM& rightRepam, const MSG_ARRAY& msgs, ITERATOR omegaIt)
   {
      auto var_access_op = [](const EqualityMessage& msg) -> INDEX { return msg.rightVar_; };
      MakeFactorUniformParallel(var_access_op, msgs, rightRepam, omegaIt);
   }

   template<typename LEFT_FACTOR, typename MSG_ARRAY, typename LEFT_REPAM, typename ITERATOR>
   static void SendMessagesToRight(const LEFT_FACTOR& leftFactor, const LEFT_REPAM& leftRepam, const MSG_ARRAY& msgs, ITERATOR omegaIt)
   {
      auto var_access_op = [](const EqualityMessage& msg) -> INDEX { return msg.leftVar_; };
      MakeFactorUniformParallel(var_access_op, msgs, leftRepam, omegaIt);
   }

   template<typename G>
   void RepamLeft(G& leftRepamPot, const REAL msg, const INDEX dim) { 
      assert(dim == 0); 
      leftRepamPot[leftVar_] += msg; 
   }
   template<typename G>
   void RepamRight(G& rightRepamPot, const REAL msg, const INDEX dim) { 
      assert(dim == 0); 
      rightRepamPot[rightVar_] += msg; 
   }

   // for implicit repam storage, not really needed, only test
   template<typename M>
   const REAL GetLeftMessage(const INDEX i, const M& msg) const
   {
      if(i==leftVar_) return msg[0];
      else return 0.0;
   }
   template<typename M>
   const REAL GetRightMessage(const INDEX i, const M& msg) const
   {
      if(i==rightVar_) return msg[0];
      else return 0.0;
   }

   /*
   void ComputeLeftFromRightPrimal(PrimalSolutionStorage::Element left, PrimalSolutionStorage::Element right) 
   {
      if(right[rightVar_] == true) { 
         left[leftVar_] = true;
         // it would be nice to set all other entries to false
      } else if(right[rightVar_] == false) {
         left[leftVar_] = false;
      }
   }
   */

   void ComputeRightFromLeftPrimal(PrimalSolutionStorage::Element left, PrimalSolutionStorage::Element right)
   {
      if(left[leftVar_] == true) { 
         right[rightVar_] = true;
         // it would be nice to set all other entries to false
      } else if(left[leftVar_] == false) {
         right[rightVar_] = false;
      }
   }
  
   // here it is checked whether labeling on left side and labeling on right side fulfill the constraints of the message
   // note: If we build an LP-model, this could be checked automatically!
   bool CheckPrimalConsistency(PrimalSolutionStorage::Element leftPrimal, PrimalSolutionStorage::Element rightPrimal) const
   {
      if(leftPrimal[leftVar_]) { return rightPrimal[rightVar_] == true; }
      if(rightPrimal[rightVar_]) { return leftPrimal[leftVar_] == true; }
      return true;
   }


private:
   //do zrobienia: possibly SHORT_INDEX or some 16 bit index (i.e. short unsigned int)
   const INDEX leftVar_, rightVar_; // variables affected 
};

} // end namespace LP_MP

#endif // LP_MP_EQUALITY_MESSAGE
