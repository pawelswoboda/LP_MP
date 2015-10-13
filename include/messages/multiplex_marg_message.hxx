#ifndef LP_MP_MULTIPLEX_MARG_MESSAGE
#define LP_MP_MULTIPLEX_MARG_MESSAGE

#include "LP_MP.h"
#include "factors_messages.hxx"
#include "message_loops.hxx"

namespace LP_MP {


// kL*(l_11 + ... + l_1n) = kR*(r_11 + ... + r_l1)
// ...
// kL*(l_k1 + ... + l_kn) = kR*(r_1k + ... + r_lk)
// all variables must be distinct
// left and right factor must be either simplex factor or multiplex factor

// template parameters: left and right shape iterator, which messages to compute, i.e. message for pairwise factors need not be computed
// note: left and right factor type must be MultiplexFactor
// do zrobienia: make functions const where applicable
// {LEFT|RIGHT}_SIDE_ACTIVE activates with SFINAE message updates for left and right factor resp.
template<class LEFT_LOOP_TYPE, class RIGHT_LOOP_TYPE, bool LEFT_SIDE_ACTIVE, bool RIGHT_SIDE_ACTIVE>
class MultiplexMargMessage
{
public:
   typedef LEFT_LOOP_TYPE LeftLoopType;
   typedef RIGHT_LOOP_TYPE RightLoopType;

   MultiplexMargMessage(LEFT_LOOP_TYPE loopLeft, RIGHT_LOOP_TYPE loopRight, const INDEX kL, const INDEX kR)
      : loopLeft_(loopLeft),
      loopRight_(loopRight),
      kL_(kL),
      kR_(kR)
   { 
      assert(kL_ > 0 && kR_ > 0);
      assert(kL_ == 1 && kR_ == 1); // do zrobienia: for now
   }

   // standard functions which take all possible arguments, to be replaced with minimal ones
   template<typename RIGHT_FACTOR, typename G1, typename G2, bool LSA = LEFT_SIDE_ACTIVE>
   typename std::enable_if<LSA,void>::type
   ReceiveMessageFromRight(RIGHT_FACTOR* const r, const G1& rightPot, G2& msg) 
   {
      //std::cout << "LEFT_SIDE_ACTIVE " << LEFT_SIDE_ACTIVE << "\n"; // jest true
      MinimizeRight(r,rightPot,msg); 
   }
   template<typename LEFT_FACTOR, typename G1, typename G2, bool RSA = RIGHT_SIDE_ACTIVE>
   typename std::enable_if<RSA,void>::type 
   ReceiveMessageFromLeft(LEFT_FACTOR* l, const G1& leftPot, G2& msg) 
   { 
      //std::cout << "RIGHT_SIDE_ACTIVE " << RIGHT_SIDE_ACTIVE << "\n"; // jest false
      MaximizeLeft(l,leftPot,msg); 
   }
   template<typename LEFT_FACTOR, typename RIGHT_FACTOR, typename G1, typename G2, typename G3, bool LSA = LEFT_SIDE_ACTIVE>
   typename std::enable_if<LSA,void>::type
   SendMessageToRight(LEFT_FACTOR* const l, RIGHT_FACTOR* const r, const G1& leftPot, const G2& rightPot, G3& msg, const REAL omega)
   { 
      //std::cout << "LEFT_SIDE_ACTIVE " << LEFT_SIDE_ACTIVE << "\n"; // jest true
      MaximizeLeft(l,leftPot,msg,omega); 
   }
   template<typename LEFT_FACTOR, typename RIGHT_FACTOR, typename G1, typename G2, typename G3, bool RSA = RIGHT_SIDE_ACTIVE>
   typename std::enable_if<RSA,void>::type
   SendMessageToLeft(LEFT_FACTOR* const l, RIGHT_FACTOR* const r, const G1& leftPot, const G2& rightPot, G3& msg, const REAL omega) 
   { 
      //std::cout << "RIGHT_SIDE_ACTIVE " << RIGHT_SIDE_ACTIVE << "\n"; // jest false
      MinimizeRight(r,rightPot,msg,omega); 
   }

   /*
   // do zrobienia: only templates, fill them
   template<typename MSG_ARRAY, typename RIGHT_REPAM, typename ITERATOR>
   static void
   SendMessagesToLeft(const typename std::enable_if<RIGHT_SIDE_ACTIVE, MSG_ARRAY>::type & msgs, const RIGHT_REPAM& leftRepam, ITERATOR omegaIt)
   {
      std::cout << "kwasikus\n";
      if(RIGHT_SIDE_ACTIVE) {
         
      }
   }

   template<typename MSG_ARRAY, typename LEFT_REPAM, typename ITERATOR>
   static void
   SendMessagesToRight(const typename std::enable_if<LEFT_SIDE_ACTIVE, MSG_ARRAY&>::type msgs, const LEFT_REPAM& leftRepam, ITERATOR omegaIt)
   { 
      std::cout << "kwasikus to right\n";
      if(LEFT_SIDE_ACTIVE) {

      }
   }
   */

   // given multiple message of the same type, process them all at once
   // do zrobienia: SendMessageToLeft
   /*
   template<typename MSG_ARRAY, typename REPAM_ARRAY, typename ITERATOR, bool LSA = LEFT_SIDE_ACTIVE>
   static typename std::enable_if<LSA,void>::type 
   SendMessagesToRight(const MSG_ARRAY& msgs, const REPAM_ARRAY& leftRepam, ITERATOR omegaIt) 
   {
      REAL delta;
      auto pot = static_cast<typename decltype(*msgs[0])::LEFT_FACTOR_TYPE*>(msgs[0]->GetLeftFactor());
      auto leftLoop = msgs[0]->msg_op_.loopLeft_;
      const REAL breakPoINDEX = msgs[0]->GetLeftFactor()->GetBreakpoINDEXCost(pot);
      const REAL k = msgs[0]->msg_op_.kL_;
      leftLoop.loop( 
         [&](const INDEX outer_idx){ delta = std::numeric_limits<REAL>::max(); 
         //std::cout << "outer_idx = " << outer_idx << "\n";
         }, 
         [&](const INDEX full_idx, const INDEX outer_idx){ 
         //std::cout << "full_idx = " << full_idx <<  ", outer_idx = " << outer_idx << "\n";
         assert(full_idx < pot.size());
         if(pot[ full_idx ] > breakPoINDEX) { delta = std::min(delta, (pot[ full_idx ] - breakPoINDEX)/REAL(k)); }
         else { delta = 0.0; }
         },
         [&](const INDEX outer_idx){ 
         //std::cout << "outer_idx = " << outer_idx << "\n";
         for(INDEX msg_idx=0; msg_idx<msgs.size(); ++msg_idx, ++omegaIt) {
            REAL new_msg = msgs[msg_idx][ outer_idx ] + delta;
            REAL old_msg = msgs[msg_idx][ outer_idx ] + 0.0;
            msgs[msg_idx][ outer_idx ] = (1.0-*omegaIt[msg_idx])*old_msg + *omegaIt[msg_idx]*new_msg;
         }
         } );
   }
   */
   // reparametrize left potential for i-th entry of msg
   template<typename G>
   void RepamLeft(G& repamPot, const REAL msg, const INDEX dim)
   { 
      loopLeft_.loop(dim, [&](const INDEX i) { repamPot[i] = repamPot[i] - msg; });
   }
   template<typename G>
   void RepamRight(G& repamPot, const REAL msg, const INDEX dim)
   {
      loopRight_.loop(dim, [&](const INDEX i) { repamPot[i] = repamPot[i] + msg; });
   }
   // get message for i-th dimension of left potential
   // do zrobienia: prawidlowo?
   // do zrobinia: currently only for 2-dim
   // for ImplicitRepamStorage
   template<typename M>
   const REAL GetLeftMessage(const INDEX i, const M& msg) const
   {
      //assert(i < loopLeft_.GetDim(0)*loopLeft_.GetDim(1));
      //return -msg[ i%loopLeft_.GetDim(0) ];
      return -msg[ loopLeft_.GetMsgIndex(i) ];
   }
   template<typename M>
   const REAL GetRightMessage(const INDEX i, const M& msg) const
   {
      //assert(i < loopRight_.GetDim(0)*loopRight_.GetDim(1));
      //return msg[ i/loopRight_.GetDim(1) ]; // do zrobienia: correct?
      return msg[ loopRight_.GetMsgIndex(i) ];
   }
   
private:
   template<class OP, typename FACTOR_TYPE, typename G1, typename G2, class LOOP>
   void Optimize(OP op, const G1& pot, G2& msg, FACTOR_TYPE* fac, const INDEX k, LOOP& l, const REAL omega); 

   template<typename LEFT_FACTOR, typename G1, typename G2>
   void MaximizeLeft(LEFT_FACTOR* const l, const G1& leftPot, G2& msg, const REAL omega = 1.0)
   {
      // do this with template lambdas
      //auto plus_op = [](decltype(msg[0]) x, const REAL y) -> REAL { return x+y; };
      ;
      Optimize(
            [](const REAL y) -> REAL { return y; },
            leftPot, msg, l, kL_, loopLeft_, omega);
   }
   template<typename RIGHT_FACTOR, typename G1, typename G2>
   void MinimizeRight(RIGHT_FACTOR* r, const G1& rightPot, G2& msg, const REAL omega = 1.0)
   {
      // do this with template lambdas
      //auto minus_op = [](decltype(msg[0]) x, const REAL y) -> REAL { return x-y; };
      Optimize(
            [](const REAL y) -> REAL { return -y; }, 
            rightPot, msg, r, kR_, loopRight_, omega);
   }

   LEFT_LOOP_TYPE loopLeft_;
   RIGHT_LOOP_TYPE loopRight_;
   // do zrobienia: replace this with templates or get from factor
   const INDEX kL_, kR_; // multipliers, see explanation of the model above
};

template<class LEFT_FACTOR_TYPE, class RIGHT_FACTOR_TYPE, bool LEFT_SIDE_ACTIVE, bool RIGHT_SIDE_ACTIVE>
//template<typename MultiplexMargMessage<LEFT_FACTOR_TYPE,RIGHT_FACTOR_TYPE>::OptimizeOp Op, class FACTOR_TYPE, typename G1, typename G2, class LOOP>
template<typename OP, class FACTOR_TYPE, typename G1, typename G2, class LOOP>
void 
MultiplexMargMessage<LEFT_FACTOR_TYPE,RIGHT_FACTOR_TYPE,LEFT_SIDE_ACTIVE,RIGHT_SIDE_ACTIVE>::Optimize(OP op, const G1& pot, G2& msg, FACTOR_TYPE* fac, const INDEX k, LOOP& l, const REAL omega)
//MultiplexMargMessage::Optimize(const G1& pot, G2& msg, FACTOR_TYPE* fac, const INDEX k, LOOP& l, const REAL omega)
{
   //const REAL breakPoINDEX = fac->GetBreakpoINDEXCost(pot);
   assert(k == 1); // do zrobienia: templatize for k such that it may be a constexpr, thus fast, possibly get k directly from factor

   REAL delta;
   l.loop( 
         [&](const INDEX outer_idx){ 
         delta = std::numeric_limits<REAL>::max(); 
         }, 
         [&](const INDEX full_idx, const INDEX outer_idx){ 
         delta = std::min(delta, pot[ full_idx ]);
         //delta = std::min(delta, (pot[ full_idx ] - breakPoINDEX)/REAL(k)); // this is very slow, makes algorithm run two times longer
         },
         [&](const INDEX outer_idx){ 
         assert( outer_idx < msg.size() );
         msg[ outer_idx ] = msg[ outer_idx ] + omega*op(delta);
         } );
}

} // end namespace LP_MP

#endif // LP_MP_MULTIPLEX_MARG_MESSAGE
