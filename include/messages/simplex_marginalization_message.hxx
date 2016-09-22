#ifndef LP_MP_SIMPLEX_MARGINALIZATION_MESSAGE_HXX
#define LP_MP_SIMPLEX_MARGINALIZATION_MESSAGE_HXX

#include "LP_MP.h"
#include "factors_messages.hxx"
#include "message_loops.hxx"

namespace LP_MP {


  // l_11 + ... + l_1n = r_11 + ... + r_l1
  // ...
  // l_k1 + ... + l_kn = r_1k + ... + r_lk
  // all variables must be distinct
  // left and right factor must be simplices and all variables of left and right simplex must be covered

  // {LEFT|RIGHT}_SIDE_ACTIVE activates by SFINAE message updates for left and right factor resp.

  // do zrobienia: rename to SimplexMarginalizationMessage. Simplex will not be supported anymore. Throw away the associated types

  // stride classes are used for cases when variables accessed in reparametrization are not 0,...,n-1, e.g. for min cost flow factor
  /*
    struct ConstantStrideOffset {
    const INDEX operator[](const INDEX i) const { return offset_ + stride_*i; }
    private:
    const INDEX offset_;
    const INDEX stride_;
    };
    struct UniformStrideOffset {
    const INDEX operator[](const INDEX i) const { return offset_ + i; }
    private:
    const INDEX offset_;
    };
  */
  struct UniformStride {
    const INDEX operator[](const INDEX i) const { return i; }
  };

  // replace all thos bools by named parameters via enum classes
  template<class LEFT_LOOP_TYPE, class RIGHT_LOOP_TYPE, bool LEFT_SIDE_ACTIVE, bool RIGHT_SIDE_ACTIVE, bool PROPAGATE_PRIMAL_TO_LEFT = false, bool PROPAGATE_PRIMAL_TO_RIGHT = false, typename LEFT_STRIDE = UniformStride, typename RIGHT_STRIDE = UniformStride>
  class SimplexMarginalizationMessage
  {
  public:
    typedef LEFT_LOOP_TYPE LeftLoopType;
    typedef RIGHT_LOOP_TYPE RightLoopType;

    SimplexMarginalizationMessage(LEFT_LOOP_TYPE loopLeft, RIGHT_LOOP_TYPE loopRight)
      : loopLeft_(loopLeft),
	loopRight_(loopRight)
    {}
    SimplexMarginalizationMessage(LEFT_LOOP_TYPE loopLeft, RIGHT_LOOP_TYPE loopRight, LEFT_STRIDE leftStride, RIGHT_STRIDE rightStride)
      : loopLeft_(loopLeft),
	loopRight_(loopRight),
	leftStride_(leftStride),
	rightStride_(rightStride)
    {}

    // standard functions which take all possible arguments, to be replaced with minimal ones
    template<typename RIGHT_FACTOR, typename G1, typename G2, bool LSA = LEFT_SIDE_ACTIVE>
    typename std::enable_if<LSA,void>::type
    ReceiveMessageFromRight(RIGHT_FACTOR* const r, const G1& rightPot, G2& msg) 
    {
      static_assert(LSA == LEFT_SIDE_ACTIVE,"");
      MinimizeRight(r,rightPot,msg); 
    }
    template<typename LEFT_FACTOR, typename G1, typename G2, bool RSA = RIGHT_SIDE_ACTIVE>
    typename std::enable_if<RSA,void>::type 
    ReceiveMessageFromLeft(LEFT_FACTOR* l, const G1& leftPot, G2& msg) 
    { 
      static_assert(RSA == RIGHT_SIDE_ACTIVE,"");
      MaximizeLeft(l,leftPot,msg); 
    }
    template<typename LEFT_FACTOR, typename G1, typename G3, bool LSA = LEFT_SIDE_ACTIVE>
    typename std::enable_if<LSA,void>::type
    SendMessageToRight(LEFT_FACTOR* const l, const G1& leftPot, G3& msg, const REAL omega)
    { 
      static_assert(LSA == LEFT_SIDE_ACTIVE,"");
      MaximizeLeft(l,leftPot,msg,omega); 
    }
    template<typename RIGHT_FACTOR, typename G2, typename G3, bool RSA = RIGHT_SIDE_ACTIVE>
    typename std::enable_if<RSA,void>::type
    SendMessageToLeft(RIGHT_FACTOR* const r, const G2& rightPot, G3& msg, const REAL omega) 
    { 
      static_assert(RSA == RIGHT_SIDE_ACTIVE,"");
      MinimizeRight(r,rightPot,msg,omega); 
    }

    // for primal computation as in TRW-S, we need to compute restricted messages as well
    template<typename RIGHT_FACTOR, typename G1, typename G2, bool LSA = LEFT_SIDE_ACTIVE>
    typename std::enable_if<LSA,void>::type
    ReceiveRestrictedMessageFromRight(RIGHT_FACTOR* const r, const G1& rightPot, G2& msg, typename PrimalSolutionStorage::Element rightPrimal) 
    {
       static_assert(LSA == LEFT_SIDE_ACTIVE,"");
       OptimizeRestricted(rightPot, msg, r, loopRight_, rightPrimal);
    }
    template<typename LEFT_FACTOR, typename G1, typename G2, bool RSA = RIGHT_SIDE_ACTIVE>
    typename std::enable_if<RSA,void>::type 
    ReceiveRestrictedMessageFromLeft(LEFT_FACTOR* l, const G1& leftPot, G2& msg, typename PrimalSolutionStorage::Element leftPrimal)
    { 
       static_assert(RSA == RIGHT_SIDE_ACTIVE,"");
       OptimizeRestricted(leftPot, msg, l, loopLeft_, leftPrimal);
    }
   

    /*
    // do zrobienia: only templates, fill them
    template<typename MSG_ARRAY, typename RIGHT_REPAM, typename ITERATOR>
    static void
    SendMessagesToLeft(const typename std::enable_if<RIGHT_SIDE_ACTIVE, MSG_ARRAY>::type & msgs, const RIGHT_REPAM& leftRepam, ITERATOR omegaIt)
    {
    if(RIGHT_SIDE_ACTIVE) {
         
    }
    }

    template<typename MSG_ARRAY, typename LEFT_REPAM, typename ITERATOR>
    static void
    SendMessagesToRight(const typename std::enable_if<LEFT_SIDE_ACTIVE, MSG_ARRAY&>::type msgs, const LEFT_REPAM& leftRepam, ITERATOR omegaIt)
    { 
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
      const REAL breakPoint = msgs[0]->GetLeftFactor()->GetBreakpoINDEXCost(pot);
      const REAL k = msgs[0]->msg_op_.kL_;
      leftLoop.loop( 
      [&](const INDEX outer_idx){ delta = std::numeric_limits<REAL>::max(); 
      //std::cout << "outer_idx = " << outer_idx << "\n";
      }, 
      [&](const INDEX full_idx, const INDEX outer_idx){ 
      //std::cout << "full_idx = " << full_idx <<  ", outer_idx = " << outer_idx << "\n";
      assert(full_idx < pot.size());
      if(pot[ full_idx ] > breakPoint) { delta = std::min(delta, (pot[ full_idx ] - breakPoint)/REAL(k)); }
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
    // do zrobienia: put strides in here and below
    template<typename G>
    void RepamLeft(G& repamPot, const REAL msg, const INDEX msg_dim)
    {
      REAL msgn = 0;
      if( std::isfinite(msg) ){ msgn = msg; }
      else{ msgn = std::numeric_limits<REAL>::infinity(); }
      loopLeft_.loop(msg_dim, [&](const INDEX i) { repamPot[leftStride_[i]] += msgn; });
    }
    template<typename G>
    void RepamRight(G& repamPot, const REAL msg, const INDEX dim)
    {
      REAL msgn = 0;
      if( std::isfinite(msg) ){ msgn = msg; }
      else{ msgn = std::numeric_limits<REAL>::infinity(); }
      loopRight_.loop(dim, [&](const INDEX i) { repamPot[rightStride_[i]] += msgn; });
    }
    // get message for i-th dimension of left potential
    // do zrobienia: prawidlowo?
    // do zrobinia: currently only for 2-dim
    // used in conjunction with ImplicitRepamStorage
    template<typename MSG>
    const REAL GetLeftMessage(const INDEX pot_dim, const MSG& msg) const
    {
      assert(false); // sign not correct
      return -msg[ loopLeft_.GetMsgIndex(pot_dim) ];
    }
    template<typename MSG>
    const REAL GetRightMessage(const INDEX i, const MSG& msg) const
    {
      assert(false); // sign not correct
      return msg[ loopRight_.GetMsgIndex(i) ];
    }

    // note: the two functions below only make sense, if loop type is pairwise.
    // in general, labels are propagated as follows: 
    // initally, the primal is set to all unknown vector. The left and right zero out the vector, and exactly one entry should remain unknown, which is then forced to true
    template<bool PROPAGATE_PRIMAL_TO_LEFT_TMP = PROPAGATE_PRIMAL_TO_LEFT, typename LEFT_FACTOR, typename RIGHT_FACTOR>
    typename std::enable_if<PROPAGATE_PRIMAL_TO_LEFT_TMP,void>::type
    ComputeLeftFromRightPrimal(const typename PrimalSolutionStorage::Element left, LEFT_FACTOR* l, typename PrimalSolutionStorage::Element right, RIGHT_FACTOR* r)
    {
       static_assert(PROPAGATE_PRIMAL_TO_LEFT_TMP == PROPAGATE_PRIMAL_TO_LEFT,"");
       loopLeft_.PropagateLabel(right,left);
    }

    template<bool PROPAGATE_PRIMAL_TO_RIGHT_TMP = PROPAGATE_PRIMAL_TO_RIGHT, typename LEFT_FACTOR, typename RIGHT_FACTOR>
    typename std::enable_if<PROPAGATE_PRIMAL_TO_RIGHT_TMP,void>::type
    ComputeRightFromLeftPrimal(const typename PrimalSolutionStorage::Element left, LEFT_FACTOR* l, typename PrimalSolutionStorage::Element right, RIGHT_FACTOR* r)
    {
       static_assert(PROPAGATE_PRIMAL_TO_RIGHT_TMP == PROPAGATE_PRIMAL_TO_RIGHT,"");
       loopRight_.PropagateLabel(left,right);
    }
   
   
    template<class LEFT_FACTOR_TYPE,class RIGHT_FACTOR_TYPE>
    void CreateConstraints(LpInterfaceAdapter* lp,LEFT_FACTOR_TYPE* LeftFactor,RIGHT_FACTOR_TYPE* RightFactor) const
    { 
      INDEX no_constraints = 0;
      for(auto i=0; i<LeftFactor->size(); ++i){ 
        no_constraints = std::max(loopLeft_.GetMsgIndex(i), no_constraints);
      } 
      no_constraints++; 
    
      std::vector<std::vector<INDEX>> leftIndices(no_constraints), rightIndices(no_constraints);
      
      for(auto i=0; i<LeftFactor->size(); ++i){
        leftIndices[loopLeft_.GetMsgIndex(i)].push_back(i); 
      }
      for(auto i=0; i<RightFactor->size(); ++i){
        rightIndices[loopRight_.GetMsgIndex(i)].push_back(i); 
      }
          
      for(auto i=0; i<no_constraints ; i++){
        LinExpr lhs = lp->CreateLinExpr();
        LinExpr rhs = lp->CreateLinExpr();
        for(auto l=0; l<leftIndices[i].size(); l++){
          lhs += lp->GetLeftVariable(leftIndices[i][l]);
        }
        for(auto r=0; r<rightIndices[i].size(); r++){
          rhs += lp->GetRightVariable(rightIndices[i][r]);
        }
        lp->addLinearEquality(lhs,rhs);
      }
    }
  
  private:
    template<typename LEFT_FACTOR, typename G1, typename G2>
    void MaximizeLeft(LEFT_FACTOR* const l, const G1& leftPot, G2& msg, const REAL omega = 1.0)
    {
      Optimize(
	       leftPot, msg, l, loopLeft_, omega);
    }
    template<typename RIGHT_FACTOR, typename G1, typename G2>
    void MinimizeRight(RIGHT_FACTOR* r, const G1& rightPot, G2& msg, const REAL omega = 1.0)
    {
      Optimize(
	       rightPot, msg, r, loopRight_, omega);
    }

    // do zrobienia: derive from the four classes below, to enable empty base class optimization.
    LEFT_LOOP_TYPE loopLeft_;
    RIGHT_LOOP_TYPE loopRight_;

    LEFT_STRIDE leftStride_;
    RIGHT_STRIDE rightStride_;

    template<class FACTOR_TYPE, typename G1, typename G2, class LOOP>
    void Optimize(const G1& pot, G2& msg, FACTOR_TYPE* fac, LOOP& l, const REAL omega)
    {
      REAL delta;
      // do zrobienia: normalize by minimum value
      l.loop( 
	     [&](const INDEX outer_idx){ 
	       delta = std::numeric_limits<REAL>::infinity(); 
	     }, 
	     [&](const INDEX full_idx, const INDEX outer_idx){ 
	       delta = std::min(delta, pot[ full_idx ]);
	     },
	     [&](const INDEX outer_idx){ 
	       assert( outer_idx < msg.size() );
	       msg[ outer_idx ] -= omega*delta;
	     } );
    }


    template<class FACTOR_TYPE, typename G1, typename G2, class LOOP>
    void OptimizeRestricted(const G1& pot, G2& msg, FACTOR_TYPE* fac, LOOP& l, const typename PrimalSolutionStorage::Element primal)
    {
      // only take into account those entries which are not false
      REAL delta;
      l.loop(
	     [&](const INDEX outer_idx) {
	       delta = std::numeric_limits<REAL>::infinity();
	     },
	     [&](const INDEX full_idx, const INDEX outer_idx){ 
	       if(primal[ full_idx ] == unknownState || primal[ full_idx ] == true) { // i.e. may be true or unknown. do zrobienia: if true, then go back, set delta to infinity for all other coordinates and this coordinate let be -infinity
		 delta = std::min(delta, pot[ full_idx ]);
	       } 
	       //if(primal[ full_idx ] == true) {
	       //assert(false); // need to implement the above described behaviour
	       //delta = -std::numeric_limits<REAL>::infinity();
	       //}
	     },
	     [&](const INDEX outer_idx){ 
	       assert(delta > -std::numeric_limits<REAL>::infinity());
	       assert( outer_idx < msg.size() );
	       msg[ outer_idx ] -= delta;
	     } );
    }
  };

} // end namespace LP_MP

#endif // LP_MP_SIMPLEX_MARGINALIZATION_MESSAGE_HXX
