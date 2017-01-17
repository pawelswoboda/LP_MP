#ifndef LP_MP_SIMPLEX_MARGINALIZATION_MESSAGE_HXX
#define LP_MP_SIMPLEX_MARGINALIZATION_MESSAGE_HXX

#include "LP_MP.h"
#include "factors_messages.hxx"
#include "vector.hxx"
#include "message_loops.hxx"
#include "memory_allocator.hxx"

namespace LP_MP {

// specialized messages between UnarySimplexFactor and PairwiseSimplexFactor
  template<MessageSendingType TYPE,  bool PROPAGATE_PRIMAL_TO_LEFT = false, bool PROPAGATE_PRIMAL_TO_RIGHT = false, bool SUPPORT_INFINITY = true> 
  class UnaryPairwiseMessageLeft {
    public:
      UnaryPairwiseMessageLeft(const INDEX i1, const INDEX i2) : i1_(i1), i2_(i2) {} // the pairwise factor size
      // standard functions which take all possible arguments, to be replaced with minimal ones
      template<typename RIGHT_FACTOR, typename G2, bool ENABLE = TYPE == MessageSendingType::SRMP>
        typename std::enable_if<ENABLE,void>::type
        ReceiveMessageFromRight(const RIGHT_FACTOR& r, G2& msg) 
        {
          MinimizeRight(r,msg); 
        }
      template<typename LEFT_FACTOR, typename G2, bool ENABLE = TYPE == MessageSendingType::MPLP>
        typename std::enable_if<ENABLE,void>::type 
        ReceiveMessageFromLeft(const LEFT_FACTOR& l, G2& msg) 
        { 
          MaximizeLeft(l,msg); 
        }
      template<typename LEFT_FACTOR, typename G3, bool ENABLE = TYPE == MessageSendingType::SRMP>
        typename std::enable_if<ENABLE,void>::type
        SendMessageToRight(const LEFT_FACTOR& l, G3& msg, const REAL omega)
        { 
          MaximizeLeft(l,msg,omega); 
        }
      template<typename RIGHT_FACTOR, typename G3, bool ENABLE = TYPE == MessageSendingType::MPLP>
        typename std::enable_if<ENABLE,void>::type
        SendMessageToLeft(const RIGHT_FACTOR& r, G3& msg, const REAL omega) 
        { 
          MinimizeRight(r,msg,omega); 
        }

      // for primal computation as in TRW-S, we need to compute restricted messages as well
      template<typename RIGHT_FACTOR, typename G2, bool ENABLE = TYPE == MessageSendingType::SRMP>
        typename std::enable_if<ENABLE,void>::type
        ReceiveRestrictedMessageFromRight(const RIGHT_FACTOR& r, G2& msg) 
        {
           // we assume that only r.right_primal was assigned, r.left_primal not
           assert(r.left_primal_ == i1_);
           if(r.right_primal_ < i2_) {
              vector msgs(i1_);
              for(INDEX x1=0; x1<i1_; ++x1) {
                 msgs[x1] = r(x1,r.right_primal_);
              }
              msg -= msgs;
           }
        }

    // reparametrize left potential for i-th entry of msg
    // do zrobienia: put strides in here and below
    template<typename G>
    void RepamLeft(G& r, const REAL msg, const INDEX msg_dim)
    {
       assert(!std::isnan(msg));
       if(SUPPORT_INFINITY) {
          r[msg_dim] += normalize( msg );
       } else {
          r[msg_dim] += msg;
       }
    }
    template<typename A1, typename A2>
    void RepamRight(A1& r, const A2& msgs)
    {
       for(INDEX x1=0; x1<r.dim1(); ++x1) {
          assert(!std::isnan(msgs[x1]));
           if(SUPPORT_INFINITY) {
              r.left(x1) += normalize( msgs[x1] );
           } else {
              r.left(x1) += msgs[x1];
           }
       }
    }
    template<typename G>
    void RepamRight(G& r, const REAL msg, const INDEX dim)
    {
       assert(!std::isnan(msg));
       if(SUPPORT_INFINITY) {
          r.left(dim) += normalize( msg );
       } else {
          r.left(dim) += msg;
       }
    }

    template<bool ENABLE = TYPE == MessageSendingType::SRMP, typename LEFT_FACTOR, typename RIGHT_FACTOR>
    //typename std::enable_if<ENABLE,void>::type
    void
    ComputeRightFromLeftPrimal(const LEFT_FACTOR& l, RIGHT_FACTOR& r)
    {
       assert(l.primal() < i1_);
       assert(r.left_primal_ == i1_);
       r.left_primal_ = l.primal();
    }


    /*
    template<class LEFT_FACTOR_TYPE,class RIGHT_FACTOR_TYPE>
      void CreateConstraints(LpInterfaceAdapter* lp,LEFT_FACTOR_TYPE* LeftFactor,RIGHT_FACTOR_TYPE* RightFactor) const
      { 
        for(auto x1=0; x1<i1_ ; x1++){
          LinExpr lhs = lp->CreateLinExpr() + 0.0;
          LinExpr rhs = lp->CreateLinExpr() + 0.0;
          if(lp->IsLeftObjective(x1)){
            lhs += lp->GetLeftVariable(x1);
          }
          for(auto x2=0; x2<i2_; x2++){
            if(lp->IsRightObjective(x2*i1_ + x1)){
              rhs += lp->GetRightVariable(x2*i1_ + x1);
            }
          }
          lp->addLinearEquality(lhs,rhs);
        }
      }
      */

  private:
    template<typename LEFT_FACTOR, typename G2>
    void MaximizeLeft(const LEFT_FACTOR& l, G2& msg, const REAL omega = 1.0)
    {
      for(INDEX x1=0; x1<i1_; ++x1) {
        msg[x1] -= omega*l[x1];
      }
    }
    template<typename RIGHT_FACTOR, typename G2>
    void MinimizeRight(const RIGHT_FACTOR& r, G2& msg, const REAL omega = 1.0)
    {
       vector msgs(i1_,std::numeric_limits<REAL>::infinity());
       //std::vector<REAL> msgs(i1_,std::numeric_limits<REAL>::infinity());
       //REAL msgs[i1_];
       //for(INDEX x1=0; x1<i1_; ++x1) {
       //  msgs[x1] = std::numeric_limits<REAL>::infinity();
       //}
       for(INDEX x1=0; x1<i1_; ++x1) {
          for(INDEX x2=0; x2<i2_; ++x2) {
             msgs[x1] = std::min(msgs[x1],r(x1,x2));
             //msgs[x1] = std::min(msgs[x1],omega*r[x2*i1_ + x1]);
          }
       }
       msg -= omega*msgs;
    }

    const INDEX i1_,i2_;
  };


  template<MessageSendingType TYPE,  bool PROPAGATE_PRIMAL_TO_LEFT = false, bool PROPAGATE_PRIMAL_TO_RIGHT = false, bool SUPPORT_INFINITY = true> 
  class UnaryPairwiseMessageRight {
    public:
      UnaryPairwiseMessageRight(const INDEX i1, const INDEX i2) : i1_(i1), i2_(i2) {} // the pairwise factor size
      // standard functions which take all possible arguments, to be replaced with minimal ones
      template<typename RIGHT_FACTOR, typename G2, bool ENABLE = TYPE == MessageSendingType::SRMP>
        typename std::enable_if<ENABLE,void>::type
        ReceiveMessageFromRight(const RIGHT_FACTOR& r, G2& msg) 
        {
          MinimizeRight(r,msg); 
        }
      template<typename LEFT_FACTOR, typename G2, bool ENABLE = TYPE == MessageSendingType::MPLP>
        typename std::enable_if<ENABLE,void>::type 
        ReceiveMessageFromLeft(const LEFT_FACTOR& l, G2& msg) 
        { 
          MaximizeLeft(l,msg); 
        }
      template<typename LEFT_FACTOR, typename G3, bool ENABLE = TYPE == MessageSendingType::SRMP>
        typename std::enable_if<ENABLE,void>::type
        SendMessageToRight(const LEFT_FACTOR& l, G3& msg, const REAL omega)
        { 
          MaximizeLeft(l,msg,omega); 
        }
      template<typename RIGHT_FACTOR, typename G3, bool ENABLE = TYPE == MessageSendingType::MPLP>
        typename std::enable_if<ENABLE,void>::type
        SendMessageToLeft(const RIGHT_FACTOR& r, G3& msg, const REAL omega) 
        { 
          MinimizeRight(r,msg,omega); 
        }

      // for primal computation as in TRW-S, we need to compute restricted messages as well
      template<typename RIGHT_FACTOR, typename G2, bool ENABLE = TYPE == MessageSendingType::SRMP>
        typename std::enable_if<ENABLE,void>::type
        ReceiveRestrictedMessageFromRight(const RIGHT_FACTOR& r, G2& msg) 
        {
           assert(r.right_primal_ == i1_);
           if(r.left_primal_ < i2_) {
              vector msgs(i2_);
              for(INDEX x2=0; x2<i2_; ++x2) {
                 msgs[x2] = r(r.left_primal_,x2);
              }
              msg -= msgs;
           }
        }

    // reparametrize left potential for i-th entry of msg
    // do zrobienia: put strides in here and below
    template<typename G>
    void RepamLeft(G& l, const REAL msg, const INDEX msg_dim)
    {
       assert(!std::isnan(msg));
       if(SUPPORT_INFINITY) {
          l[msg_dim] += normalize( msg );
       } else {
          l[msg_dim] += msg;
       }
    }
    template<typename A1, typename A2>
    void RepamRight(A1& r, const A2& msgs)
    {
      for(INDEX x2=0; x2<i2_; ++x2) {
       assert(!std::isnan(msgs[x2]));
         if(SUPPORT_INFINITY) {
            r.right(x2) += normalize( msgs[x2] );
         } else {
            r.right(x2) += msgs[x2];
         }
      }
    }
    template<typename G>
    void RepamRight(G& r, const REAL msg, const INDEX dim)
    {
       assert(!std::isnan(msg));
       if(SUPPORT_INFINITY) {
          r.right(dim) += normalize( msg );
       } else {
          r.right(dim) += msg;
       }
    }

    template<bool ENABLE = TYPE == MessageSendingType::SRMP, typename LEFT_FACTOR, typename RIGHT_FACTOR>
    //typename std::enable_if<ENABLE,void>::type
    void
    ComputeRightFromLeftPrimal(const LEFT_FACTOR& l, RIGHT_FACTOR& r)
    {
       assert(l.primal() < i2_);
       assert(r.right_primal_ == i2_);
       r.right_primal_ = l.primal();
    }


    /*
    template<class LEFT_FACTOR_TYPE,class RIGHT_FACTOR_TYPE>
      void CreateConstraints(LpInterfaceAdapter* lp,LEFT_FACTOR_TYPE* LeftFactor,RIGHT_FACTOR_TYPE* RightFactor) const
      { 
        for(auto x2=0; x2<i2_ ; x2++){
          LinExpr lhs = lp->CreateLinExpr() + 0.0;
          LinExpr rhs = lp->CreateLinExpr() + 0.0;
          if(lp->IsLeftObjective(x2)){
            lhs += lp->GetLeftVariable(x2);
          }
          for(auto x1=0; x1<i1_; x1++){
            if(lp->IsRightObjective(x2*i1_ + x1)){
              rhs += lp->GetRightVariable(x2*i1_ + x1);
            }
          }
          lp->addLinearEquality(lhs,rhs);
        }
      }
      */

  private:
    template<typename LEFT_FACTOR, typename G2>
    void MaximizeLeft(const LEFT_FACTOR& l, G2& msg, const REAL omega = 1.0)
    {
      for(INDEX x2=0; x2<i2_; ++x2) {
        msg[x2] -= omega*l[x2];
      }
    }
    template<typename RIGHT_FACTOR, typename G2>
    void MinimizeRight(const RIGHT_FACTOR& r, G2& msg, const REAL omega = 1.0)
    {
      vector msgs(i2_, std::numeric_limits<REAL>::infinity());//std::vector<REAL, stack_allocator<REAL>>;
      //vector_type msgs(i2_,std::numeric_limits<REAL>::infinity(), global_real_stack_allocator);
      //std::vector<REAL> msgs(i2_,std::numeric_limits<REAL>::infinity());
      for(INDEX x1=0; x1<i1_; ++x1) {
         for(INDEX x2=0; x2<i2_; ++x2) {
          msgs[x2] = std::min(msgs[x2],omega*r(x1,x2));
          //msgs[x2] = std::min(msgs[x2],omega*r[x2*i1_ + x1]);
        }
      }
      msg -= msgs;
    }

    const INDEX i1_,i2_; // do zrobienia: these values are not needed, as they can be obtained from the factors whenever they are used
  };

  // specialized messages for pairwise/triplet marginalization
  template<MessageSendingType TYPE,  bool PROPAGATE_PRIMAL_TO_LEFT = false, bool PROPAGATE_PRIMAL_TO_RIGHT = false> 
  class PairwiseTripletMessage12 {
    public:
      PairwiseTripletMessage12(const INDEX i1, const INDEX i2, const INDEX i3) : i1_(i1), i2_(i2), i3_(i3) {} // the pairwise factor size
      // standard functions which take all possible arguments, to be replaced with minimal ones
      template<typename RIGHT_FACTOR, typename G2, bool ENABLE = TYPE == MessageSendingType::SRMP>
        typename std::enable_if<ENABLE,void>::type
        ReceiveMessageFromRight(const RIGHT_FACTOR& r, G2& msg) 
        {
          MinimizeRight(r,msg); 
        }
      template<typename LEFT_FACTOR, typename G2, bool ENABLE = TYPE == MessageSendingType::MPLP>
        typename std::enable_if<ENABLE,void>::type 
        ReceiveMessageFromLeft(const LEFT_FACTOR& l, G2& msg) 
        { 
          MaximizeLeft(l,msg); 
        }
      template<typename LEFT_FACTOR, typename G3, bool ENABLE = TYPE == MessageSendingType::SRMP>
        typename std::enable_if<ENABLE,void>::type
        SendMessageToRight(const LEFT_FACTOR& l, G3& msg, const REAL omega)
        { 
          MaximizeLeft(l,msg,omega); 
        }
      template<typename RIGHT_FACTOR, typename G3, bool ENABLE = TYPE == MessageSendingType::MPLP>
        typename std::enable_if<ENABLE,void>::type
        SendMessageToLeft(const RIGHT_FACTOR& r, G3& msg, const REAL omega) 
        { 
          MinimizeRight(r,msg,omega); 
        }

      // for primal computation as in TRW-S, we need to compute restricted messages as well
      template<typename RIGHT_FACTOR, typename G2, bool ENABLE = TYPE == MessageSendingType::SRMP>
        typename std::enable_if<ENABLE,void>::type
        ReceiveRestrictedMessageFromRight(const RIGHT_FACTOR& r, G2& msg, typename PrimalSolutionStorage::Element rightPrimal) 
        {
          assert(false);
        }

    // reparametrize left potential for i-th entry of msg
    // do zrobienia: put strides in here and below
      /*
    template<typename G>
    void RepamLeft(G& l, const REAL msg, const INDEX msg_dim)
    {
      const INDEX x1 = msg_dim/i2_;
      const INDEX x2 = msg_dim%i2_;
      //if(SUPPORT_INFINITY) {
         l(x1,x2) += normalize( msg );
      //} else {
      //   l(x1,x2) += msg;
      //}
    }
    */
    template<typename A1, typename A2>
    void RepamLeft(A1& l, const A2& msgs)
    {
      // do zrobienia: possibly use counter
       for(INDEX x1=0; x1<i1_; ++x1) {
          for(INDEX x2=0; x2<i2_; ++x2) {
             l(x1,x2) += normalize( msgs(x1,x2) );
          }
       }
    }
    template<typename A1, typename A2>
    void RepamRight(A1& r, const A2& msgs)
    {
      // do zrobienia: possibly use counter
       for(INDEX x1=0; x1<i1_; ++x1) {
          for(INDEX x2=0; x2<i2_; ++x2) {
             r.msg12(x1,x2) += normalize( msgs(x1,x2) );
          }
       }
    }
    /*
    template<typename G>
    void RepamRight(G& r, const REAL msg, const INDEX dim)
    {
      const INDEX x1 = dim/i2_;
      const INDEX x2 = dim%i2_;
      r.msg12(x1,x2) += normalize( msg );
    }
    */

    template<bool ENABLE = TYPE == MessageSendingType::SRMP, typename LEFT_FACTOR, typename RIGHT_FACTOR>
    typename std::enable_if<ENABLE,void>::type
    ComputeRightFromLeftPrimal(const typename PrimalSolutionStorage::Element left, LEFT_FACTOR* l, typename PrimalSolutionStorage::Element right, RIGHT_FACTOR* r)
    {
      assert(r->dim1() == i1_ && r->dim2() == i2_ && r->dim3() == i3_);
      assert(l->size() == i1_*i2_);
      for(INDEX x2=0; x2<i2_; ++x2) {
        for(INDEX x1=0; x1<i1_; ++x1) {
          if(left[x2*i1_ + x1] == false) {
            for(INDEX x3=0; x3<i3_; ++x3) {
              right[x3*i2_*i1_ + x2*i1_ + x1] = false;
            }
          } else if(left[x2*i1_ + x1] == true) {
            INDEX no_false = 0;
            for(INDEX x3=0; x3<i3_; ++x3) {
              if(right[x3*i1_*i2_ + x2*i1_ + x1] == false) {
                ++no_false;
              }
            }
            if(no_false == i3_-1) {
              for(INDEX x3=0; x3<i3_; ++x3) {
                if(right[x3*i1_*i2_ + x2*i1_ + x1] == unknownState) {
                  right[x3*i1_*i2_ + x2*i1_ + x1] = true;
                }
              }
            }
          } else {
            // may happen when only part of pairwise was labelled until now
          }
        }
      }
    }


    /*
    template<class LEFT_FACTOR_TYPE,class RIGHT_FACTOR_TYPE>
      void CreateConstraints(LpInterfaceAdapter* lp,LEFT_FACTOR_TYPE* LeftFactor,RIGHT_FACTOR_TYPE* RightFactor) const
      { 
        for(auto x1=0; x1<i1_ ; x1++){
          for(auto x2=0; x2<i2_; x2++){
            LinExpr lhs = lp->CreateLinExpr() + 0.0;
            LinExpr rhs = lp->CreateLinExpr() + 0.0;
            if(lp->IsLeftObjective(x2*i1_ + x1)){
              lhs += lp->GetLeftVariable(x2*i1_ + x1);
            }
          for(auto x3=0; x3<i3_; x3++){
            if(lp->IsRightObjective(x3*i1_*i2_ + x2*i1_ + x1)){
              rhs += lp->GetRightVariable(x3*i1_+i2_ + x2*i1_ + x1);
            }
          }
          lp->addLinearEquality(lhs,rhs);
          }
        }
      }
      */

  private:
    template<typename LEFT_FACTOR, typename G2>
    void MaximizeLeft(const LEFT_FACTOR& l, G2& msg, const REAL omega = 1.0)
    {
       msg -= omega*l;
      //for(INDEX x2=0; x2<i2_; ++x2) {
      //  for(INDEX x1=0; x1<i1_; ++x1) {
      //    msg[x2*i1_ + x1] -= omega*l[x2*i1_ + x1];
      //  }
      //}
    }
    template<typename RIGHT_FACTOR, typename G2>
    void MinimizeRight(const RIGHT_FACTOR& r, G2& msg, const REAL omega = 1.0)
    {
       matrix msgs(i1_,i2_, std::numeric_limits<REAL>::infinity());
       r.min_marginal12(msgs);
       msg -= omega*msgs;
    }

    const INDEX i1_,i2_, i3_;
  };

  template<MessageSendingType TYPE,  bool PROPAGATE_PRIMAL_TO_LEFT = false, bool PROPAGATE_PRIMAL_TO_RIGHT = false> 
  class PairwiseTripletMessage13 {
    public:
      PairwiseTripletMessage13(const INDEX i1, const INDEX i2, const INDEX i3) : i1_(i1), i2_(i2), i3_(i3) {} // the pairwise factor size
      // standard functions which take all possible arguments, to be replaced with minimal ones
      template<typename RIGHT_FACTOR, typename G2, bool ENABLE = TYPE == MessageSendingType::SRMP>
        typename std::enable_if<ENABLE,void>::type
        ReceiveMessageFromRight(const RIGHT_FACTOR& r, G2& msg) 
        {
          MinimizeRight(r,msg); 
        }
      template<typename LEFT_FACTOR, typename G2, bool ENABLE = TYPE == MessageSendingType::MPLP>
        typename std::enable_if<ENABLE,void>::type 
        ReceiveMessageFromLeft(const LEFT_FACTOR& l, G2& msg) 
        { 
          MaximizeLeft(l,msg); 
        }
      template<typename LEFT_FACTOR, typename G3, bool ENABLE = TYPE == MessageSendingType::SRMP>
        typename std::enable_if<ENABLE,void>::type
        SendMessageToRight(const LEFT_FACTOR& l, G3& msg, const REAL omega)
        { 
          MaximizeLeft(l,msg,omega); 
        }
      template<typename RIGHT_FACTOR, typename G3, bool ENABLE = TYPE == MessageSendingType::MPLP>
        typename std::enable_if<ENABLE,void>::type
        SendMessageToLeft(const RIGHT_FACTOR& r, G3& msg, const REAL omega) 
        { 
          MinimizeRight(r,msg,omega); 
        }

      // for primal computation as in TRW-S, we need to compute restricted messages as well
      template<typename RIGHT_FACTOR, typename G2, bool ENABLE = TYPE == MessageSendingType::SRMP>
        typename std::enable_if<ENABLE,void>::type
        ReceiveRestrictedMessageFromRight(const RIGHT_FACTOR& r, G2& msg, typename PrimalSolutionStorage::Element rightPrimal) 
        {
          assert(false);
        }

    // reparametrize left potential for i-th entry of msg
    // do zrobienia: put strides in here and below
    /*
    template<typename G>
    void RepamLeft(G& repamPot, const REAL msg, const INDEX msg_dim)
    {
      REAL msgn;
      if( std::isfinite(msg) ){ msgn = msg; }
      else{ msgn = std::numeric_limits<REAL>::infinity(); }
      const INDEX x1 = msg_dim/i3_;
      const INDEX x3 = msg_dim%i3_;
      repamPot(x1,x3) += msgn;
    }
    */
    template<typename A1, typename A2>
    void RepamLeft(A1& repamPot, const A2& msgs)
    {
       // do zrobienia: possibly use counter
       for(INDEX x1=0; x1<i1_; ++x1) {
          for(INDEX x3=0; x3<i3_; ++x3) {
             repamPot(x1,x3) += normalize( msgs(x1,x3) );
          }
       }
    }
    template<typename A1, typename A2>
    void RepamRight(A1& repamPot, const A2& msgs)
    {
       // do zrobienia: possibly use counter
       for(INDEX x1=0; x1<i1_; ++x1) {
          for(INDEX x3=0; x3<i3_; ++x3) {
             repamPot.msg13(x1,x3) += normalize( msgs(x1,x3) );
          }
       }
    }
    /*
    template<typename G>
    void RepamRight(G& repamPot, const REAL msg, const INDEX dim)
    {
      REAL msgn;
      const INDEX x1 = dim % i1_;
      const INDEX x3 = dim / i1_;
      if( std::isfinite(msg) ){ msgn = msg; }
      else{ msgn = std::numeric_limits<REAL>::infinity(); }
      assert(false);
      //for(INDEX x2 = 0; x2<i2_; ++x2) {
      //  repamPot[x3*i2_*i1_ + x2*i1_ + x1] += msgn;
      //}
    }
    */

    template<bool ENABLE = TYPE == MessageSendingType::SRMP, typename LEFT_FACTOR, typename RIGHT_FACTOR>
    typename std::enable_if<ENABLE,void>::type
    ComputeRightFromLeftPrimal(const typename PrimalSolutionStorage::Element left, LEFT_FACTOR* l, typename PrimalSolutionStorage::Element right, RIGHT_FACTOR* r)
    {
      assert(r->dim1() == i1_ && r->dim2() == i2_ && r->dim3() == i3_);
      assert(l->size() == i1_*i3_);
      for(INDEX x3=0; x3<i3_; ++x3) {
        for(INDEX x1=0; x1<i1_; ++x1) {
          if(left[x3*i1_ + x1] == false) {
            for(INDEX x2=0; x2<i2_; ++x2) {
              right[x3*i2_*i1_ + x2*i1_ + x1] = false;
            }
          } else if(left[x3*i1_ + x1] == true) {
            INDEX no_false = 0;
            for(INDEX x2=0; x2<i2_; ++x2) {
              if(right[x3*i1_*i2_ + x2*i1_ + x1] == false) {
                ++no_false;
              }
            }
            if(no_false == i2_-1) {
              for(INDEX x2=0; x2<i2_; ++x2) {
                if(right[x3*i1_*i2_ + x2*i1_ + x1] == unknownState) {
                  right[x3*i1_*i2_ + x2*i1_ + x1] = true;
                }
              }
            }
          } else {
            // may happen when only part of pairwise was labelled until now
          }
        }
      }
    }


    /*
    template<class LEFT_FACTOR_TYPE,class RIGHT_FACTOR_TYPE>
      void CreateConstraints(LpInterfaceAdapter* lp,LEFT_FACTOR_TYPE* LeftFactor,RIGHT_FACTOR_TYPE* RightFactor) const
      { 
        for(auto x1=0; x1<i1_ ; x1++){
          for(auto x3=0; x3<i3_; x3++){
            LinExpr lhs = lp->CreateLinExpr() + 0.0;
            LinExpr rhs = lp->CreateLinExpr() + 0.0;
            if(lp->IsLeftObjective(x3*i1_ + x1)){
              lhs += lp->GetLeftVariable(x3*i1_ + x1);
            }
          for(auto x2=0; x2<i2_; x2++){
            if(lp->IsRightObjective(x3*i1_*i2_ + x2*i1_ + x1)){
              rhs += lp->GetRightVariable(x3*i1_+i2_ + x2*i1_ + x1);
            }
          }
          lp->addLinearEquality(lhs,rhs);
          }
        }
      }
      */

  private:
    template<typename LEFT_FACTOR, typename G2>
    void MaximizeLeft(const LEFT_FACTOR& l, G2& msg, const REAL omega = 1.0)
    {
       msg -= omega*l;
    }
    template<typename RIGHT_FACTOR, typename G2>
    void MinimizeRight(const RIGHT_FACTOR& r, G2& msg, const REAL omega = 1.0)
    {
       matrix msgs(i1_,i3_, std::numeric_limits<REAL>::infinity());
       r.min_marginal13(msgs);
       msg -= omega*msgs;
    }

    const INDEX i1_,i2_, i3_;
  };

  template<MessageSendingType TYPE,  bool PROPAGATE_PRIMAL_TO_LEFT = false, bool PROPAGATE_PRIMAL_TO_RIGHT = false> 
  class PairwiseTripletMessage23 {
    public:
      PairwiseTripletMessage23(const INDEX i1, const INDEX i2, const INDEX i3) : i1_(i1), i2_(i2), i3_(i3) {} // the pairwise factor size
      // standard functions which take all possible arguments, to be replaced with minimal ones
      template<typename RIGHT_FACTOR, typename G2, bool ENABLE = TYPE == MessageSendingType::SRMP>
        typename std::enable_if<ENABLE,void>::type
        ReceiveMessageFromRight(const RIGHT_FACTOR& r, G2& msg) 
        {
          MinimizeRight(r,msg); 
        }
      template<typename LEFT_FACTOR, typename G2, bool ENABLE = TYPE == MessageSendingType::MPLP>
        typename std::enable_if<ENABLE,void>::type 
        ReceiveMessageFromLeft(const LEFT_FACTOR& l, G2& msg) 
        { 
          MaximizeLeft(l,msg); 
        }
      template<typename LEFT_FACTOR, typename G3, bool ENABLE = TYPE == MessageSendingType::SRMP>
        typename std::enable_if<ENABLE,void>::type
        SendMessageToRight(const LEFT_FACTOR& l, G3& msg, const REAL omega)
        { 
          MaximizeLeft(l,msg,omega); 
        }
      template<typename RIGHT_FACTOR, typename G3, bool ENABLE = TYPE == MessageSendingType::MPLP>
        typename std::enable_if<ENABLE,void>::type
        SendMessageToLeft(const RIGHT_FACTOR& r, G3& msg, const REAL omega) 
        { 
          MinimizeRight(r,msg,omega); 
        }

      // for primal computation as in TRW-S, we need to compute restricted messages as well
      template<typename RIGHT_FACTOR, typename G1, typename G2, bool ENABLE = TYPE == MessageSendingType::SRMP>
        typename std::enable_if<ENABLE,void>::type
        ReceiveRestrictedMessageFromRight(RIGHT_FACTOR* const r, const G1& rightPot, G2& msg, typename PrimalSolutionStorage::Element rightPrimal) 
        {
          assert(false);
        }

    // reparametrize left potential for i-th entry of msg
    // do zrobienia: put strides in here and below
      /*
    template<typename G>
    void RepamLeft(G& repamPot, const REAL msg, const INDEX msg_dim)
    {
      REAL msgn;
      if( std::isfinite(msg) ){ msgn = msg; }
      else{ msgn = std::numeric_limits<REAL>::infinity(); }
      const INDEX x2 = msg_dim/i3_;
      const INDEX x3 = msg_dim%i3_;
      repamPot(x2,x3) += msgn;
    }
    */
    template<typename A1, typename A2>
    void RepamLeft(A1& repamPot, const A2& msgs)
    {
      // do zrobienia: possibly use counter
       for(INDEX x2=0; x2<i2_; ++x2) {
          for(INDEX x3 = 0; x3<i3_; ++x3) {
             repamPot(x2,x3) += normalize( msgs(x2,x3) );
          }
       }
    }
    template<typename A1, typename A2>
    void RepamRight(A1& repamPot, const A2& msgs)
    {
      // do zrobienia: possibly use counter
       for(INDEX x2=0; x2<i2_; ++x2) {
          for(INDEX x3=0; x3<i3_; ++x3) {
             repamPot.msg23(x2,x3) += normalize( msgs(x2,x3) );
          }
       }
    }
    /*
    template<typename G>
    void RepamRight(G& repamPot, const REAL msg, const INDEX dim)
    {
      REAL msgn;
      const INDEX x2 = dim % i2_;
      const INDEX x3 = dim / i2_;
      if( std::isfinite(msg) ){ msgn = msg; }
      else{ msgn = std::numeric_limits<REAL>::infinity(); }
      assert(false);
      //for(INDEX x1 = 0; x1<i1_; ++x1) {
      //  repamPot[x3*i2_*i1_ + x2*i1_ + x1] += msgn;
      //}
    }
    */

    template<bool ENABLE = TYPE == MessageSendingType::SRMP, typename LEFT_FACTOR, typename RIGHT_FACTOR>
    typename std::enable_if<ENABLE,void>::type
    ComputeRightFromLeftPrimal(const typename PrimalSolutionStorage::Element left, LEFT_FACTOR* l, typename PrimalSolutionStorage::Element right, RIGHT_FACTOR* r)
    {
      assert(r->dim1() == i1_ && r->dim2() == i2_ && r->dim3() == i3_);
      assert(l->size() == i2_*i3_);
      for(INDEX x3=0; x3<i3_; ++x3) {
        for(INDEX x2=0; x2<i2_; ++x2) {
          if(left[x3*i2_ + x2] == false) {
            for(INDEX x1=0; x1<i1_; ++x1) {
              right[x3*i2_*i1_ + x2*i1_ + x1] = false;
            }
          } else if(left[x3*i2_ + x2] == true) {
            INDEX no_false = 0;
            for(INDEX x1=0; x1<i1_; ++x1) {
              if(right[x3*i1_*i2_ + x2*i1_ + x1] == false) {
                ++no_false;
              }
            }
            if(no_false == i1_-1) {
              for(INDEX x1=0; x1<i1_; ++x1) {
                if(right[x3*i1_*i2_ + x2*i1_ + x1] == unknownState) {
                  right[x3*i1_*i2_ + x2*i1_ + x1] = true;
                }
              }
            }
          } else {
            // may happen when only part of pairwise was labelled until now
          }
        }
      }
    }


    /*
    template<class LEFT_FACTOR_TYPE,class RIGHT_FACTOR_TYPE>
      void CreateConstraints(LpInterfaceAdapter* lp,LEFT_FACTOR_TYPE* LeftFactor,RIGHT_FACTOR_TYPE* RightFactor) const
      { 
        for(auto x2=0; x2<i2_ ; x2++){
          for(auto x3=0; x3<i3_; x3++){
            LinExpr lhs = lp->CreateLinExpr() + 0.0;
            LinExpr rhs = lp->CreateLinExpr() + 0.0;
            if(lp->IsLeftObjective(x3*i2_ + x2)){
              lhs += lp->GetLeftVariable(x3*i2_ + x2);
            }
          for(auto x1=0; x1<i1_; x1++){
            if(lp->IsRightObjective(x3*i1_*i2_ + x2*i1_ + x1)){
              rhs += lp->GetRightVariable(x3*i1_+i2_ + x2*i1_ + x1);
            }
          }
          lp->addLinearEquality(lhs,rhs);
          }
        }
      }
      */

  private:
    template<typename LEFT_FACTOR, typename G2>
    void MaximizeLeft(const LEFT_FACTOR& l, G2& msg, const REAL omega = 1.0)
    {
       msg -= omega*l;
    }
    template<typename RIGHT_FACTOR, typename G2>
    void MinimizeRight(const RIGHT_FACTOR& r, G2& msg, const REAL omega = 1.0)
    {
       matrix msgs(i2_,i3_, std::numeric_limits<REAL>::infinity());
       r.min_marginal23(msgs);
       msg -= omega*msgs;
    }

    const INDEX i1_,i2_, i3_;
  };

} // end namespace LP_MP

#endif // LP_MP_SIMPLEX_MARGINALIZATION_MESSAGE_HXX
