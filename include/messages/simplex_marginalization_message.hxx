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
        ReceiveRestrictedMessageFromRight(const RIGHT_FACTOR& r, G2& msg, typename PrimalSolutionStorage::Element rightPrimal) 
        {
           //std::cout << "before restricted allocation\n";
           vector msgs(i1_,std::numeric_limits<REAL>::infinity());
          //std::vector<REAL> msgs(i1_,std::numeric_limits<REAL>::infinity());
          for(INDEX x2=0; x2<i2_; ++x2) {
            for(INDEX x1=0; x1<i1_; ++x1) {
              if(rightPrimal[x2*i1_ + x1] == unknownState) {
                msgs[x1] = std::min(msgs[x1],r[x2*i1_ + x1]);
              } else if(rightPrimal[x2*i1_ + x1] == true) {
                msgs[x1] = -std::numeric_limits<REAL>::infinity();
              }
            }
          }
          msg -= msgs;
           //std::cout << "before restricted deallocation\n";
        }

    // reparametrize left potential for i-th entry of msg
    // do zrobienia: put strides in here and below
    template<typename G>
    void RepamLeft(G& r, const REAL msg, const INDEX msg_dim)
    {
       assert(!std::isnan(msg));
       if(SUPPORT_INFINITY) {
          if(std::isfinite(msg)) {
             r[msg_dim] += msg;
          } else {
             r[msg_dim] = std::numeric_limits<REAL>::infinity();
          }
          //r[msg_dim] += std::isfinite(msg) ? msg : std::numeric_limits<REAL>::infinity();
          //repamPot[msg_dim] += std::isfinite(msg) ? msg : std::numeric_limits<REAL>::infinity();
       } else {
          r[msg_dim] += msg;
          //repamPot[msg_dim] += msg;
       }
    }
    template<typename A1, typename A2>
    void RepamRight(A1& r, const A2& msgs)
    {
       for(INDEX x1=0; x1<r.dim1(); ++x1) {
          assert(!std::isnan(msgs[x1]));
          //for(INDEX x2=0; x2<r.dim2(); ++x2) {
           if(SUPPORT_INFINITY) {
              if(std::isfinite(msgs[x1])) {
                 r.left(x1) += msgs[x1];
              } else {
                 r.left(x1) = std::numeric_limits<REAL>::infinity();
              }
              //repamPot[x2*i1_ + x1]
              //r.left(x1) += std::isfinite(msgs[x1]) ? msgs[x1] : std::numeric_limits<REAL>::infinity(); 
           } else {
              r.left(x1) += msgs[x1];
              //repamPot[x2*i1_ + x1] += msgs[x1];
           }
         //}
       }
    }
    template<typename G>
    void RepamRight(G& r, const REAL msg, const INDEX dim)
    {
       assert(!std::isnan(msg));
       if(SUPPORT_INFINITY) {
          if(std::isfinite(msg)) {
             r.left(dim) += msg;
          } else {
             r.left(dim) = std::numeric_limits<REAL>::infinity();
          }
          //REAL msgn = std::isfinite(msg) ? msg : std::numeric_limits<REAL>::infinity();
          //r.left(dim) += msgn;
          //for(INDEX x2=0; x2<i2_; ++x2) {
             //repamPot[x2*i1_ + dim] += msgn;
          //}
       } else {
          r.left(dim) += msg;
          //for(INDEX x2=0; x2<i2_; ++x2) {
             //repamPot[x2*i1_ + dim] += msg;
          //}
       }
    }

    template<bool ENABLE = TYPE == MessageSendingType::SRMP, typename LEFT_FACTOR, typename RIGHT_FACTOR>
    //typename std::enable_if<ENABLE,void>::type
    void
    ComputeRightFromLeftPrimal(const typename PrimalSolutionStorage::Element left, const LEFT_FACTOR& l, typename PrimalSolutionStorage::Element right, const RIGHT_FACTOR& r)
    {
      for(INDEX x1=0; x1<i1_; ++x1) {
        if(left[x1] == false) {
          for(INDEX x2=0; x2<i2_; ++x2) {
            right[x2*i1_ + x1] = false;
          }
        } else if(left[x1] == true) {
          INDEX no_false = 0;
          for(INDEX x2=0; x2<i2_; ++x2) {
            if(right[x2*i1_ + x1] == false) {
              ++no_false;
            }
          }
          if(no_false == i1_-1) {
            for(INDEX x2=0; x2<i2_; ++x2) {
              if(right[x2*i1_ + x1] == unknownState) {
                right[x2*i1_ + x1] = true;
              }
            }
          }
        } else {
          assert(false); // unary should have been labelled
        }
      }
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
             msgs[x1] = std::min(msgs[x1],omega*r(x1,x2));
             //msgs[x1] = std::min(msgs[x1],omega*r[x2*i1_ + x1]);
          }
       }
       msg -= msgs;
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
        ReceiveRestrictedMessageFromRight(const RIGHT_FACTOR& r, G2& msg, typename PrimalSolutionStorage::Element rightPrimal) 
        {
           vector msgs(i2_,std::numeric_limits<REAL>::infinity());
          //std::vector<REAL> msgs(i2_,std::numeric_limits<REAL>::infinity());
          for(INDEX x2=0; x2<i2_; ++x2) {
            for(INDEX x1=0; x1<i1_; ++x1) {
              if(rightPrimal[x2*i1_ + x1] == unknownState) {
                msgs[x2] = std::min(msgs[x2],r[x2*i1_ + x1]);
              } else if(rightPrimal[x2*i1_ + x1] == true) {
                msgs[x2] = -std::numeric_limits<REAL>::infinity();
              }
            }
          }
          msg -= msgs;
          //assert(false);
        }

    // reparametrize left potential for i-th entry of msg
    // do zrobienia: put strides in here and below
    template<typename G>
    void RepamLeft(G& l, const REAL msg, const INDEX msg_dim)
    {
       assert(!std::isnan(msg));
       if(SUPPORT_INFINITY) {
          if(std::isfinite(msg)) {
             l[msg_dim] += msg;
          } else {
             l[msg_dim] = std::numeric_limits<REAL>::infinity();
          }
          //repamPot[msg_dim] += std::isfinite(msg) ? msg : std::numeric_limits<REAL>::infinity();
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
            if(std::isfinite(msgs[x2])) {
               r.right(x2) += msgs[x2];
            } else {
               r.right(x2) = std::numeric_limits<REAL>::infinity();
            }
            //const REAL msgn = std::isfinite(msgs[x2]) ? msgs[x2] : std::numeric_limits<REAL>::infinity();
            //r.right(x2) += msgn;
            //for(INDEX x1=0; x1<i1_; ++x1) {
            //   repamPot[x2*i1_ + x1] += msgn;
            //}
         } else {
            //for(INDEX x1=0; x1<i1_; ++x1) {
            //   repamPot[x2*i1_ + x1] += msgs[x2];
            //}
            r.right(x2) += msgs[x2];
         }
      }
    }
    template<typename G>
    void RepamRight(G& r, const REAL msg, const INDEX dim)
    {
       assert(!std::isnan(msg));
       if(SUPPORT_INFINITY) {
          if(std::isfinite(msg)) {
             r.right(dim) += msg;
          } else {
             r.right(dim) = std::numeric_limits<REAL>::infinity();
          }
          //REAL msgn = std::isfinite(msg) ? msg : std::numeric_limits<REAL>::infinity();
          //r.right(dim) += msgn;
          //for(INDEX x1=0; x1<i1_; ++x1) {
          //   repamPot[dim*i1_ + x1] += msgn;
          //}
       } else {
          //for(INDEX x1=0; x1<i1_; ++x1) {
          //   repamPot[dim*i1_ + x1] += msg;
          //}
          r.right(dim) += msg;
       }
    }

    template<bool ENABLE = TYPE == MessageSendingType::SRMP, typename LEFT_FACTOR, typename RIGHT_FACTOR>
    //typename std::enable_if<ENABLE,void>::type
    void
    ComputeRightFromLeftPrimal(const typename PrimalSolutionStorage::Element left, const LEFT_FACTOR& l, typename PrimalSolutionStorage::Element right, const RIGHT_FACTOR& r)
    {
      for(INDEX x2=0; x2<i2_; ++x2) {
        if(left[x2] == false) {
          for(INDEX x1=0; x1<i1_; ++x1) {
            right[x2*i1_ + x1] = false;
          }
        } else if(left[x2] == true) {
          INDEX no_false = 0;
          for(INDEX x1=0; x1<i1_; ++x1) {
            if(right[x2*i1_ + x1] == false) {
              ++no_false;
            }
          }
          if(no_false == i2_-1) {
            for(INDEX x1=0; x1<i1_; ++x1) {
              if(right[x2*i1_ + x1] == unknownState) {
                right[x2*i1_ + x1] = true;
              }
            }
          }
        } else {
          assert(false); // unary should have been labelled
        }
      }
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

    const INDEX i1_,i2_;
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
    template<typename G>
    void RepamRight(G& r, const REAL msg, const INDEX dim)
    {
      const INDEX x1 = dim/i2_;
      const INDEX x2 = dim%i2_;
      r.msg12(x1,x2) += normalize( msg );
    }

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
      for(INDEX x2=0; x2<i2_; ++x2) {
        for(INDEX x1=0; x1<i1_; ++x1) {
          msg[x2*i1_ + x1] -= omega*l[x2*i1_ + x1];
        }
      }
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
    template<typename G>
    void RepamLeft(G& repamPot, const REAL msg, const INDEX msg_dim)
    {
      REAL msgn;
      if( std::isfinite(msg) ){ msgn = msg; }
      else{ msgn = std::numeric_limits<REAL>::infinity(); }
      const INDEX x1 = msg_dim/i3_;
      const INDEX x2 = msg_dim%i3_;
      repamPot(x1,x2) += msgn;
    }
    template<typename A1, typename A2>
    void RepamRight(A1& repamPot, const A2& msgs)
    {
      // do zrobienia: possibly use counter
      for(INDEX x3 = 0; x3<i3_; ++x3) {
        for(INDEX x2=0; x2<i2_; ++x2) {
          for(INDEX x1=0; x1<i1_; ++x1) {
            REAL msgn;
            if( std::isfinite(msgs[x3*i1_ + x1]) ){ msgn = msgs[x3*i1_ + x1]; }
            else{ msgn = std::numeric_limits<REAL>::infinity(); }
            repamPot[x3*i1_*i2_ + x2*i1_ + x1] += msgn;
          }
        }
      }
    }
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
      for(INDEX x3=0; x3<i3_; ++x3) {
        for(INDEX x1=0; x1<i1_; ++x1) {
          msg[x3*i1_ + x1] -= omega*l[x3*i1_ + x1];
        }
      }
    }
    template<typename RIGHT_FACTOR, typename G2>
    void MinimizeRight(const RIGHT_FACTOR& r, G2& msg, const REAL omega = 1.0)
    {
      std::vector<REAL> msgs(i1_*i3_, std::numeric_limits<REAL>::infinity());
      for(INDEX x3=0; x3<i3_; ++x3) {
        for(INDEX x2=0; x2<i2_; ++x2) {
          for(INDEX x1=0; x1<i1_; ++x1) {
            msgs[x3*i1_ + x1] = std::min(msgs[x3*i1_ + x1],omega*r[x3*i1_*i2_ + x2*i1_ + x1]);
          }
        }
      }
      assert(false);
      //msg -= msgs;
    }

    const INDEX i1_,i2_, i3_;
  };

  template<MessageSendingType TYPE,  bool PROPAGATE_PRIMAL_TO_LEFT = false, bool PROPAGATE_PRIMAL_TO_RIGHT = false> 
  class PairwiseTripletMessage23 {
    public:
      PairwiseTripletMessage23(const INDEX i1, const INDEX i2, const INDEX i3) : i1_(i1), i2_(i2), i3_(i3) {} // the pairwise factor size
      // standard functions which take all possible arguments, to be replaced with minimal ones
      template<typename RIGHT_FACTOR, typename G1, typename G2, bool ENABLE = TYPE == MessageSendingType::SRMP>
        typename std::enable_if<ENABLE,void>::type
        ReceiveMessageFromRight(RIGHT_FACTOR* const r, const G1& rightPot, G2& msg) 
        {
          MinimizeRight(r,rightPot,msg); 
        }
      template<typename LEFT_FACTOR, typename G1, typename G2, bool ENABLE = TYPE == MessageSendingType::MPLP>
        typename std::enable_if<ENABLE,void>::type 
        ReceiveMessageFromLeft(LEFT_FACTOR* l, const G1& leftPot, G2& msg) 
        { 
          MaximizeLeft(l,leftPot,msg); 
        }
      template<typename LEFT_FACTOR, typename G1, typename G3, bool ENABLE = TYPE == MessageSendingType::SRMP>
        typename std::enable_if<ENABLE,void>::type
        SendMessageToRight(LEFT_FACTOR* const l, const G1& leftPot, G3& msg, const REAL omega)
        { 
          MaximizeLeft(l,leftPot,msg,omega); 
        }
      template<typename RIGHT_FACTOR, typename G2, typename G3, bool ENABLE = TYPE == MessageSendingType::MPLP>
        typename std::enable_if<ENABLE,void>::type
        SendMessageToLeft(RIGHT_FACTOR* const r, const G2& rightPot, G3& msg, const REAL omega) 
        { 
          MinimizeRight(r,rightPot,msg,omega); 
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
    template<typename G>
    void RepamLeft(G& repamPot, const REAL msg, const INDEX msg_dim)
    {
      REAL msgn;
      if( std::isfinite(msg) ){ msgn = msg; }
      else{ msgn = std::numeric_limits<REAL>::infinity(); }
      repamPot[msg_dim] += msgn;
    }
    template<typename A1, typename A2>
    void RepamRight(A1& repamPot, const A2& msgs)
    {
      // do zrobienia: possibly use counter
      for(INDEX x3 = 0; x3<i3_; ++x3) {
        for(INDEX x2=0; x2<i2_; ++x2) {
          REAL msgn;
          if( std::isfinite(msgs[x3*i2_ + x2]) ){ msgn = msgs[x3*i2_ + x2]; }
          else{ msgn = std::numeric_limits<REAL>::infinity(); }
          for(INDEX x1=0; x1<i1_; ++x1) {
            repamPot[x3*i1_*i2_ + x2*i1_ + x1] += msgn;
          }
        }
      }
    }
    template<typename G>
    void RepamRight(G& repamPot, const REAL msg, const INDEX dim)
    {
      REAL msgn;
      const INDEX x2 = dim % i2_;
      const INDEX x3 = dim / i2_;
      if( std::isfinite(msg) ){ msgn = msg; }
      else{ msgn = std::numeric_limits<REAL>::infinity(); }
      for(INDEX x1 = 0; x1<i1_; ++x1) {
        repamPot[x3*i2_*i1_ + x2*i1_ + x1] += msgn;
      }
    }

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
    template<typename LEFT_FACTOR, typename G1, typename G2>
    void MaximizeLeft(LEFT_FACTOR* const l, const G1& leftPot, G2& msg, const REAL omega = 1.0)
    {
      for(INDEX x3=0; x3<i3_; ++x3) {
        for(INDEX x2=0; x2<i2_; ++x2) {
          msg[x3*i2_ + x2] -= omega*leftPot[x3*i2_ + x2];
        }
      }
    }
    template<typename RIGHT_FACTOR, typename G1, typename G2>
    void MinimizeRight(RIGHT_FACTOR* r, const G1& rightPot, G2& msg, const REAL omega = 1.0)
    {
      for(INDEX x3=0; x3<i3_; ++x3) {
        for(INDEX x2=0; x2<i2_; ++x2) {
          REAL msg_val = std::numeric_limits<REAL>::infinity();
          for(INDEX x1=0; x1<i1_; ++x1) {
            msg_val = std::min(msg_val,omega*rightPot[x3*i1_*i2_ + x2*i1_ + x1]);
          }
          msg[x3*i2_ + x2] -= msg_val;
        }
      }
    }

    const INDEX i1_,i2_, i3_;
  };



























































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
      REAL msgn;
      if( std::isfinite(msg) ){ msgn = msg; }
      else{ msgn = std::numeric_limits<REAL>::infinity(); }
      loopLeft_.loop(msg_dim, [&](const INDEX i) { repamPot[leftStride_[i]] += msgn; });
    }
    template<typename G>
    void RepamRight(G& repamPot, const REAL msg, const INDEX dim)
    {
      REAL msgn;
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
   
   
    /*
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
        LinExpr lhs = lp->CreateLinExpr() + 0.0;
        LinExpr rhs = lp->CreateLinExpr() + 0.0;
        bool least = false;
        for(auto l=0; l<leftIndices[i].size(); l++){
          if(lp->IsLeftObjective(leftIndices[i][l])){
            lhs += lp->GetLeftVariable(leftIndices[i][l]);
            least = true;
          }
        }
        for(auto r=0; r<rightIndices[i].size(); r++){
          if(lp->IsRightObjective(rightIndices[i][r])){
            rhs += lp->GetRightVariable(rightIndices[i][r]);
            least = true;
          }
        }
        if(least){
          lp->addLinearEquality(lhs,rhs);
        }
      }
    }
    */
  
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
	       if(primal[ full_idx ] == unknownState) { // i.e. may be true or unknown. do zrobienia: if true, then go back, set delta to infinity for all other coordinates and this coordinate let be -infinity
          delta = std::min(delta, pot[ full_idx ]);
          } else if(primal[ full_idx ] == true) {
          delta = -std::numeric_limits<REAL>::infinity();
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
