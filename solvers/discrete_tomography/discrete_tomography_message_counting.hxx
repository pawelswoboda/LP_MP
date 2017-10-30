#ifndef LP_MP_DT_COUNTING_MESSAGE_HXX
#define LP_MP_DT_COUNTING_MESSAGE_HXX

#include "LP_MP.h"

#include <math.h> 

#include "tropical_convolution.hxx"
//#include "discrete_tomography_algorithms.hxx"
#include "vector.hxx"

namespace LP_MP {
    
  //using MinConv = MinConv<std::function<REAL(INDEX)>,REAL,INDEX>;

  enum class DIRECTION {left,right};

  //template<DIRECTION DR>
  class DiscreteTomographyMessageCounting2 {
   public:
      DiscreteTomographyMessageCounting2(const DIRECTION dr) : dr_(dr) {}

      // do zrobienia: make the member functions static

    //DiscreteTomographyMessageCounting2(INDEX numberOfLabels,INDEX numberOfVarsLeft,INDEX numberOfVarsRight,INDEX SumBound);

    // RIGHT -> TOP Ternary
    // LEFT  -> BOTTOM Ternary

    /* repam.GetFactor()->&f  */
      
    //template<typename RIGHT_FACTOR, typename G1, typename G2>
    //void ReceiveMessageFromRight(RIGHT_FACTOR* const f_right, const G1& repam_right, G2& msg){
    //   MakeRightFactorUniform(f_right, repam_right, msg, 1.0);
    //}

      /*
    template<typename RIGHT_FACTOR, typename MSG_ARRAY, typename ITERATOR>
    static void
    SendMessagesToLeft(const RIGHT_FACTOR& rightFactor, MSG_ARRAY msg_begin, MSG_ARRAY msg_end, ITERATOR omegaIt)
    {
       auto it = msg_begin;
       ++it;
       if(it == msg_end) { // one message
          auto& msg = *msg_begin;
          msg.GetMessageOp().MakeRightFactorUniform(rightFactor, msg, *omegaIt);
       } else {
          auto& first_msg = *msg_begin;
          auto& second_msg = *it;
          assert(first_msg.GetMessageOp().dr_ == DIRECTION::left && second_msg.GetMessageOp().dr_ == DIRECTION::right);
          const REAL omega = *omegaIt + *std::next(omegaIt);
          
          auto factor_copy = rightFactor;
          // first reparametrize half to left, then rest to right, then rest to left
          tensor3 msg_left(factor_copy.no_left_labels(), factor_copy.no_center_left_labels(), factor_copy.left_sum_size(),0.0);
          tensor3 msg_right(factor_copy.no_center_right_labels(), factor_copy.no_right_labels(), factor_copy.right_sum_size(),0.0);
          first_msg.GetMessageOp().MakeRightFactorUniform(factor_copy, msg_left , 0.5*omega);
          first_msg.GetMessageOp().RepamRight(factor_copy, msg_left);
          second_msg.GetMessageOp().MakeRightFactorUniform(factor_copy, msg_right, omega);
          second_msg.GetMessageOp().RepamRight(factor_copy, msg_right);
          first_msg.GetMessageOp().MakeRightFactorUniform(factor_copy, msg_left, omega);

          first_msg -= -msg_left;
          second_msg -= -msg_right;

          ++it;
          assert(it == msg_end);
       }

       //std::cout << "testetset\n";
    }
    */
    template<typename RIGHT_FACTOR, typename G2>
    void SendMessageToLeft(const RIGHT_FACTOR& f_right, G2& msg, const REAL omega){
       //std::cout << "omega = " << omega << "\n";
       MakeRightFactorUniform(f_right, msg, omega);
    }

    //template<typename LEFT_FACTOR, typename G1, typename G3>
    //void SendMessageToRight(LEFT_FACTOR* const f_left, const G1& repam_left, G3& msg, const REAL omega){
    //   MakeLeftFactorUniform(f_left, repam_left, msg, omega);
    //}

    template<typename LEFT_FACTOR, typename G3>
    void ReceiveMessageFromLeft(const LEFT_FACTOR& f_left, G3& msg){
       MakeLeftFactorUniform(f_left, msg, 1.0);
    }

    template<typename RIGHT_FACTOR, typename MSG>
    void MakeRightFactorUniform(const RIGHT_FACTOR& f_right, MSG& msg, const REAL omega) const
    {
       if( dr_ == DIRECTION::left ){
          tensor3<REAL> msg_tmp(f_right.no_left_labels(), f_right.no_center_left_labels(), f_right.left_sum_size());

          f_right.MessageCalculation_Naive_Left(msg_tmp);

          //for(auto it=msg_tmp.begin(); it!=msg_tmp.end(); ++it) {
          //   *it = omega*(*it);
          //}
          msg -= omega*msg_tmp;
       } else if(dr_ == DIRECTION::right) {
          tensor3<REAL> msg_tmp(f_right.no_center_right_labels(), f_right.no_right_labels(), f_right.right_sum_size());

          f_right.MessageCalculation_Naive_Right(msg_tmp);

          //for(auto it=msg_tmp.begin(); it!=msg_tmp.end(); ++it) {
          //   *it = omega*(*it);
          //}
          msg -= omega*msg_tmp;
       } else {
          assert(false);
       } 
    }

    template<typename LEFT_FACTOR, typename MSG>
    void MakeLeftFactorUniform(const LEFT_FACTOR& f_left, MSG& msg, const REAL omega) const
    {
       tensor3<REAL> msg_tmp(f_left.no_left_labels(), f_left.no_right_labels(), f_left.up_sum_size());
       f_left.MessageCalculation_Up(msg_tmp);
       //for(auto it=msg_tmp.begin(); it!=msg_tmp.end(); ++it) {
       //   *it = omega*(*it);
       //}
       msg -= omega*msg_tmp;
    }
    /*------*/

    //template<typename LEFT_FACTOR, typename G1, typename G2>
    //void ReceiveMessageFromLeft(LEFT_FACTOR* const f_left, const G1& repam_left, G2& msg);

    //template<typename LEFT_FACTOR, typename RIGHT_FACTOR, typename G2, typename G3>
    //void SendMessageToLeft(RIGHT_FACTOR* const f_right, const G2& repam_right, G3& msg, const REAL omega);

    /*------*/
      
    template<typename LEFT_FACTOR, typename MSG>
    void RepamLeft(LEFT_FACTOR& l, const MSG& msg) const {
       auto& up = l.up();
       assert(up.size() == msg.size());
       for(INDEX x_l=0; x_l<l.no_left_labels(); ++x_l) {
          for(INDEX x_r=0; x_r<l.no_right_labels(); ++x_r) {
             for(INDEX sum=0; sum<l.up_sum_size(); ++sum) {
                up(x_l,x_r,sum) += normalize( msg(x_l,x_r,sum) );
             }
          }
       } 
       //assert(l.LowerBound() < std::numeric_limits<REAL>::infinity()); 
    }

    template<typename RIGHT_FACTOR, typename MSG>
    void RepamRight(RIGHT_FACTOR& r, const MSG& msg) const {
       //assert(r.LowerBound() < std::numeric_limits<REAL>::infinity()); 
       if(dr_ == DIRECTION::left) {
          auto& left = r.left();
          assert(left.size() == msg.size());
          for(INDEX x_l=0; x_l<r.no_left_labels(); ++x_l) {
             for(INDEX x_r=0; x_r<r.no_center_left_labels(); ++x_r) {
                for(INDEX sum=0; sum<r.left_sum_size(); ++sum) {
                   left(x_l,x_r,sum) += normalize( msg(x_l,x_r,sum) );
                }
             }
          }
       } else if(dr_ == DIRECTION::right) {
          auto& right = r.right();
          assert(right.size() == msg.size());
          for(INDEX x_l=0; x_l<r.no_center_right_labels(); ++x_l) {
             for(INDEX x_r=0; x_r<r.no_right_labels(); ++x_r) {
                for(INDEX sum=0; sum<r.right_sum_size(); ++sum) {
                   right(x_l,x_r,sum) += normalize( msg(x_l,x_r,sum) );
                }
             }
          } 
       } else {
          assert(false);
       }
       //REAL min_val = std::numeric_limits<REAL>::infinity();
       //for(INDEX i=0; i<msg.size(); ++i) { min_val = std::min(min_val, msg[i]); }
       //assert(min_val < std::numeric_limits<REAL>::infinity());
       //assert(*std::min_element(r.up().begin(), r.up().end()) < std::numeric_limits<REAL>::infinity());
       //assert(r.LowerBound() < std::numeric_limits<REAL>::infinity()); 
    }

    /*------*/

    // not used currently, as primal rounding does not give good results
    //template<typename LEFT_FACTOR, typename RIGHT_FACTOR>
    //void ComputeRightFromLeftPrimal(PrimalSolutionStorage::Element left, LEFT_FACTOR* l,
    //                                PrimalSolutionStorage::Element right, RIGHT_FACTOR* r);

    
     private:
    const DIRECTION dr_;
  };



}

#endif // LP_MP_DT_COUNTING_MESSAGE_HXX
