#ifndef LP_MP_DT_COUNTING_MESSAGE_HXX
#define LP_MP_DT_COUNTING_MESSAGE_HXX

#include "LP_MP.h"

#include <math.h> 

#include "minConv.hxx"
#include "discrete_tomography_algorithms.hxx"
#include "vector.hxx"

namespace LP_MP {
    
  //using MinConv = MinConv<std::function<REAL(INDEX)>,REAL,INDEX>;

  enum class DIRECTION {left,right};

  template<DIRECTION DR>
  class DiscreteTomographyMessageCounting2 {
   public:

      // do zrobienia: make the member functions static

    //DiscreteTomographyMessageCounting2(INDEX numberOfLabels,INDEX numberOfVarsLeft,INDEX numberOfVarsRight,INDEX SumBound);

    // RIGHT -> TOP Ternary
    // LEFT  -> BOTTOM Ternary

    /* repam.GetFactor()->&f  */
      
    //template<typename RIGHT_FACTOR, typename G1, typename G2>
    //void ReceiveMessageFromRight(RIGHT_FACTOR* const f_right, const G1& repam_right, G2& msg){
    //   MakeRightFactorUniform(f_right, repam_right, msg, 1.0);
    //}

    template<typename RIGHT_FACTOR, typename G2>
    void SendMessageToLeft(const RIGHT_FACTOR& f_right, G2& msg, const REAL omega){
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
    void MakeRightFactorUniform(const RIGHT_FACTOR& f_right, MSG& msg, const REAL omega)
    {
       if( DR == DIRECTION::left ){
          tensor3 msg_tmp(f_right.no_labels(), f_right.no_labels(), f_right.left_sum_size());

          f_right.MessageCalculation_Naive_Left(msg_tmp);

          //for(auto it=msg_tmp.begin(); it!=msg_tmp.end(); ++it) {
          //   *it = omega*(*it);
          //}
          msg -= omega*msg_tmp;
       } else if(DR == DIRECTION::right) {
          tensor3 msg_tmp(f_right.no_labels(), f_right.no_labels(), f_right.right_sum_size());

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
    void MakeLeftFactorUniform(const LEFT_FACTOR& f_left, MSG& msg, const REAL omega)
    {
       tensor3 msg_tmp(f_left.no_labels(), f_left.no_labels(), f_left.up_sum_size());
       f_left.MessageCalculation_Naive_Up(msg_tmp);
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
    void RepamLeft(LEFT_FACTOR& l, const MSG msg){
       auto& up = l.up();
       assert(up.size() == msg.size());
       for(INDEX x_l=0; x_l<l.no_labels(); ++x_l) {
          for(INDEX x_r=0; x_r<l.no_labels(); ++x_r) {
             for(INDEX sum=0; sum<l.up_sum_size(); ++sum) {
                up(x_l,x_r,sum) += normalize( msg(x_l,x_r,sum) );
             }
          }
       } 
       assert(l.LowerBound() < std::numeric_limits<REAL>::infinity()); 
    }

    template<typename RIGHT_FACTOR, typename MSG>
    void RepamRight(RIGHT_FACTOR& r, const MSG msg){
       assert(r.LowerBound() < std::numeric_limits<REAL>::infinity()); 
       if(DR == DIRECTION::left) {
          auto& left = r.left();
          assert(left.size() == msg.size());
          for(INDEX x_l=0; x_l<r.no_labels(); ++x_l) {
             for(INDEX x_r=0; x_r<r.no_labels(); ++x_r) {
                for(INDEX sum=0; sum<r.left_sum_size(); ++sum) {
                   left(x_l,x_r,sum) += normalize( msg(x_l,x_r,sum) );
                }
             }
          }
       } else if(DR == DIRECTION::right) {
          auto& right = r.right();
          assert(right.size() == msg.size());
          for(INDEX x_l=0; x_l<r.no_labels(); ++x_l) {
             for(INDEX x_r=0; x_r<r.no_labels(); ++x_r) {
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

    
  };


  template<DIRECTION DR>
  class DiscreteTomographyMessageCounting{

  public:

    DiscreteTomographyMessageCounting(INDEX numberOfLabels,INDEX numberOfVarsLeft,INDEX numberOfVarsRight,INDEX SumBound);

    // RIGHT -> TOP Ternary
    // LEFT  -> BOTTOM Ternary

    template<class LEFT_FACTOR_TYPE,class RIGHT_FACTOR_TYPE>
    void CreateConstraints(LpInterfaceAdapter* lp,LEFT_FACTOR_TYPE* LeftFactor,RIGHT_FACTOR_TYPE* RightFactor) const;
    
    /* repam.GetFactor()->&f  */
      
    //template<typename RIGHT_FACTOR, typename G1, typename G2>
    //void ReceiveMessageFromRight(RIGHT_FACTOR* const f_right, const G1& repam_right, G2& msg){
    //   MakeRightFactorUniform(f_right, repam_right, msg, 1.0);
    //}

    template<typename RIGHT_FACTOR, typename G1, typename G2>
    void SendMessageToLeft(RIGHT_FACTOR* const f_right, const G1& repam_right, G2& msg, const REAL omega){
       MakeRightFactorUniform(f_right, repam_right, msg, omega);
    }

    //template<typename LEFT_FACTOR, typename G1, typename G3>
    //void SendMessageToRight(LEFT_FACTOR* const f_left, const G1& repam_left, G3& msg, const REAL omega){
    //   MakeLeftFactorUniform(f_left, repam_left, msg, omega);
    //}

    template<typename LEFT_FACTOR, typename G1, typename G3>
    void ReceiveMessageFromLeft(LEFT_FACTOR* const f_left, const G1& repam_left, G3& msg){
       MakeLeftFactorUniform(f_left, repam_left, msg, 1.0);
    }

    template<typename RIGHT_FACTOR, typename REPAM_ARRAY, typename MSG>
    void MakeRightFactorUniform(RIGHT_FACTOR* const f_right, const REPAM_ARRAY& repam_right, MSG& msg, const REAL omega);

    template<typename LEFT_FACTOR, typename REPAM_ARRAY, typename MSG>
    void MakeLeftFactorUniform(LEFT_FACTOR* const f_left, const REPAM_ARRAY& repam_left, MSG& msg, const REAL omega);
    /*------*/
      
    //template<typename LEFT_FACTOR, typename G1, typename G2>
    //void ReceiveMessageFromLeft(LEFT_FACTOR* const f_left, const G1& repam_left, G2& msg);

    //template<typename LEFT_FACTOR, typename RIGHT_FACTOR, typename G2, typename G3>
    //void SendMessageToLeft(RIGHT_FACTOR* const f_right, const G2& repam_right, G3& msg, const REAL omega);

    /*------*/
      
    template<typename G>
    void RepamLeft(G& repam, const REAL msg, const INDEX msg_dim);
    
    template<typename G>
    void RepamRight(G& repam, const REAL msg, const INDEX msg_dim);

    /*------*/

    // not used currently, as primal rounding does not give good results
    //template<typename LEFT_FACTOR, typename RIGHT_FACTOR>
    //void ComputeRightFromLeftPrimal(PrimalSolutionStorage::Element left, LEFT_FACTOR* l,
    //                                PrimalSolutionStorage::Element right, RIGHT_FACTOR* r);

    INDEX size() const {
       if(DR == DIRECTION::left) {
          return pow(numberOfLabels_,2)*std::min(sumBound_, numberOfVarsLeft_*(numberOfLabels_-1)+1);
       } else if(DR == DIRECTION::right) {
          return pow(numberOfLabels_,2)*std::min(sumBound_, numberOfVarsRight_*(numberOfLabels_-1)+1);
       } else {
          assert(false);
       }
    }
    
  private:
    const INDEX numberOfLabels_,numberOfVarsLeft_,numberOfVarsRight_,sumBound_;
    INDEX upSize_,leftSize_,rightSize_,regSize_;
  };

  template<DIRECTION DR>
  DiscreteTomographyMessageCounting<DR>::DiscreteTomographyMessageCounting(INDEX numberOfLabels,INDEX numberOfVarsLeft,INDEX numberOfVarsRight,INDEX SumBound)
    : numberOfLabels_(numberOfLabels),numberOfVarsLeft_(numberOfVarsLeft),numberOfVarsRight_(numberOfVarsRight),sumBound_(SumBound) {
    assert(numberOfLabels_ > 1);
    assert(numberOfVarsLeft_ > 0);
    assert(numberOfVarsRight_ > 0);
    assert(sumBound_ > 0);
    
    upSize_ = pow(numberOfLabels_,2)*std::min(((numberOfVarsLeft_+numberOfVarsRight_)*(numberOfLabels_-1)+1),SumBound);
    leftSize_ = pow(numberOfLabels_,2)*std::min((numberOfVarsLeft_*(numberOfLabels_-1)+1),SumBound);
    rightSize_ = pow(numberOfLabels_,2)*std::min((numberOfVarsRight_*(numberOfLabels_-1)+1),SumBound);
    
    regSize_ = pow(numberOfLabels_,2);
  }

  /*
  template<DIRECTION DR>
  template<typename LEFT_FACTOR, typename G1, typename G3>
  void DiscreteTomographyMessageCounting<DR>::SendMessageToRight(LEFT_FACTOR* const f_left, const G1& repam_left, G3& msg, const REAL omega){    
     MakeLeftFactorUniform(f_left, repam_left, msg, omega);
  }

  template<DIRECTION DR>
  template<typename RIGHT_FACTOR, typename G1, typename G2>
  void DiscreteTomographyMessageCounting<DR>::ReceiveMessageFromRight(RIGHT_FACTOR* const f_right, const G1& repam_right, G2& msg){
      MakeRightFactorUniform(f_right, repam_right, msg, 1.0);
  }
  */

  template<DIRECTION DR>
  template<class LEFT_FACTOR_TYPE,class RIGHT_FACTOR_TYPE>
  void DiscreteTomographyMessageCounting<DR>::CreateConstraints(LpInterfaceAdapter* lp,LEFT_FACTOR_TYPE* LeftFactor,RIGHT_FACTOR_TYPE* RightFactor) const {
    INDEX upSize = (*RightFactor).getSize(DiscreteTomographyFactorCounting::NODE::up);
    INDEX leftSize = (*RightFactor).getSize(DiscreteTomographyFactorCounting::NODE::left);
    INDEX rightSize = (*RightFactor).getSize(DiscreteTomographyFactorCounting::NODE::right);

    INDEX noVars = (*LeftFactor).getSize(DiscreteTomographyFactorCounting::NODE::up);
    auto CoupleConstraints = [&](INDEX offset){
      for(INDEX i=0;i<noVars;i++){
        LinExpr lhs = lp->CreateLinExpr();
        LinExpr rhs = lp->CreateLinExpr();
        lhs += lp->GetLeftVariable(i);
        rhs += lp->GetRightVariable(offset + i);
        lp->addLinearEquality(lhs,rhs);
      }
    };
    if( DR == DIRECTION::left ){
      assert((*LeftFactor).getSize(DiscreteTomographyFactorCounting::NODE::up) == leftSize );
      CoupleConstraints(upSize);
    } else {
      assert((*LeftFactor).getSize(DiscreteTomographyFactorCounting::NODE::up) == rightSize );
      CoupleConstraints(upSize+leftSize);
    }
  }

  template<DIRECTION DR>
  template<typename RIGHT_FACTOR, typename REPAM_ARRAY, typename MSG>
  void DiscreteTomographyMessageCounting<DR>::MakeRightFactorUniform(RIGHT_FACTOR* const f_right, const REPAM_ARRAY& repam_right, MSG& msg, const REAL omega)
  {
    assert(repam_right.size() == ((*f_right).getSize(DiscreteTomographyFactorCounting::NODE::up) +
                                  (*f_right).getSize(DiscreteTomographyFactorCounting::NODE::left) +
                                  (*f_right).getSize(DiscreteTomographyFactorCounting::NODE::right) +
                                  (*f_right).getSize(DiscreteTomographyFactorCounting::NODE::reg)));

    INDEX up_size = (*f_right).getSize(DiscreteTomographyFactorCounting::NODE::up)/pow(numberOfLabels_,2);
    INDEX right_size = (*f_right).getSize(DiscreteTomographyFactorCounting::NODE::right)/pow(numberOfLabels_,2);
    INDEX left_size = (*f_right).getSize(DiscreteTomographyFactorCounting::NODE::left)/pow(numberOfLabels_,2);
     
    if( DR == DIRECTION::left ){
      assert(msg.size() == (*f_right).getSize(DiscreteTomographyFactorCounting::NODE::left));

      std::vector<REAL> msg_v(left_size*pow(numberOfLabels_,2),std::numeric_limits<REAL>::infinity());
      if( DiscreteTomo::AlgorithmThreshold < msg_v.size() ){
        DiscreteTomo::MessageCalculation_MinConv_Left(f_right, repam_right, msg_v,numberOfLabels_);
      } else {
        DiscreteTomo::MessageCalculation_Naive_Left(f_right, repam_right, msg_v,numberOfLabels_);
      }
      for(INDEX i=0;i<msg_v.size();i++){
        assert(msg_v[i] > -eps);
        msg[i] -= omega*msg_v[i];
      }
       
    }
    else{
      assert(msg.size() == (*f_right).getSize(DiscreteTomographyFactorCounting::NODE::right));

      std::vector<REAL> msg_v(right_size*pow(numberOfLabels_,2),std::numeric_limits<REAL>::infinity());
      if( DiscreteTomo::AlgorithmThreshold < msg_v.size() ){
        DiscreteTomo::MessageCalculation_MinConv_Right(f_right, repam_right, msg_v,numberOfLabels_);
      } else {
        DiscreteTomo::MessageCalculation_Naive_Right(f_right, repam_right, msg_v,numberOfLabels_);
      }
      for(INDEX i=0;i<msg_v.size();i++){
        assert(msg_v[i] > -eps);
        msg[i] -= omega*msg_v[i];
      }
       
    }      
  }

  template<DIRECTION DR>
  template<typename LEFT_FACTOR, typename REPAM_ARRAY, typename MSG>
  void DiscreteTomographyMessageCounting<DR>::MakeLeftFactorUniform(LEFT_FACTOR* const f_left, const REPAM_ARRAY& repam_left, MSG& msg, const REAL omega)
  {
    assert(msg.size() == (*f_left).getSize(DiscreteTomographyFactorCounting::NODE::up));
    assert(repam_left.size() == ((*f_left).getSize(DiscreteTomographyFactorCounting::NODE::up) +
                                 (*f_left).getSize(DiscreteTomographyFactorCounting::NODE::left) +
                                 (*f_left).getSize(DiscreteTomographyFactorCounting::NODE::right) +
                                 (*f_left).getSize(DiscreteTomographyFactorCounting::NODE::reg)));

    INDEX up_size = (*f_left).getSize(DiscreteTomographyFactorCounting::NODE::up)/pow(numberOfLabels_,2);
    INDEX right_size = (*f_left).getSize(DiscreteTomographyFactorCounting::NODE::right)/pow(numberOfLabels_,2);
    INDEX left_size = (*f_left).getSize(DiscreteTomographyFactorCounting::NODE::left)/pow(numberOfLabels_,2);
    
    std::vector<REAL> msg_v(up_size*pow(numberOfLabels_,2),std::numeric_limits<REAL>::infinity());
    if( DiscreteTomo::AlgorithmThreshold < msg_v.size() ){
      DiscreteTomo::MessageCalculation_MinConv_Up(f_left, repam_left, msg_v,numberOfLabels_);
    } else {
      DiscreteTomo::MessageCalculation_Naive_Up(f_left, repam_left, msg_v,numberOfLabels_);
    }
    
    assert(msg.size() == (*f_left).getSize(DiscreteTomographyFactorCounting::NODE::up));
    assert(repam_left.size() == ((*f_left).getSize(DiscreteTomographyFactorCounting::NODE::up) +
                                 (*f_left).getSize(DiscreteTomographyFactorCounting::NODE::left) +
                                 (*f_left).getSize(DiscreteTomographyFactorCounting::NODE::right) +
                                 (*f_left).getSize(DiscreteTomographyFactorCounting::NODE::reg)));
      
    for(INDEX i=0;i<msg_v.size();i++){
      assert(msg_v[i] > -eps);
      msg[i] -= omega*msg_v[i];
    }
 
  }

  template<DIRECTION DR>
  template<typename G>
  void DiscreteTomographyMessageCounting<DR>::RepamLeft(G& repam, const REAL msg, const INDEX msg_dim){
    
    assert(!std::isnan(msg));
    auto f = repam.GetFactor();
    assert( repam.size() == (f->getSize(DiscreteTomographyFactorCounting::NODE::left) +
                             f->getSize(DiscreteTomographyFactorCounting::NODE::right) +
                             f->getSize(DiscreteTomographyFactorCounting::NODE::up) +
                             f->getSize(DiscreteTomographyFactorCounting::NODE::reg)));

     assert(msg_dim < f->getSize(DiscreteTomographyFactorCounting::NODE::up));
     assert(repam[msg_dim] > -std::numeric_limits<REAL>::max());
     if( std::isfinite(msg) ){ 
        repam[msg_dim] += msg; 
        assert(repam[msg_dim] > -std::numeric_limits<REAL>::max() );
     }
     else{ repam[msg_dim] = std::numeric_limits<REAL>::infinity(); }

  }

  template<DIRECTION DR>    
  template<typename G>
  void DiscreteTomographyMessageCounting<DR>::RepamRight(G& repam, const REAL msg, const INDEX msg_dim){

     assert(!std::isnan(msg));
     auto f = repam.GetFactor();
     assert( repam.size() == (f->getSize(DiscreteTomographyFactorCounting::NODE::left) +
              f->getSize(DiscreteTomographyFactorCounting::NODE::right) +
              f->getSize(DiscreteTomographyFactorCounting::NODE::up) +
              f->getSize(DiscreteTomographyFactorCounting::NODE::reg)));

     if( DR == DIRECTION::left ){
        assert(msg_dim < f->getSize(DiscreteTomographyFactorCounting::NODE::left));
        assert(repam[f->getSize(DiscreteTomographyFactorCounting::NODE::up) + msg_dim] > -std::numeric_limits<REAL>::max());
        if( std::isfinite(msg) ){
           repam[f->getSize(DiscreteTomographyFactorCounting::NODE::up) + msg_dim] +=  msg; 
           assert(repam[f->getSize(DiscreteTomographyFactorCounting::NODE::up) + msg_dim] > -std::numeric_limits<REAL>::max());
        }
        else{
           repam[f->getSize(DiscreteTomographyFactorCounting::NODE::up) + msg_dim] = std::numeric_limits<REAL>::infinity();
        }
     }
     else{
        assert(msg_dim < f->getSize(DiscreteTomographyFactorCounting::NODE::right));
        assert(repam[f->getSize(DiscreteTomographyFactorCounting::NODE::up) + f->getSize(DiscreteTomographyFactorCounting::NODE::left) + msg_dim] > -std::numeric_limits<REAL>::max());
        if( std::isfinite(msg) ){      
           repam[f->getSize(DiscreteTomographyFactorCounting::NODE::up) + f->getSize(DiscreteTomographyFactorCounting::NODE::left) + msg_dim] += msg; 
           assert(repam[f->getSize(DiscreteTomographyFactorCounting::NODE::up) + f->getSize(DiscreteTomographyFactorCounting::NODE::left) + msg_dim] > -std::numeric_limits<REAL>::max());
        }
        else{
           repam[f->getSize(DiscreteTomographyFactorCounting::NODE::up) + f->getSize(DiscreteTomographyFactorCounting::NODE::left) + msg_dim] = std::numeric_limits<REAL>::infinity();
        }
     }
  }

  /* 
  template<DIRECTION DR>  
  template<typename LEFT_FACTOR, typename RIGHT_FACTOR>
  void DiscreteTomographyMessageCounting<DR>::ComputeRightFromLeftPrimal(PrimalSolutionStorage::Element left, LEFT_FACTOR* leftFactor,
                                                                         PrimalSolutionStorage::Element right, RIGHT_FACTOR* rightFactor){
    assert(rightFactor->getSize(DiscreteTomographyFactorCounting::NODE::up) == upSize_);
    assert(rightFactor->getSize(DiscreteTomographyFactorCounting::NODE::left) == leftSize_);
    assert(rightFactor->getSize(DiscreteTomographyFactorCounting::NODE::right) == rightSize_);
    assert(rightFactor->getSize(DiscreteTomographyFactorCounting::NODE::reg) == regSize_);
    
    auto copyPrimal = [&](INDEX s,INDEX t){
      INDEX noTrue = 0;
      INDEX noUnkwn = 0;
      INDEX opt = 0;
      for(INDEX i=0;i<s;i++){
        if( left[i] == true ){
          noTrue++;
          assert(right[upSize_ + t + i] != false); // do not overwrite results!
        }
        if( left[i] == unknownState ){ noUnkwn++; }
        right[upSize_ + t + i] = left[i];
      }
      assert(noTrue <=1);
      assert(noTrue != 1 || noUnkwn == 0);
      assert(noUnkwn != 0 || noTrue == 1 );
      assert(noUnkwn != 0 || noTrue != 0);
    };
     
    if( DR == DIRECTION::left ){
      assert( leftFactor->getSize(DiscreteTomographyFactorCounting::NODE::up) == leftSize_ );
      copyPrimal(leftSize_,0);
    }
    else{
      assert( leftFactor->getSize(DiscreteTomographyFactorCounting::NODE::up) == rightSize_ );
      copyPrimal(rightSize_,leftSize_);
    }

  }
  */

}

#endif // LP_MP_DT_COUNTING_MESSAGE_HXX
