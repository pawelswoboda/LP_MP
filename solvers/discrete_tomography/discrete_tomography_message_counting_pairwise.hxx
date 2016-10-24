#ifndef LP_MP_DT_COUNTING_PAIRWISE_MESSAGE_HXX
#define LP_MP_DT_COUNTING_PAIRWISE_MESSAGE_HXX

#include "LP_MP.h"
#include "minConv.hxx"
#include <math.h> 

namespace LP_MP {
    
  //using MinConv = MinConv<REAL,INDEX>;

  class DiscreteTomographyMessageCountingPairwise{

  public:

    DiscreteTomographyMessageCountingPairwise(INDEX numberOfLabels,INDEX numberOfVarsLeft,INDEX numberOfVarsRight,INDEX SumBound);

    // RIGHT -> TOP Ternary
    // LEFT  -> Pairwise

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


    template<typename LEFT_FACTOR, typename REPAM_ARRAY, typename MSG>
    void MakeLeftFactorUniform(LEFT_FACTOR* f_left, const REPAM_ARRAY& repam_left, MSG& msg, const REAL omega);

    template<typename RIGHT_FACTOR, typename REPAM_ARRAY, typename MSG>
    void MakeRightFactorUniform(RIGHT_FACTOR* f_right, const REPAM_ARRAY& repam_right, MSG& msg, const REAL omega);

    // Update repam with message for each factor ?
      
    template<typename G>
    void RepamLeft(G& repam, const REAL msg, const INDEX msg_dim);
    
    template<typename G>
    void RepamRight(G& repam, const REAL msg, const INDEX msg_dim);

    // not used currently, as primal rounding does not give good results
    //template<typename LEFT_FACTOR, typename RIGHT_FACTOR>
    //void ComputeRightFromLeftPrimal(PrimalSolutionStorage::Element left, LEFT_FACTOR* l,
    //                                PrimalSolutionStorage::Element right, RIGHT_FACTOR* r);

    
  private:
    const INDEX numberOfLabels_,numberOfVarsLeft_,numberOfVarsRight_,SumBound_;
    INDEX upSize_,leftSize_,rightSize_,regSize_;
  };

  inline DiscreteTomographyMessageCountingPairwise::DiscreteTomographyMessageCountingPairwise(INDEX numberOfLabels,INDEX numberOfVarsLeft,INDEX numberOfVarsRight,INDEX SumBound)
    : numberOfLabels_(numberOfLabels),numberOfVarsLeft_(numberOfVarsLeft),numberOfVarsRight_(numberOfVarsRight),SumBound_(SumBound) {
    assert(numberOfLabels_ > 1);
    assert(numberOfVarsLeft_ > 0);
    assert(numberOfVarsRight_ > 0);
    assert(SumBound_ > 0);
        
    upSize_ = pow(numberOfLabels_,2)*std::min(((numberOfVarsLeft_+numberOfVarsRight_)*(numberOfLabels_-1)+1),SumBound);
    leftSize_ = pow(numberOfLabels_,2)*std::min((numberOfVarsLeft_*(numberOfLabels_-1)+1),SumBound);
    rightSize_ = pow(numberOfLabels_,2)*std::min((numberOfVarsRight_*(numberOfLabels_-1)+1),SumBound);
    
    regSize_ = pow(numberOfLabels_,2);
  }

  template<class LEFT_FACTOR_TYPE,class RIGHT_FACTOR_TYPE>
  void DiscreteTomographyMessageCountingPairwise::CreateConstraints(LpInterfaceAdapter* lp,LEFT_FACTOR_TYPE* LeftFactor,RIGHT_FACTOR_TYPE* RightFactor) const {
    for(INDEX i=0;i<regSize_;i++){
      LinExpr lhs = lp->CreateLinExpr() + 0.0;
      LinExpr rhs = lp->CreateLinExpr() + 0.0;
      assert(lp->GetLeftFactorSize() == regSize_);
      assert(lp->GetRightFactorSize() == upSize_+leftSize_+rightSize_+regSize_);
      bool least = false;
      if(lp->IsLeftObjective(i)){
        lhs += lp->GetLeftVariable(i);
        least = true;
      }
      if(lp->IsRightObjective(upSize_+leftSize_+rightSize_+i)){
        rhs += lp->GetRightVariable(upSize_+leftSize_+rightSize_+i);
        least = true;
      }
      if(least){
        lp->addLinearEquality(lhs,rhs);
      }
    }
  }
  
  template<typename LEFT_FACTOR, typename REPAM_ARRAY, typename MSG>
  void DiscreteTomographyMessageCountingPairwise::MakeLeftFactorUniform(LEFT_FACTOR* f_left, const REPAM_ARRAY& repam_left, MSG& msg, const REAL omega)
  {
     assert(msg.size() == pow(numberOfLabels_,2));
     assert(repam_left.size() == pow(numberOfLabels_,2));

     for(INDEX i=0;i<pow(numberOfLabels_,2);i++){
        msg[i] -= omega*repam_left[i];
     }   
  }

  template<typename RIGHT_FACTOR, typename REPAM_ARRAY, typename MSG>
  void DiscreteTomographyMessageCountingPairwise::MakeRightFactorUniform(RIGHT_FACTOR* f_right, const REPAM_ARRAY& repam_right, MSG& msg, const REAL omega)
  {
    assert(msg.size() == pow(numberOfLabels_,2));
    assert(repam_right.size() == ((*f_right).getSize(DiscreteTomographyFactorCounting::NODE::up) +
                                  (*f_right).getSize(DiscreteTomographyFactorCounting::NODE::left) +
                                  (*f_right).getSize(DiscreteTomographyFactorCounting::NODE::right) +
                                  (*f_right).getSize(DiscreteTomographyFactorCounting::NODE::reg)));

    std::vector<REAL> msg_v(pow(numberOfLabels_,2),std::numeric_limits<REAL>::infinity());
    if( DiscreteTomo::AlgorithmThreshold < (*f_right).getSize(DiscreteTomographyFactorCounting::NODE::up) ){
      DiscreteTomo::MessageCalculation_MinConv_Reg(f_right,repam_right, msg_v,numberOfLabels_);
    } else {
      DiscreteTomo::MessageCalculation_Naive_Reg(f_right,repam_right, msg_v,numberOfLabels_);
    }

    for(INDEX i=0;i<msg_v.size();i++){
      msg[i] -= omega*msg_v[i];
    }
    
  }

  template<typename G>
  void DiscreteTomographyMessageCountingPairwise::RepamLeft(G& repam, const REAL msg, const INDEX msg_dim){

    auto f = repam.GetFactor();
    assert( repam.size() == pow(numberOfLabels_,2));

    assert(msg_dim < pow(numberOfLabels_,2));
    assert(repam[msg_dim] > -std::numeric_limits<REAL>::max());
    if( std::isfinite(msg) ){ repam[msg_dim] += msg; } 
    
    else{ repam[msg_dim] = std::numeric_limits<REAL>::infinity(); }
    assert(repam[msg_dim] > -1.0e-02);
  }

  template<typename G>
  void DiscreteTomographyMessageCountingPairwise::RepamRight(G& repam, const REAL msg, const INDEX msg_dim){

    auto f = repam.GetFactor();
    assert( repam.size() == (f->getSize(DiscreteTomographyFactorCounting::NODE::up) +
                             f->getSize(DiscreteTomographyFactorCounting::NODE::left) +
                             f->getSize(DiscreteTomographyFactorCounting::NODE::right) +
                             f->getSize(DiscreteTomographyFactorCounting::NODE::reg))
            );
    
    assert(msg_dim < pow(numberOfLabels_,2));
    assert(repam[f->getSize(DiscreteTomographyFactorCounting::NODE::up) +
                 f->getSize(DiscreteTomographyFactorCounting::NODE::left) +
                 f->getSize(DiscreteTomographyFactorCounting::NODE::right) +
                 msg_dim] > -std::numeric_limits<REAL>::max());
    
    if( std::isfinite(msg) ){
      repam[f->getSize(DiscreteTomographyFactorCounting::NODE::up) +
            f->getSize(DiscreteTomographyFactorCounting::NODE::left) +
            f->getSize(DiscreteTomographyFactorCounting::NODE::right) +
            msg_dim] +=  msg;
      assert(repam[f->getSize(DiscreteTomographyFactorCounting::NODE::up) +
                   f->getSize(DiscreteTomographyFactorCounting::NODE::left) +
                   f->getSize(DiscreteTomographyFactorCounting::NODE::right) +
                   msg_dim] > -std::numeric_limits<REAL>::max());
    }
    else{
      repam[f->getSize(DiscreteTomographyFactorCounting::NODE::up) +
            f->getSize(DiscreteTomographyFactorCounting::NODE::left) +
            f->getSize(DiscreteTomographyFactorCounting::NODE::right) +
            msg_dim] = std::numeric_limits<REAL>::infinity();
    }
  }

  /*
  template<typename LEFT_FACTOR, typename RIGHT_FACTOR>
  void DiscreteTomographyMessageCountingPairwise::ComputeRightFromLeftPrimal(PrimalSolutionStorage::Element left, LEFT_FACTOR* leftFactor,
                                                                             PrimalSolutionStorage::Element right, RIGHT_FACTOR* rightFactor){
    assert(rightFactor->getSize(DiscreteTomographyFactorCounting::NODE::up) == upSize_);
    assert(rightFactor->getSize(DiscreteTomographyFactorCounting::NODE::left) == leftSize_);
    assert(rightFactor->getSize(DiscreteTomographyFactorCounting::NODE::right) == rightSize_);
    assert(rightFactor->getSize(DiscreteTomographyFactorCounting::NODE::reg) == regSize_);

    
    INDEX opt = 0;
    INDEX count = 0;
    INDEX noUnkwn = 0;

    // Check if there is one "true" label for the pairwise term
    for(INDEX i=0;i<pow(numberOfLabels_,2);i++){
      if( left[i] == true ){
        opt = i; count++;
        assert(right[upSize_ + leftSize_ + rightSize_ + i] != 0); // do not overwrite results!
      }
      if( left[i] == unknownState ){ opt = i; noUnkwn++; };
      right[upSize_ + leftSize_ + rightSize_ + i] = left[i]; // copy pairwise primal to ternary
    }
    assert(count <= 1);
    assert(count != 1 || noUnkwn == 0);
    assert(noUnkwn != 0 || count == 1);
    assert(noUnkwn != 0 || count != 0);
    
    // Calculate left and right label 
    INDEX b = opt % numberOfLabels_;
    INDEX c = ((opt - b)/numberOfLabels_) % numberOfLabels_;

    // Check for leafs
    if(noUnkwn == 1 && count == 0){ count = 1; }

    // If the left variable is a leaf
    INDEX lIdx = b + b*numberOfLabels_ + b*pow(numberOfLabels_,2);
    if(numberOfVarsLeft_ == 1 && count == 1 && lIdx < leftSize_){
      for(INDEX i=0;i<leftSize_;i++){
        right[upSize_ + i]=false;
      }
      
      right[upSize_ + lIdx]=true;
    }

    // If the right variable is a leaf
    INDEX rIdx = c + c*numberOfLabels_ + c*pow(numberOfLabels_,2);
    if(numberOfVarsRight_ == 1 && count == 1 && rIdx < rightSize_){
      for(INDEX i=0;i<rightSize_;i++){
        right[upSize_ + leftSize_ + i]=false;
      }

      right[upSize_ + leftSize_ + rIdx]=true;
    }
     
  }
  */

}


#endif // LP_MP_DT_COUNTING_MESSAGE_HXX
