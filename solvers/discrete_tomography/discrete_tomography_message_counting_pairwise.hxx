#ifndef LP_MP_DT_COUNTING_PAIRWISE_MESSAGE_HXX
#define LP_MP_DT_COUNTING_PAIRWISE_MESSAGE_HXX

#include "LP_MP.h"
#include "minConv.hxx"
#include <math.h> 
#include "static_if.hxx"

namespace LP_MP {
    
//using MinConv = MinConv<REAL,INDEX>;

// regularizer to counting factor, when counting factor has intermediate sum=0, i.e. no intermediate node
enum class CountingPairwiseMessageType { left, center, right };
template<CountingPairwiseMessageType TYPE>
class DiscreteTomographyMessageCountingPairwise2{
public:
   DiscreteTomographyMessageCountingPairwise2(const bool transpose) : transpose_(transpose) {
      assert(transpose == false);
   }

   template<typename LEFT_FACTOR, typename MSG>
   void send_message_to_right(const LEFT_FACTOR& f_left, MSG& msg, const REAL omega){
      //std::cout << "\nbefore:\n";
      for(INDEX x1=0; x1<f_left.dim1(); ++x1) {
         for(INDEX x2=0; x2<f_left.dim2(); ++x2) {
            //std::cout << f_left(x1,x2) << ", ";
            assert(std::isnan(f_left(x1,x2)) == false);
         }
      }
      //std::cout << "\nafter, omega=" << omega << ":\n";
      msg -= omega*f_left;
      for(INDEX x1=0; x1<f_left.dim1(); ++x1) {
         for(INDEX x2=0; x2<f_left.dim2(); ++x2) {
            //std::cout << f_left(x1,x2) << ", ";
            assert(std::isnan(f_left(x1,x2)) == false);
         }
      }
      //std::cout << "\n--------------------------\n";
   }

   template<typename RIGHT_FACTOR, typename MSG>
   void send_message_to_left(const RIGHT_FACTOR& f_right, MSG& msg, const REAL omega){
      if(TYPE == CountingPairwiseMessageType::center) {
         matrix<REAL> msg_tmp(f_right.no_center_left_labels(), f_right.no_center_right_labels());
         f_right.MessageCalculation_Reg(msg_tmp);
         msg -= omega*msg_tmp;
      } else if(TYPE == CountingPairwiseMessageType::left) {
         matrix<REAL> msg_tmp(f_right.no_left_labels(), f_right.no_center_left_labels());
         f_right.MessageCalculation_Naive_Left(msg_tmp);
         msg -= omega*msg_tmp;
      } else if(TYPE == CountingPairwiseMessageType::right) {
         matrix<REAL> msg_tmp(f_right.no_center_right_labels(), f_right.no_right_labels());
         f_right.MessageCalculation_Naive_Right(msg_tmp);
         msg -= omega*msg_tmp;
      } else {
         assert(false);
      }
   }

   template<typename LEFT_FACTOR, typename MSG>
   void RepamLeft(LEFT_FACTOR& f_left, const MSG& msg)
   {
      for(INDEX x1=0; x1<f_left.dim1(); ++x1) {
         for(INDEX x2=0; x2<f_left.dim2(); ++x2) {
            assert(std::isnan(msg(x1,x2)) == false);
         }
      }

      for(INDEX x1=0; x1<f_left.dim1(); ++x1) {
         for(INDEX x2=0; x2<f_left.dim2(); ++x2) {
            f_left.cost(x1,x2) += normalize( msg(x1,x2) );
            //static_if<TYPE == CountingPairwiseMessageType::center>([&](auto f){
            //});
            //static_if<TYPE == CountingPairwiseMessageType::left>([&](auto f){
            //      f_left(x1,x2) += f(msg)(x1,x2);
            //});
            //static_if<TYPE == CountingPairwiseMessageType::right>([&](auto f){
            //      f_left(x1,x2) += f(msg)(x1,x2);
            //});
         }
      }
      //std::cout << "\n";
      //assert(f_left.LowerBound() < std::numeric_limits<REAL>::infinity());
   }

   template<typename RIGHT_FACTOR, typename MSG>
   void RepamRight(RIGHT_FACTOR& f_right, const MSG& msg)
   {
      if(TYPE == CountingPairwiseMessageType::center) {
         auto& reg = f_right.reg();
         if(!transpose_) {
            for(INDEX x1=0; x1<f_right.no_center_left_labels(); ++x1) {
               for(INDEX x2=0; x2<f_right.no_center_right_labels(); ++x2) {
                  reg(x1,x2) += normalize( msg(x1,x2) );
               }
            }
         } else {
            assert(false);
         } 
      }

      if(TYPE == CountingPairwiseMessageType::left) {
         auto& left = f_right.left();
         assert(left.dim3() == 1);
         if(!transpose_) {
            for(INDEX x1=0; x1<f_right.no_left_labels(); ++x1) {
               for(INDEX x2=0; x2<f_right.no_center_left_labels(); ++x2) {
                  left(x1,x2,0) += normalize( msg(x1,x2) );
               }
            }
         } else {
            assert(false);
         } 
      }

      if(TYPE == CountingPairwiseMessageType::right) {
         auto& right = f_right.right();
         assert(right.dim3() == 1);
         if(!transpose_) {
            for(INDEX x1=0; x1<f_right.no_center_right_labels(); ++x1) {
               for(INDEX x2=0; x2<f_right.no_right_labels(); ++x2) {
                  right(x1,x2,0) += normalize( msg(x1,x2) );
               }
            }
         } else {
            assert(false);
         } 
      }
   }
   private:
   const bool transpose_; // pairwise potential in discrete tomography factor might be transposed
};

// regularizer to counting factor, when left or right counting factor has intermediate sum=0,...,no_labels i.e. exactly one intermediate node
// the first template parameter specifies whether message applies to left or to right counting subfactor in DiscreteCountiangFactor.
// The second template parameter specifies whether the message applies to left and middle or middle and right labels in counting factor
template<Chirality COUNTING_FACTOR, Chirality DIRECTION>
class DiscreteTomographyMessageCountingPairwise3{
public:
   DiscreteTomographyMessageCountingPairwise3(const bool transpose) : transpose_(transpose) {
      assert(transpose == false);
   }
   ~DiscreteTomographyMessageCountingPairwise3() {
      static_assert(COUNTING_FACTOR == Chirality::left || COUNTING_FACTOR == Chirality::right,"");
      static_assert(DIRECTION == Chirality::left || DIRECTION == Chirality::right,"");
   }

   template<typename LEFT_FACTOR, typename MSG>
   void send_message_to_right(const LEFT_FACTOR& f_left, MSG& msg, const REAL omega){
      msg -= omega*f_left;
      //for(INDEX x1=0;x1<f_left.dim1();x1++){
      //   for(INDEX x2=0;x2<f_left.dim2();x2++){
      //      msg[i] -= omega*repam_left[i];
      //   }   
      //}
   }


   template<typename RIGHT_FACTOR, typename MSG>
   void send_message_to_left(const RIGHT_FACTOR& f_right, MSG& msg, const REAL omega){
      std::array<INDEX,3> dim;
      if(COUNTING_FACTOR == Chirality::left) {
         dim = {f_right.no_left_labels(), f_right.no_center_left_labels(), f_right.left_sum_size()};
      } else {
         dim = {f_right.no_center_right_labels(), f_right.no_right_labels(), f_right.right_sum_size()};
      }

      tensor3<REAL> msg_tmp(dim[0], dim[1], dim[2]);
      
      if(COUNTING_FACTOR == Chirality::left) {
         assert(f_right.no_left_labels() ==  f_right.left_sum_size());
         assert(f_right.no_center_left_labels() ==  f_right.left_sum_size());
         f_right.MessageCalculation_Naive_Left(msg_tmp);
      } else {
         assert(f_right.no_right_labels() ==  f_right.right_sum_size());
         assert(f_right.no_center_right_labels() ==  f_right.right_sum_size());
         f_right.MessageCalculation_Naive_Right(msg_tmp);
      }

      if(COUNTING_FACTOR == Chirality::left) {
         if(DIRECTION == Chirality::left) {
            dim = {f_right.no_left_labels(), f_right.left_sum_size(), 0};
         } else {
            dim = {f_right.left_sum_size(), f_right.no_center_left_labels(), 0}; 
         }
      } else {
         if(DIRECTION == Chirality::left) {
            dim = {f_right.no_center_right_labels(), f_right.right_sum_size(), 0}; 
         } else {
            dim = {f_right.right_sum_size(), f_right.no_right_labels(), 0};
         } 
      }

      matrix<REAL> msg_marg(dim[0], dim[1], std::numeric_limits<REAL>::infinity());

      if(!transpose_) {
         for(INDEX x1=0; x1<msg_tmp.dim1(); ++x1) {
            for(INDEX x2=0; x2<msg_tmp.dim2(); ++x2) {
               for(INDEX sum=0;sum<msg_tmp.dim3(); ++sum) {
                  if(DIRECTION == Chirality::left) {
                     msg_marg(x1,sum) = std::min(msg_marg(x1,sum), msg_tmp(x1,x2,sum));
                  } else {
                     msg_marg(sum,x2) = std::min(msg_marg(sum,x2), msg_tmp(x1,x2,sum));
                  }
               }
            }
         }
      } else {
         assert(false);
      }
      msg -= omega*msg_marg;
   }

   template<typename LEFT_FACTOR, typename MSG>
   void RepamLeft(LEFT_FACTOR& f_left, MSG& msg)
   {
      for(INDEX x1=0; x1<f_left.dim1(); ++x1) {
         for(INDEX x2=0; x2<f_left.dim2(); ++x2) {
            f_left.cost(x1,x2) += normalize( msg(x1,x2) );
         }
      }
   }

   template<typename RIGHT_FACTOR, typename MSG>
   void RepamRight(RIGHT_FACTOR& f_right, MSG& msg)
   {
      auto& f = (COUNTING_FACTOR == Chirality::left) ? f_right.left() : f_right.right();
      //assert(f.dim3() == f_right.no_labels());

      //if(COUNTING_FACTOR == Chirality::left) {
      //   assert(f_right.left_sum_size() == f_right.no_labels());
      //} else {
      //   assert(f_right.right_sum_size() == f_right.no_labels());
      //}

      if(!transpose_) {
         if(COUNTING_FACTOR == Chirality::left) {
            for(INDEX x1=0; x1<f_right.no_left_labels(); ++x1) {
               for(INDEX x2=0; x2<f_right.no_center_left_labels(); ++x2) {
                  for(INDEX sum=0;sum<f_right.left_sum_size(); ++sum) {
                     if(DIRECTION == Chirality::left) {
                        f(x1,x2,sum) += normalize( msg(x1,sum) );
                     } else {
                        f(x1,x2,sum) += normalize( msg(sum,x2) );
                     }
                  }
               }
            }
         } else {
            for(INDEX x1=0; x1<f_right.no_center_right_labels(); ++x1) {
               for(INDEX x2=0; x2<f_right.no_right_labels(); ++x2) {
                  for(INDEX sum=0;sum<f_right.right_sum_size(); ++sum) {
                     if(DIRECTION == Chirality::left) {
                        f(x1,x2,sum) += normalize( msg(x1,sum) );
                     } else {
                        f(x1,x2,sum) += normalize( msg(sum,x2) );
                     }
                  }
               }

            }
         }
      } else {
         assert(false);
      } 
   }

private:
   bool transpose_;
};












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
      LinExpr lhs = lp->CreateLinExpr();
      LinExpr rhs = lp->CreateLinExpr();
      assert(lp->GetLeftFactorSize() == regSize_);
      assert(lp->GetRightFactorSize() == upSize_+leftSize_+rightSize_+regSize_);
      lhs += lp->GetLeftVariable(i);
      rhs += lp->GetRightVariable(upSize_+leftSize_+rightSize_+i);
      lp->addLinearEquality(lhs,rhs);
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
