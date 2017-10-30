#ifndef LP_MP_DT_COUNTING_PAIRWISE_MESSAGE_HXX
#define LP_MP_DT_COUNTING_PAIRWISE_MESSAGE_HXX

#include "LP_MP.h"
#include "tropical_convolution.hxx"
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










}


#endif // LP_MP_DT_COUNTING_MESSAGE_HXX
