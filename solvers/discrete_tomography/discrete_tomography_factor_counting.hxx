#ifndef LPMP_DTOMOGRAPHY_FACTOR_COUNTING_HXX
#define LPMP_DTOMOGRAPHY_FACTOR_COUNTING_HXX

#include "LP_MP.h"
#include "tropical_convolution.hxx"
#include <math.h> 

namespace LP_MP{

  constexpr static INDEX MinSumConvolutionThreshold = 1000000;

  class DiscreteTomographyFactorCounting2{
  public:

     DiscreteTomographyFactorCounting2(const INDEX no_left_labels, const INDEX left_sum_size, const INDEX no_center_left_labels, const INDEX no_center_right_labels, const INDEX right_sum_size, const INDEX no_right_labels, const INDEX up_sum_size)
        : reg_(no_center_left_labels, no_center_right_labels, 0.0),
        up_(no_left_labels, no_right_labels, up_sum_size,0.0), 
        left_(no_left_labels, no_center_left_labels, left_sum_size,0.0), 
        right_(no_center_right_labels, no_right_labels, right_sum_size,0.0)
     {
        assert(no_left_labels > 0);
        assert(no_center_left_labels > 0);
        assert(no_center_right_labels > 0);
        assert(no_right_labels > 0);
        assert(left_sum_size > 0);
        assert(right_sum_size > 0);
        assert(up_sum_size > 0);
        //assert(up_sum_size <= left_sum_size + no_center_left_labels + no_center_right_labels + right_sum_size);
        assert(no_center_left_labels-1 + no_center_right_labels-1 + left_sum_size-1 + right_sum_size-1 >= up_sum_size-1);
        //assert(min_conv_lower_bound() < 0.000000001);
     }
     DiscreteTomographyFactorCounting2(const INDEX no_labels, const INDEX left_sum_size, const INDEX right_sum_size, const INDEX up_sum_size)
        : DiscreteTomographyFactorCounting2(no_labels, left_sum_size, no_labels, no_labels, right_sum_size, no_labels, up_sum_size)
     {}

     REAL eval(const INDEX x_left, const INDEX left_sum, const INDEX x_center_left, const INDEX x_center_right, const INDEX right_sum, const INDEX x_right) const {
        assert(x_left < no_left_labels());
        assert(x_center_left < no_center_left_labels());
        assert(x_center_right < no_center_right_labels());
        assert(x_right < no_right_labels());
        assert(left_sum + x_center_left + x_center_right + right_sum < up_sum_size());
        return up_(x_left, x_right, x_center_left + x_center_right + left_sum + right_sum) 
           + left_(x_left, x_center_left, left_sum) 
           + right_(x_center_right, x_right, right_sum) 
           + reg_(x_center_left, x_center_right);
     }

     template<typename ARRAY>
     void summation_cost(const ARRAY& cost) {
        assert(cost.size() <= up_sum_size());
        for(INDEX x1=0; x1<up_.dim1(); ++x1) {
           for(INDEX x2=0; x2<up_.dim2(); ++x2) {
              for(INDEX i=0; i<up_.dim3(); ++i) {
                 if(x1+x2+i < cost.size()) {
                    up_(x1,x2,i) = cost[x1+x2+i];
                 } else {
                    up_(x1,x2,i) = std::numeric_limits<REAL>::infinity();
                 }
              }
           } 
        }
     }

     // explicitly iterate over all possibilities
     REAL naive_lower_bound() const {
        REAL bound = std::numeric_limits<REAL>::infinity();

        for_each_label_sum([&](const INDEX x_l, const INDEX lz, const INDEX x_cl, const INDEX x_cr, const INDEX rz, const INDEX x_r) {
              bound = std::min(bound, eval(x_l, lz, x_cl, x_cr, rz, x_r) );
              });
        return bound;
     }

     // use min sum convolution to replace inner loops for cardinality variables
     REAL min_conv_lower_bound() const {
        REAL bound = std::numeric_limits<REAL>::infinity();

        for_each_label([&](const INDEX x_l, const INDEX x_cl, const INDEX x_cr, const INDEX x_r) {
              assert(up_sum_size()-1 >= x_cl + x_cr);
              // set max_sum_size as sum_max as below! Possibly make do not give left_sum_size and right_sum_size to mc but min(left_sum_size(), max_sum_size()) etc.
              const INDEX max_sum_size = std::min(up_sum_size(), left_sum_size() + right_sum_size()-1 );
              assert(max_sum_size >= left_sum_size() && max_sum_size >= right_sum_size());
              auto op = [=](INDEX i,INDEX j){ return std::min(i+j, max_sum_size); };
              auto z_left = [&](INDEX k){ return left_(x_l,x_cl, k); };
              auto z_right = [&](INDEX k){ return right_(x_cr, x_r, k); }; 

              vector<REAL> mc(max_sum_size);
              tropical_concolution::min_conv(z_left.begin(), z_left.end(), z_right.begin(), z_right.end(), mc.begin(), mc.end());

              const INDEX sum_max = std::min(max_sum_size, up_sum_size() - x_cl - x_cr);
              for(INDEX sum=0; sum<sum_max; sum++){
               assert(sum == op(mc.getIdxA(sum),mc.getIdxB(sum)));
               assert(sum == (mc.getIdxA(sum) + mc.getIdxB(sum)));

               const REAL val = mc.getConv(sum) + up_(x_l, x_r, sum + x_cl + x_cr) + reg_(x_cl, x_cr);
               bound = std::min(bound, val);
              }
              });
        return bound;
    }
    
     /*
        for_each_label_min_conv([&](const INDEX x_l, const INDEX x_cl, const INDEX x_cr, const INDEX x_r, const MinConv& mc) {

              for(INDEX sum=0;sum<up_sum_size()-x_cl-x_cr;sum++){
               const REAL val = mc.getConv(sum+x_cl+x_cr) + up_(x_l, x_r, sum) + reg_(x_cl, x_cr);
               bound = std::min(bound, val);
              }
              });
        return bound; 
     }
*/
     template<typename LAMBDA>
     void for_each_label_sum(LAMBDA f) const {
        for(INDEX x_l=0; x_l<no_left_labels(); ++x_l) {
           for(INDEX x_cl=0; x_cl<no_center_left_labels(); ++x_cl) {
              for(INDEX x_cr=0; x_cr<no_center_right_labels(); ++x_cr) {
                 for(INDEX x_r=0; x_r<no_right_labels(); ++x_r) {
                    for(INDEX lz=0;lz<left_sum_size();lz++){
                       for(INDEX rz=0;rz<right_sum_size() && rz+lz+x_cl+x_cr<up_sum_size();rz++){
                          f(x_l, lz, x_cl, x_cr, rz, x_r);
                       }
                    }
                 }
              }
           }
        }
     }

     template<typename LAMBDA>
     void for_each_label(LAMBDA f) const {
        for(INDEX x_l=0; x_l<no_left_labels(); ++x_l) {
           for(INDEX x_cl=0; x_cl<no_center_left_labels(); ++x_cl) {
              for(INDEX x_cr=0; x_cr<no_center_right_labels(); ++x_cr) {
                 for(INDEX x_r=0; x_r<no_right_labels(); ++x_r) {
                    f(x_l,x_cl, x_cr, x_r);
                 }
              }
           }
        }
     }

     template<typename LAMBDA>
     void for_each_label_min_conv(LAMBDA f) const
     {
        for(INDEX x_l=0; x_l<no_left_labels(); ++x_l) {
           for(INDEX x_cl=0; x_cl<no_center_left_labels(); ++x_cl) {
              for(INDEX x_cr=0; x_cr<no_center_right_labels(); ++x_cr) {
                 for(INDEX x_r=0; x_r<no_right_labels(); ++x_r) {
                    assert(up_sum_size() >= x_cl + x_cr);
                    auto op = [&](INDEX i,INDEX j){ return std::min(i+j, up_sum_size() - x_cl - x_cr); };
                    auto z_left = [&](INDEX k){ return left_(x_l,x_cl, k); };
                    auto z_right = [&](INDEX k){ return right_(x_cr, x_r, k); }; 

                    MinConv mc(z_left,z_right,left_sum_size(),right_sum_size(),up_sum_size() - x_cl - x_cr);
                    mc.CalcConv(op,z_left,z_right);

                    f(x_l,x_cl, x_cr, x_r, mc);
                 }
              }
           }
        }
     }

     REAL LowerBound() const {
        if(up_.dim3() < MinSumConvolutionThreshold) { // check constant!
           return naive_lower_bound();
        } else {
           return min_conv_lower_bound();
        }
     }

     REAL EvaluatePrimal() const {
        return std::numeric_limits<REAL>::infinity();
     }

     INDEX size() const { return up_.size() + left_.size() + right_.size() + reg_.size(); }
     INDEX no_left_labels() const { return left_.dim1(); }
     INDEX no_center_left_labels() const { return left_.dim2(); }
     INDEX no_center_right_labels() const { return right_.dim1(); }
     INDEX no_right_labels() const { return right_.dim2(); }
     INDEX up_sum_size() const { return up_.dim3(); }
     INDEX left_sum_size() const { return left_.dim3(); }
     INDEX right_sum_size() const { return right_.dim3(); }
     REAL& reg(const INDEX x_cl, const INDEX x_cr) { return reg_(x_cl, x_cr); }
     REAL& right(const INDEX x_cr, const INDEX x_r, const INDEX sum) { return right_(x_cr, x_r, sum); }
     REAL& left(const INDEX x_l, const INDEX x_cl, const INDEX sum) { return left_(x_l, x_cl, sum); }
     REAL& up(const INDEX x_left, const INDEX x_right, const INDEX sum) { return up_(x_left, x_right, sum); }
     matrix<REAL>& reg() { return reg_; }
     tensor3<REAL>& right() { return right_; }
     tensor3<REAL>& left() { return left_; }
     tensor3<REAL>& up() { return up_; }



     // marginalization operations
     template<typename MSG>
     void MessageCalculation_Up(MSG& msg) const {
        if(up_sum_size() > MinSumConvolutionThreshold) {
           MessageCalculation_MinConv_Up(msg);
        } else {
           MessageCalculation_Naive_Up(msg);
        } 
     }
     template<typename MSG>
     void MessageCalculation_Left(MSG& msg) const {
        if(up_sum_size() > MinSumConvolutionThreshold) {
           MessageCalculation_MinConv_Left(msg);
        } else {
           MessageCalculation_Naive_Left(msg);
        } 
     }
     template<typename MSG>
     void MessageCalculation_Right(MSG& msg) const {
        if(up_sum_size() > MinSumConvolutionThreshold) {
           MessageCalculation_MinConv_Right(msg);
        } else {
           MessageCalculation_Naive_Right(msg);
        } 
     }
     template<typename MSG>
     void MessageCalculation_Reg(MSG& msg) const {
        if(up_sum_size() > MinSumConvolutionThreshold) {
           MessageCalculation_MinConv_Reg(msg);
        } else {
           MessageCalculation_Naive_Reg(msg);
        } 
     }

    template<typename MSG>
    void MessageCalculation_Naive_Up(MSG& msg) const {
       std::fill(msg.begin(), msg.end(), std::numeric_limits<REAL>::infinity());

        for_each_label_sum([&](const INDEX x_l, const INDEX lz, const INDEX x_cl, const INDEX x_cr, const INDEX rz, const INDEX x_r) {
           const REAL value = eval(x_l, lz, x_cl, x_cr, rz, x_r);
           const INDEX sum = lz + x_cl + x_cr + rz;
           msg(x_l, x_r, sum) = std::min(msg(x_l, x_r, sum), value);
           });
    }

    template<typename MSG>
    void MessageCalculation_Naive_Left(MSG& msg) const {
        std::fill(msg.begin(), msg.end(), std::numeric_limits<REAL>::infinity());

        for_each_label_sum([&](const INDEX x_l, const INDEX lz, const INDEX x_cl, const INDEX x_cr, const INDEX rz, const INDEX x_r) {
           const REAL value = eval(x_l, lz, x_cl, x_cr, rz, x_r);
           msg(x_l, x_cl, lz) = std::min(msg(x_l, x_cl, lz), value);
           });
    }

    template<typename MSG>
    void MessageCalculation_Naive_Right(MSG& msg) const {
        std::fill(msg.begin(), msg.end(), std::numeric_limits<REAL>::infinity());

        for_each_label_sum([&](const INDEX x_l, const INDEX lz, const INDEX x_cl, const INDEX x_cr, const INDEX rz, const INDEX x_r) {
           const REAL value = eval(x_l, lz, x_cl, x_cr, rz, x_r);
           msg(x_cr, x_r, rz) = std::min(msg(x_cr, x_r, rz), value);
           });
    }


    template<typename MSG>
    void MessageCalculation_Naive_Reg(MSG& msg) const {
        std::fill(msg.begin(), msg.end(), std::numeric_limits<REAL>::infinity());

        for_each_label_sum([&](const INDEX x_l, const INDEX lz, const INDEX x_cl, const INDEX x_cr, const INDEX rz, const INDEX x_r) {
           const REAL value = eval(x_l, lz, x_cl, x_cr, rz, x_r);
           msg(x_cl, x_cr) = std::min(msg(x_cl, x_cr), value);
           });
    }
    
    template<typename MSG>
    void MessageCalculation_MinConv_Up(MSG& msg) const {
        std::fill(msg.begin(), msg.end(), std::numeric_limits<REAL>::infinity());

        for_each_label([&](const INDEX x_l, const INDEX x_cl, const INDEX x_cr, const INDEX x_r) {
              assert(up_sum_size()-1 >= x_cl + x_cr);
              // set max_sum_size as sum_max as below! Possibly make do not give left_sum_size and right_sum_size to mc but min(left_sum_size(), max_sum_size()) etc.
              const INDEX max_sum_size = std::min(up_sum_size(), left_sum_size() + right_sum_size()-1 );
              assert(max_sum_size >= left_sum_size() && max_sum_size >= right_sum_size());
              auto op = [=](INDEX i,INDEX j){ return i+j; };
              auto z_left = [&](INDEX k){ return left_(x_l,x_cl, k); };
              auto z_right = [&](INDEX k){ return right_(x_cr, x_r, k); }; 

              MinConv mc(z_left, z_right, left_sum_size(), right_sum_size(), max_sum_size);
              mc.CalcConv(op,z_left,z_right);

              const INDEX sum_max = std::min(max_sum_size, up_sum_size() - x_cl - x_cr);
              for(INDEX sum=0;  sum<sum_max; sum++) {
               assert(sum == op(mc.getIdxA(sum),mc.getIdxB(sum)));
               assert(sum == (mc.getIdxA(sum) + mc.getIdxB(sum)));

               const REAL val = mc.getConv(sum) + up_(x_l, x_r, sum + x_cl + x_cr) + reg_(x_cl, x_cr);
               msg(x_l, x_r, sum+x_cl+x_cr) = std::min(msg(x_l, x_r, sum+x_cl+x_cr), val);
               }
               });
    }
    
    template<typename MSG>
    void MessageCalculation_MinConv_Left(MSG& msg) const {
        std::fill(msg.begin(), msg.end(), std::numeric_limits<REAL>::infinity());

        for_each_label([&](const INDEX x_l, const INDEX x_cl, const INDEX x_cr, const INDEX x_r) {
              assert(up_sum_size() >= x_cl + x_cr);
              assert(up_sum_size() >= x_cl + x_cr);
              const INDEX left_size = up_sum_size();// - x_cl - x_cr;
              // total_sum = left_sum + x_cl + x_cr + right_sum => left_sum = total_sum - x_cl - x_cr - right_sum
              //const INDEX min_conv_size = left_sum_size();//std::min(left_size, left_sum_size());
              const INDEX min_conv_size = left_sum_size();//std::min(left_sum_size(), up_sum_size() - right_sum_size());
              //auto op = [=](INDEX i, INDEX j) { return i-x_cl-x_cr - j; }; 
              auto op = [&](INDEX i,INDEX j){ // 0 <= i-j < left_size
              if( i  < j + x_cl + x_cr ){ // i-j < 0
              return min_conv_size;
              }
              else{
              return (i - x_cl - x_cr - j < min_conv_size) ? i-x_cl-x_cr-j : min_conv_size;
              }
              };
              
              const INDEX shift = 0; //by which number to shift the result.
              auto z_left = [x_l,x_r,x_cl,x_cr,this](INDEX k){ return up_(x_l, x_r, k ); };
              auto z_right = [x_cr,x_r,this](INDEX k){ return right_(x_cr, x_r, k); }; 

              MinConv mc(z_left,z_right, left_size, right_sum_size(), min_conv_size);
              mc.CalcConv(op,z_left,z_right);

              for(INDEX left_sum=0; left_sum<min_conv_size; left_sum++) {
                 //REAL test_val = std::numeric_limits<REAL>::infinity();
                 //for(INDEX right_sum=0; right_sum<right_sum_size(); ++right_sum) {
                 //   if(left_sum + right_sum + x_cl + x_cr < up_sum_size()) {
                 //      test_val = std::min(test_val, up_(x_l, x_r, left_sum + right_sum + x_cl + x_cr) + right_(x_cr, x_r, right_sum));
                 //   }
                 //}
                 //const REAL val = test_val + left_(x_l, x_cl, left_sum) + reg_(x_cl, x_cr);
                 //assert(test_val == mc.getConv(left_sum));
                 const REAL val = mc.getConv(left_sum) + left_(x_l, x_cl, left_sum) + reg_(x_cl, x_cr);
                 msg(x_l, x_cl, left_sum) = std::min(msg(x_l, x_cl, left_sum), val);
              }
        });
    }

    template<typename MSG>
    void MessageCalculation_MinConv_Right(MSG& msg) const {
        std::fill(msg.begin(), msg.end(), std::numeric_limits<REAL>::infinity());

        for_each_label([&](const INDEX x_l, const INDEX x_cl, const INDEX x_cr, const INDEX x_r) {
              assert(up_sum_size() >= x_cl + x_cr);
              assert(up_sum_size() >= x_cl + x_cr);
              const INDEX left_size = up_sum_size();// - x_cl - x_cr;
              // total_sum = left_sum + x_cl + x_cr + right_sum => right_sum = total_sum - x_cl - x_cr - left_sum
              //const INDEX min_conv_size = left_sum_size();//std::min(left_size, left_sum_size());
              const INDEX min_conv_size = right_sum_size();//std::min(left_sum_size(), up_sum_size() - right_sum_size());
              //auto op = [=](INDEX i, INDEX j) { return i-x_cl-x_cr - j; }; 
              auto op = [&](INDEX i,INDEX j){ // 0 <= i-j < left_size
              if( i  < j + x_cl + x_cr ){ // i-j < 0
              return min_conv_size;
              }
              else{
              return (i - x_cl - x_cr - j < min_conv_size) ? i-x_cl-x_cr-j : min_conv_size;
              }
              };
              
              const INDEX shift = 0; //by which number to shift the result.
              auto z_left = [x_l,x_r,x_cl,x_cr,this](INDEX k){ return up_(x_l, x_r, k ); };
              auto z_right = [x_l,x_cl,this](INDEX k){ return left_(x_l, x_cl, k); }; 

              MinConv mc(z_left,z_right, left_size, left_sum_size(), min_conv_size);
              mc.CalcConv(op,z_left,z_right);

              for(INDEX right_sum=0; right_sum<min_conv_size; right_sum++) {
                 //REAL test_val = std::numeric_limits<REAL>::infinity();
                 //for(INDEX right_sum=0; right_sum<right_sum_size(); ++right_sum) {
                 //   if(left_sum + right_sum + x_cl + x_cr < up_sum_size()) {
                 //      test_val = std::min(test_val, up_(x_l, x_r, left_sum + right_sum + x_cl + x_cr) + right_(x_cr, x_r, right_sum));
                 //   }
                 //}
                 //const REAL val = test_val + left_(x_l, x_cl, left_sum) + reg_(x_cl, x_cr);
                 //assert(test_val == mc.getConv(left_sum));
                 const REAL val = mc.getConv(right_sum) + right_(x_cr, x_r, right_sum) + reg_(x_cl, x_cr);
                 msg(x_cr, x_r, right_sum) = std::min(msg(x_cr, x_r, right_sum), val);
              }
        });


        /*
        for_each_label([&](const INDEX x_l, const INDEX x_cl, const INDEX x_cr, const INDEX x_r) {
              assert(up_sum_size() >= x_cl + x_cr);
              auto op = [&](INDEX i,INDEX j){ // 0 <= i-j <= right_size
              if( i < j ){ 
              return right_sum_size();
              }
              else{
              return (i-j < right_sum_size()) ? i-j : right_sum_size();
              }
              };
              auto z_left = [&](INDEX k){ return up_(x_l,x_cl, k); };
              auto z_right = [&](INDEX k){ return left_(x_cr, x_r, k); }; 

              MinConv mc(z_left,z_right,up_sum_size(),left_sum_size(),right_sum_size());
              mc.CalcConv(op,z_left,z_right);

              for(INDEX right_sum=0;right_sum<right_sum_size();right_sum++){
               const REAL val = mc.getConv(right_sum+x_cl+x_cr) + right_(x_l, x_cl, right_sum) + reg_(x_cl, x_cr);
               msg(x_cr, x_r, right_sum) = std::min(msg(x_l, right_sum, x_cl), val);
               }
               });
               */
    }

    template<typename MSG>
    void MessageCalculation_MinConv_Reg(MSG& msg) const {
       assert(msg.dim1() == no_center_left_labels());
       assert(msg.dim2() == no_center_right_labels());
       std::fill(msg.begin(), msg.end(), std::numeric_limits<REAL>::infinity());

        for_each_label([&](const INDEX x_l, const INDEX x_cl, const INDEX x_cr, const INDEX x_r) {
              assert(up_sum_size()-1 >= x_cl + x_cr);
              // set max_sum_size as sum_max as below! Possibly make do not give left_sum_size and right_sum_size to mc but min(left_sum_size(), max_sum_size()) etc.
              const INDEX max_sum_size = std::min(up_sum_size(), left_sum_size() + right_sum_size()-1 );
              assert(max_sum_size >= left_sum_size() && max_sum_size >= right_sum_size());
              auto op = [=](INDEX i,INDEX j){ return i+j; };
              auto z_left = [&](INDEX k){ return left_(x_l,x_cl, k); };
              auto z_right = [&](INDEX k){ return right_(x_cr, x_r, k); }; 

              MinConv mc(z_left, z_right, left_sum_size(), right_sum_size(), max_sum_size);
              mc.CalcConv(op,z_left,z_right);

              const INDEX sum_max = std::min(max_sum_size, up_sum_size() - x_cl - x_cr);
              for(INDEX sum=0; sum<sum_max; sum++){
               assert(sum == op(mc.getIdxA(sum),mc.getIdxB(sum)));
               assert(sum == (mc.getIdxA(sum) + mc.getIdxB(sum)));

               const REAL val = mc.getConv(sum) + up_(x_l, x_r, sum + x_cl + x_cr) + reg_(x_cl, x_cr);
               msg(x_cl, x_cr) = std::min(msg(x_l, x_r), val);
              }
              });
    }

   void init_primal() { 
      primal_.left_sum = left_sum_size();
      primal_.right_sum = right_sum_size();
      primal_.left_label = no_left_labels();
      primal_.center_left_label = no_center_left_labels();
      primal_.center_right_label = no_center_right_labels();
      primal_.right_label = no_right_labels();
   }
   template<class ARCHIVE> void serialize_primal(ARCHIVE& ar) { ar( primal_.left_sum, primal_.right_sum, primal_.left_label, primal_.center_left_label, primal_.center_right_label, primal_.right_label ); }
   template<class ARCHIVE> void serialize_dual(ARCHIVE& ar) { ar( reg_, up_, left_, right_ ); }
  private:
    matrix<REAL> reg_;
    tensor3<REAL> up_, left_, right_; 

    struct { INDEX left_sum, right_sum, left_label, center_left_label, center_right_label, right_label; } primal_;
  };
  
  class DiscreteTomographyFactorCounting{

  public:

    enum class NODE {left,right,up,reg};
    
    DiscreteTomographyFactorCounting(INDEX numberOfLabels,INDEX numberOfVarsLeft,INDEX numberOfVarsRight,INDEX SumBound);
    
    //template<typename REPAM_ARRAY>
    //INDEX ComputeOptimalLabeling(const REPAM_ARRAY& repam) const { }; //--required

    template<typename REPAM_ARRAY>
    REAL LowerBound(const REPAM_ARRAY& repam) const; //--required

    //template<typename REPAM_ARRAY>
    //void MaximizePotentialAndComputePrimal(const REPAM_ARRAY& repam, PrimalSolutionStorage::Element primal) const;
    
    //template<typename REPAM_ARRAY>
    //void NarrowPrimal(const REPAM_ARRAY& repam, const REAL epsilon, typename PrimalSolutionStorage::Element primal);
    
    void PropagatePrimal(PrimalSolutionStorage::Element primal) const;
    
    template<typename REPAM_ARRAY>
    REAL EvaluatePrimal(const REPAM_ARRAY& repam, const PrimalSolutionStorage::Element primal) const; //--required

    INDEX size() const;
    
    //void WritePrimal(const PrimalSolutionStorage::Element, std::ofstream& fs) const; //

    // Lp Interface related
    INDEX GetNumberOfAuxVariables() const { return leftSize_*rightSize_; } 
    void CreateConstraints(LpInterfaceAdapter* lp) const;

    template<typename REPAM_ARRAY>
    void ReduceLp(LpInterfaceAdapter* lp,const REPAM_ARRAY& repam) const;
    
    INDEX getSize(NODE) const;

    template<typename REPAM_ARRAY>
    REAL eval(INDEX,INDEX,INDEX,const REPAM_ARRAY& repam);
    
  private:

    const INDEX numberOfLabels_,numberOfVarsLeft_,numberOfVarsRight_,SumBound_;
    INDEX upSize_,leftSize_,rightSize_,regSize_;
  };
  
  template<typename REPAM_ARRAY>
  REAL DiscreteTomographyFactorCounting::eval(INDEX up,INDEX left,INDEX right,const REPAM_ARRAY& repam){
    assert(up < upSize_);
    assert(left < leftSize_);
    assert(right < rightSize_);
    assert(repam.size() == (upSize_ + leftSize_ + rightSize_ + regSize_));

    auto xa = [&](INDEX idx){ return idx % numberOfLabels_;  };
    auto xb = [&](INDEX idx){ idx = (idx - xa(idx))/numberOfLabels_; return xa(idx); };
    auto z = [&](INDEX idx){ idx = (idx - xa(idx))/numberOfLabels_; return (idx - xa(idx))/numberOfLabels_; };
    
    if( xa(up) == xa(left) &&
        xb(up) == xb(right)&&
        z(left) + z(right) == z(up) )
      { return
          repam[up] +
          repam[upSize_ + left] +
          repam[upSize_ + leftSize_ + right] +
          repam[upSize_ + leftSize_ + rightSize_ + xb(left) + xa(right)*pow(numberOfLabels_,2)]; }
    else
      { return std::numeric_limits<REAL>::infinity(); }
    
  }

  inline INDEX DiscreteTomographyFactorCounting::size() const {
    return getSize(DiscreteTomographyFactorCounting::NODE::up) 
      + getSize(DiscreteTomographyFactorCounting::NODE::left) 
      + getSize(DiscreteTomographyFactorCounting::NODE::right) 
      + getSize(DiscreteTomographyFactorCounting::NODE::reg); 
  }

  
  inline INDEX DiscreteTomographyFactorCounting::getSize(NODE n) const {
    if( n == NODE::left ){ return leftSize_;  }
    if( n == NODE::right ){ return rightSize_;  }
    if( n == NODE::up ){ return upSize_;  }
    if( n == NODE::reg ){ return regSize_;  }
    return 0;
  }
   
  inline DiscreteTomographyFactorCounting::DiscreteTomographyFactorCounting(INDEX numberOfLabels,INDEX numberOfVarsLeft,INDEX numberOfVarsRight,INDEX SumBound)
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


  template<typename REPAM_ARRAY>
  inline void DiscreteTomographyFactorCounting::ReduceLp(LpInterfaceAdapter* lp,const REPAM_ARRAY& repam) const {
    assert(repam.size() == (upSize_+leftSize_+rightSize_+regSize_));

    auto xa = [&](INDEX idx){ return idx % numberOfLabels_;  };
    auto xb = [&](INDEX idx){ idx = (idx - xa(idx))/numberOfLabels_; return xa(idx); };
    auto xz = [&](INDEX idx){ idx = (idx - xa(idx))/numberOfLabels_; return (idx - xa(idx))/numberOfLabels_; };

    INDEX z_up_size = upSize_/pow(numberOfLabels_,2);
    REAL epsi = lp->GetEpsilon();
    REAL lb = LowerBound(repam);
    
    for(INDEX i=0;i<leftSize_;i++){
      INDEX lz = xz(i);
      for(INDEX j=0;j<rightSize_ && lz + xz(j) < z_up_size;j++){
        INDEX a = xa(i);
        INDEX b = xb(i);
        INDEX c = xa(j);
        INDEX d = xb(j);

        INDEX rz = xz(j);
        INDEX uz = lz+rz;

        INDEX l = a + b*numberOfLabels_ + lz*pow(numberOfLabels_,2);
        INDEX r = c + d*numberOfLabels_ + rz*pow(numberOfLabels_,2);
        INDEX u = a + d*numberOfLabels_ + uz*pow(numberOfLabels_,2);
        INDEX p = b + c*numberOfLabels_;
        
        assert(l < leftSize_);
        assert(r < rightSize_);
        assert(u < upSize_);
        assert(p < pow(numberOfLabels_,2));
        
        REAL value = repam[u] + repam[upSize_+l] + repam[upSize_+leftSize_+r] + repam[upSize_+leftSize_+rightSize_+p];
        if(value >= lb + epsi){
          lp->SetVariableBound(lp->GetAuxVariable(i + j*leftSize_),0.0,0.0);
        }        
      }
    }    
  }
      
  inline void DiscreteTomographyFactorCounting::CreateConstraints(LpInterfaceAdapter* lp) const {
    REAL inf = std::numeric_limits<REAL>::infinity();
    
    auto xa = [&](INDEX idx){ return idx % numberOfLabels_;  };
    auto xb = [&](INDEX idx){ idx = (idx - xa(idx))/numberOfLabels_; return xa(idx); };
    auto xz = [&](INDEX idx){ idx = (idx - xa(idx))/numberOfLabels_; return (idx - xa(idx))/numberOfLabels_; };
   
    LinExpr lhs_all = lp->CreateLinExpr();
    std::vector<LinExpr> lhs_up(upSize_);//,lp->CreateLinExpr());
    std::vector<LinExpr> lhs_left(leftSize_);//,lp->CreateLinExpr());
    std::vector<LinExpr> lhs_right(rightSize_);//,lp->CreateLinExpr());
    std::vector<LinExpr> lhs_reg(regSize_);//,lp->CreateLinExpr());
    
    auto InitVector = [&](std::vector<LinExpr>& v){
      for(INDEX i=0;i<v.size();i++){
        v[i] = lp->CreateLinExpr();
      }
    };
    InitVector(lhs_up);
    InitVector(lhs_left);
    InitVector(lhs_right);
    InitVector(lhs_reg);
    
    INDEX z_max = upSize_/pow(numberOfLabels_,2);
    
    for(INDEX i=0;i<leftSize_;i++){
      for(INDEX j=0;j<rightSize_;j++){
        LpVariable var;
        {
          // up variable constraint
          // sum_{ i(z) + j(z) = k(z) && i(a) = k(a) && j(b) = k(b) } eta(i,j) = mu_u(k)
          INDEX z = xz(i) + xz(j);
          if( z < z_max ){
            var = lp->GetAuxVariable(i + j*leftSize_);
            lp->SetVariableBound(var,0.0,1.0);
            
            INDEX a = xa(i);
            INDEX b = xb(j);
            INDEX idx = a + b*numberOfLabels_ + z*pow(numberOfLabels_,2);
            assert(idx < upSize_);
            lhs_up[idx] += var;
          } else {
            continue;
          }
        }
        {
          // sum_j eta(i,j) = mu_l(i) 
          lhs_left[i] += var;
          
          // sum_i eta(i,j) = mu_r(j) 
          lhs_right[j] += var;
        }
        {
          // sum_{i,j} eta(i,j) = 1
          lhs_all += var;
        }
        {
          // pairwise potential
          // sum_{i(b) = k(a) && j(a) = k(b)} eta(i,j) = mu_p(k)
          INDEX a = xb(i);
          INDEX b = xa(j);
          
          lhs_reg[a+b*numberOfLabels_] += var;
        }
      }
    }
    
    auto AddAllConstraints = [&]( std::vector<LinExpr> lhs,
                                  INDEX offset){
      for(INDEX i=0;i<lhs.size();i++){
        LinExpr rhs = lp->CreateLinExpr();
        rhs += lp->GetVariable(offset + i);
        lp->addLinearEquality(lhs[i],rhs);
      }
    };
  
    LinExpr rhs_all = lp->CreateLinExpr();
    rhs_all += 1;
    lp->addLinearEquality(lhs_all,rhs_all);
    
    AddAllConstraints(lhs_up,0);
    AddAllConstraints(lhs_left,upSize_);
    AddAllConstraints(lhs_right,upSize_+leftSize_);
    AddAllConstraints(lhs_reg,upSize_+leftSize_+rightSize_);
    
  }
  
  template<typename REPAM_ARRAY>
  REAL DiscreteTomographyFactorCounting::LowerBound(const REPAM_ARRAY& repam) const{
    assert(repam.size() == (upSize_ + leftSize_ + rightSize_ + regSize_));

    const INDEX noLabels = numberOfLabels_;
    const INDEX up_size = upSize_/pow(noLabels,2);
    const INDEX left_size = leftSize_/pow(noLabels,2);
    const INDEX right_size = rightSize_/pow(noLabels,2);
    
    REAL bound = std::numeric_limits<REAL>::infinity();
    
    if( 5000 < upSize_ ){
      auto op = [&](INDEX i,INDEX j){ return (i+j < up_size) ? i+j : up_size;  }; // 0 <= i+j < up_size
      
      for(INDEX i=0;i<pow(noLabels,4);i++){
        INDEX idx = i;
        INDEX a = idx % noLabels;
        idx = ( idx - a )/noLabels;
        INDEX b = idx % noLabels;
        idx = ( idx - b )/noLabels;
        INDEX c = idx % noLabels;
        idx = ( idx - c )/noLabels;
        INDEX d = idx % noLabels;

        auto z_up = [&](INDEX k){
          return repam[a + d*noLabels + k*pow(noLabels,2)];  };
        auto z_left = [&](INDEX k){
          return repam[up_size*pow(noLabels,2) + a + b*noLabels + k*pow(noLabels,2)];  };
        auto z_right = [&](INDEX k){
          return repam[up_size*pow(noLabels,2) + left_size*pow(noLabels,2) + c + d*noLabels + k*pow(noLabels,2)];  };

        REAL reg = repam[up_size*pow(noLabels,2) + left_size*pow(noLabels,2) + right_size*pow(noLabels,2) + b + c*noLabels];
        assert(reg > -std::numeric_limits<REAL>::max());

        MinConv mc(z_left,z_right,left_size,right_size,up_size);
        mc.CalcConv(op,z_left,z_right);

        for(INDEX k=0;k<up_size;k++){
          assert(k == op(mc.getIdxA(k),mc.getIdxB(k)));
          assert(k == (mc.getIdxA(k) + mc.getIdxB(k)));
          assert(!std::isnan(reg));

          REAL val = mc.getConv(k) + z_up(k)  + reg;
          
          assert(!std::isnan(val));
          assert(val > -eps);
          bound = std::min(bound,val);
          assert(bound > -eps);
        }
      }
    } else{
      for(INDEX lz=0;lz<left_size;lz++){
        for(INDEX rz=0;rz<right_size && rz+lz<up_size;rz++){
          for(INDEX i=0;i<pow(noLabels,4);i++){
            INDEX idx = i;
            INDEX a = idx % noLabels;
            idx = (idx - a)/noLabels;
            INDEX b = idx % noLabels;
            idx = (idx - b)/noLabels;
            INDEX c = idx % noLabels;
            INDEX d = (idx - c)/noLabels;

            assert(a + d*noLabels + (lz+rz)*pow(noLabels,2) < up_size*pow(noLabels,2));
            assert(a + b*noLabels + lz*pow(noLabels,2) < left_size*pow(noLabels,2));
            assert(c + d*noLabels + rz*pow(noLabels,2) < right_size*pow(noLabels,2));
            assert(b + c*noLabels < pow(noLabels,2));
            
            REAL value = repam[a + d*noLabels + (lz+rz)*pow(noLabels,2)];
            value += repam[up_size*pow(noLabels,2) + a + b*noLabels + lz*pow(noLabels,2)];
            value += repam[up_size*pow(noLabels,2) + left_size*pow(noLabels,2) + c + d*noLabels + rz*pow(noLabels,2)];
            value += repam[up_size*pow(noLabels,2) + left_size*pow(noLabels,2) + right_size*pow(noLabels,2) + b + c*noLabels];

            assert(value > -eps);
            bound = std::min(bound,value);
          }
        }
      }
    }
    return bound;
  }

  /*
  template<typename REPAM_ARRAY>
  inline void DiscreteTomographyFactorCounting::NarrowPrimal(const REPAM_ARRAY& repam, const REAL epsilon, typename PrimalSolutionStorage::Element primal)
  {
     const REAL min_val = LowerBound(repam);

     // go over all labelings and set to false all entries that lead to labelings with value greater than min_val + epsilon
     // do it as follows: given some entry in repam, go over all combinations such that it is used, and see whether it always leads to labeling costs larger than min_val + epsilon. If so, set primal in the corresponding place to false

     INDEX sum_up = std::accumulate(repam.begin(),repam.begin()+upSize_,0);
     INDEX sum_left = std::accumulate(repam.begin()+upSize_,repam.begin()+upSize_+leftSize_,0);
     INDEX sum_right = std::accumulate(repam.begin()+upSize_+leftSize_,repam.begin()+upSize_+leftSize_+rightSize_,0);
     INDEX sum_reg = std::accumulate(repam.begin()+upSize_+leftSize_+rightSize_,repam.end(),0);

     // (i) up
     for(INDEX i=0; i<upSize_; ++i) {

     }
     // (ii) left
     // (iii) right
     // (iv) reg
     INDEX z_up_size = upSize_/pow(numberOfLabels_,2);
     INDEX z_left_size = leftSize_/pow(numberOfLabels_,2);
     INDEX z_right_size = rightSize_/pow(numberOfLabels_,2);

     auto op = [&](INDEX i,INDEX j){ return (i+j < z_up_size) ? i+j : z_up_size;  };
     for( INDEX i=0;i<pow(numberOfLabels_,4);i++ ){
        INDEX idx = i;
        INDEX a = idx % numberOfLabels_;
        idx = ( idx - a )/numberOfLabels_;
        INDEX b = idx % numberOfLabels_;
        idx = ( idx - b )/numberOfLabels_;
        INDEX c = idx % numberOfLabels_;
        idx = ( idx - c )/numberOfLabels_;
        INDEX d = idx % numberOfLabels_;

        auto z_up = [&](INDEX k){ return repam[a + d*numberOfLabels_ + k*pow(numberOfLabels_,2)]; };
        auto z_left = [&](INDEX k){ return repam[upSize_ + a + b*numberOfLabels_ + k*pow(numberOfLabels_,2)];  };
        auto z_right = [&](INDEX k){ return repam[upSize_ + leftSize_ + c + d*numberOfLabels_ + k*pow(numberOfLabels_,2)];  };

        REAL reg = repam[upSize_ + leftSize_ + rightSize_ + b + c*numberOfLabels_];
        assert(!std::isnan(reg));
        assert(reg > -std::numeric_limits<REAL>::max() );

        MinConv mc(z_left,z_right,z_left_size,z_right_size,z_up_size);
        mc.CalcConv(op,z_left,z_right);

        REAL m_new = 0;
        for( INDEX j=0;j<z_up_size;j++ ){
           assert(j == op(mc.getIdxA(j),mc.getIdxB(j)));
           assert(z_up(j) > -std::numeric_limits<REAL>::max() );
           assert(!std::isnan(z_up(j)));
	
        m_new = mc.getConv(j)+z_up(j)+reg;
        if(m_new > min_val + epsilon) {
           //primal[] = false;
        }
      }
    }
  }
  */

  inline void DiscreteTomographyFactorCounting::PropagatePrimal(PrimalSolutionStorage::Element primal) const{
    
    struct IdxLbl {
      IdxLbl(INDEX x,INDEX y,INDEX l) : idx(x),c(y),noLabels(l) {
        INDEX i=idx;
        a = i % noLabels;
        i = (i-a) / noLabels;
        b = i % noLabels;
        z = (i-b) / noLabels;
      }
      
      const INDEX idx,c,noLabels;
      INDEX a = 0;
      INDEX b = 0;
      INDEX z = 0;
    };

    auto getOptLabel = [&](INDEX s,INDEX t){
      INDEX noTrue = 0;
      INDEX noUnkwn = 0;
      INDEX opt = 0;
      for(INDEX i=0;i<s;i++){
        if( primal[t + i] == true ){ noTrue++; opt=i; }
        if( primal[t + i] == unknownState ){ noUnkwn++; opt=i; }
      }
      assert(noTrue <=1);
      assert(noTrue != 1 || noUnkwn == 0);
      assert(noTrue == 0 || noUnkwn == 0);

      if( noTrue == 1 ){ return IdxLbl(opt,1,numberOfLabels_); }
      else if( noUnkwn == 1 && noTrue == 0  ){
        primal[t+opt]=true;
        return IdxLbl(opt,1,numberOfLabels_);
      }
      else if( noUnkwn > 1 && noTrue == 0 ){
        return IdxLbl(opt,0,numberOfLabels_);
      }
      else{ assert(false); return IdxLbl(opt,0,numberOfLabels_);  } // not possible!
    };
    
    /* Get primal label */
    auto up    = getOptLabel(upSize_   ,0);
    auto left  = getOptLabel(leftSize_ ,upSize_);
    auto right = getOptLabel(rightSize_,upSize_ + leftSize_);
    auto reg   = getOptLabel(regSize_  ,upSize_ + leftSize_ + rightSize_);

    auto CalcIdx = [&](INDEX a,INDEX b,INDEX z){ return a + b*numberOfLabels_ + z*((INDEX)pow(numberOfLabels_,2)); };
    auto Set2False = [&](INDEX s,INDEX t){ for(INDEX i=0;i<s;i++){primal[t+i]=false;} };
    INDEX count = up.c + left.c + right.c + reg.c;
    
    /* If exact one label is missing, we can calculate it */
    if( count == 3 ){
      if(up.c == 0 && left.b == reg.a && right.a == reg.b
         && CalcIdx(left.a,right.b,left.z+right.z) < upSize_ ){   // up + consistency
        Set2False(upSize_,0);
        primal[CalcIdx(left.a,right.b,left.z+right.z)]=true; }
      if(left.c == 0 && up.b == right.b && right.a == reg.b
         && up.z >= right.z
         && CalcIdx(up.a,reg.a,up.z-right.z) < leftSize_){ // left + consistency
        Set2False(leftSize_,upSize_);
        primal[upSize_ + CalcIdx(up.a,reg.a,up.z-right.z)]=true; }
      if(right.c == 0 && up.a == left.a && left.b == reg.a   // right + consistency
         && CalcIdx(up.b,reg.b,up.z-left.z)
         && up.z >= right.z ){
        Set2False(rightSize_,upSize_ + leftSize_);
        primal[upSize_ + leftSize_ + CalcIdx(up.b,reg.b,up.z-left.z)]=true; }
      if(reg.c == 0 && up.a == left.a && up.b == right.b 
         && up.z == (left.z + right.z) ){                    // pairwise + consistency
        Set2False(regSize_,upSize_ + leftSize_ + rightSize_);
        primal[upSize_ + leftSize_ + rightSize_ + CalcIdx(left.b,right.a,0)]=true; }
    } 
    /* Just for debugging */
    if( count == 4){
      assert(left.a == up.a && up.b == right.b);
      assert(left.b == reg.a && reg.b == right.a);
      assert(up.z == (left.z + right.z));
    }
    
  }

  template<typename REPAM_ARRAY>
  REAL DiscreteTomographyFactorCounting::EvaluatePrimal(const REPAM_ARRAY& repam, const PrimalSolutionStorage::Element primal) const{
    assert(repam.size() == (upSize_ + leftSize_ + rightSize_ + regSize_));
    REAL val = 0;
    INDEX count = 0;
    INDEX up = 0;
    INDEX left = 0;
    INDEX right = 0;
    INDEX reg = 0;

    auto xa = [&](INDEX idx){ return idx % numberOfLabels_;  };
    auto xb = [&](INDEX idx){ idx = (idx - xa(idx))/numberOfLabels_; return xa(idx); };
    auto z = [&](INDEX idx){ idx = (idx - xa(idx))/numberOfLabels_; return (idx - xa(idx))/numberOfLabels_; };
    
    auto updateVal = [&](INDEX s,INDEX t,INDEX& idx){
      for(INDEX i=0;i<s;i++){
        if(primal[t+i] == true){
          val += repam[t+i];
          idx = i;
          count++;
        }
      }
    };

    INDEX size = 0;      updateVal(upSize_,size,up);
    size += upSize_;     updateVal(leftSize_,size,left);
    size += leftSize_;   updateVal(rightSize_,size,right);
    size += rightSize_;  updateVal(regSize_,size,reg);
    
    assert(count <= 4);

    if( xa(up) == xa(left) &&
        xb(up) == xb(right)&&
        z(left) + z(right) == z(up) &&
        xa(reg) == xb(left) &&
        xb(reg) == xa(right) &&
        count == 4 )
      {  return val; }
    else{
      return std::numeric_limits<REAL>::infinity();
    }
      
  }
    
} // end namespace LP_MP

#endif // LPMP_DTOMOGRAPHY_FACTOR_COUNTING_HXX

