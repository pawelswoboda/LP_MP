#ifndef LPMP_DTOMOGRAPHY_FACTOR_COUNTING_HXX
#define LPMP_DTOMOGRAPHY_FACTOR_COUNTING_HXX

#include "LP_MP.h"
#include "minConv.hxx"
#include <math.h> 

namespace LP_MP{

  using MinConv = discrete_tomo::MinConv<REAL,INDEX>;
  
  class DiscreteTomographyFactorCounting{

  public:

    enum class NODE {left,right,up,reg};
    
    DiscreteTomographyFactorCounting(INDEX numberOfLabels,INDEX numberOfVarsLeft,INDEX numberOfVarsRight,INDEX SumBound);
    
    template<typename REPAM_ARRAY>
    INDEX ComputeOptimalLabeling(const REPAM_ARRAY& repam) const { }; //--required

    template<typename REPAM_ARRAY>
    REAL LowerBound(const REPAM_ARRAY& repam) const; //--required

    //template<typename REPAM_ARRAY>
    //void MaximizePotentialAndComputePrimal(const REPAM_ARRAY& repam, PrimalSolutionStorage::Element primal) const;
    
    void PropagatePrimal(PrimalSolutionStorage::Element primal) const;

    template<typename REPAM_ARRAY>
    REAL EvaluatePrimal(const REPAM_ARRAY& repam, const PrimalSolutionStorage::Element primal) const; //--required

    //void WritePrimal(const PrimalSolutionStorage::Element, std::ofstream& fs) const; //

    // --------------------
    // custom public methods 
    
    INDEX getSize(NODE);

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
    auto z = [&](INDEX idx){ idx = (idx - xa(idx))/numberOfLabels_; return xb(idx); };
    
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
  
  INDEX DiscreteTomographyFactorCounting::getSize(NODE n){
    if( n == NODE::left ){ return leftSize_;  }
    if( n == NODE::right ){ return rightSize_;  }
    if( n == NODE::up ){ return upSize_;  }
    if( n == NODE::reg ){ return regSize_;  }
    return 0;
  }
   
  DiscreteTomographyFactorCounting::DiscreteTomographyFactorCounting(INDEX numberOfLabels,INDEX numberOfVarsLeft,INDEX numberOfVarsRight,INDEX SumBound)
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
  REAL DiscreteTomographyFactorCounting::LowerBound(const REPAM_ARRAY& repam) const{
    assert(repam.size() == (upSize_ + leftSize_ + rightSize_ + regSize_));
    REAL m = std::numeric_limits<REAL>::infinity();
    
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
        m = std::min(m,m_new);
        assert(m > -1.0e-02);
      }
    }
    return m;
  }


  void DiscreteTomographyFactorCounting::PropagatePrimal(PrimalSolutionStorage::Element primal) const{

    struct IdxLbl {
      IdxLbl(INDEX x,INDEX y) : idx(x),c(y) {
        INDEX i=idx;
        a = i & numberOfLabels_;
        i = (i-a) / numberOfLabels_;
        b = i % numberOfLabels_;
        z = (i-b) / numberOfLabels_;
      }
      
      const INDEX idx,c;
      INDEX a = 0;
      INDEX b = 0;
      INDEX z = 0;
    };
    
    auto getOptLabel = [&](INDEX s,INDEX t){
      INDEX noTrue = 0;
      INDEX noUnkwn = 0;
      INDEX opt = 0;
      for(INDEX i=0;i<s;i++){
        if( primal[t + i] == true ){ noTrue++; opt=i;  }
        if( primal[t + i] == unknownState ){ noUnkwn++; }
      }
      assert(noTrue <=1);
      assert(noTrue != 1 || noUnkwn == 0);
      assert(noUnkwn != 0 || noTrue == 1 );
      return IdxLbl(opt,noTrue);
    };
    
    auto up    = getOptLabel(upSize_   ,0);
    auto left  = getOptLabel(leftSize_ ,upSize_);
    auto right = getOptLabel(rightSize_,upSize_ + leftSize_);
    auto reg   = getOptLabel(regSize_  ,upSize_ + leftSize_ + rightSize_);

    INDEX count = up.c + left.c + right.c + reg.c;
    if( count == 3 ){
      //TODO
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

    auto xa = [&](INDEX idx){ return idx % numberOfLabels_;  };
    auto xb = [&](INDEX idx){ idx = (idx - xa(idx))/numberOfLabels_; return xa(idx); };
    auto z = [&](INDEX idx){ idx = (idx - xa(idx))/numberOfLabels_; return xb(idx); };

    if( xa(up) == xa(left) &&
        xb(up) == xb(right)&&
        z(left) + z(right) == z(up) &&
        xa(reg) == xb(left) &&
        xb(reg) == xa(right) &&
        count == 4 )
      { return val; }
    else{
      return std::numeric_limits<REAL>::infinity();
    }
      
  }
    
} // end namespace LP_MP

#endif // LPMP_DTOMOGRAPHY_FACTOR_COUNTING_HXX

