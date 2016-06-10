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

    template<typename REPAM_ARRAY>
    void MaximizePotentialAndComputePrimal(const REPAM_ARRAY& repam, PrimalSolutionStorage::Element primal) const; 
    
    template<typename REPAM_ARRAY>
    REAL EvaluatePrimal(const REPAM_ARRAY& repam, const PrimalSolutionStorage::Element primal) const; //--required

    void WritePrimal(const PrimalSolutionStorage::Element, std::ofstream& fs) const; //--required

    // --------------------
    // custom public methods 
    
    INDEX getSize(NODE);
    REAL eval(INDEX,INDEX,INDEX);
    
  private:

    const INDEX numberOfLabels_,numberOfVarsLeft_,numberOfVarsRight_,SumBound_;
    INDEX upSize_,leftSize_,rightSize_,regSize_;
  };

  REAL DiscreteTomographyFactorCounting::eval(INDEX up,INDEX left,INDEX right){
    assert(up < upSize_);
    assert(left < leftSize_);
    assert(right < rightSize_);

    auto xa = [&](INDEX idx){ return idx % numberOfLabels_;  };
    auto xb = [&](INDEX idx){ idx = (idx - xa(idx))/numberOfLabels_; return xa(idx); };
    auto z = [&](INDEX idx){ idx = (idx - xa(idx))/numberOfLabels_; return xb(idx); };
    
    if( xa(up) == xa(left) &&
	xb(up) == xb(right)&&
	z(left) + z(right) == z(up) )
      { return 0; }
    else
      { return std::numeric_limits<REAL>::max(); }
    
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
    REAL m = std::numeric_limits<REAL>::max();

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

      auto z_up = [&](INDEX k){ return repam[a + d*numberOfLabels_ + k*pow(numberOfLabels_,2)];  };
      auto z_left = [&](INDEX k){ return repam[upSize_ + a + b*numberOfLabels_ + k*pow(numberOfLabels_,2)];  };
      auto z_right = [&](INDEX k){ return repam[upSize_ + leftSize_ + c + d*numberOfLabels_ + k*pow(numberOfLabels_,2)];  };
      REAL reg = repam[upSize_ + leftSize_ + rightSize_ + b + c*numberOfLabels_];

      MinConv mc(z_left,z_right,z_left_size,z_right_size,z_up_size);
      mc.CalcConv(op,z_left,z_right);

      REAL m_new = 0;
      for( INDEX j=0;j<z_up_size;j++ ){
	//std::cout << j << " " << mc.getIdxA(j) << " " << mc.getIdxB(j) << std::endl;
	assert(j == op(mc.getIdxA(j),mc.getIdxB(j)));
	m_new = mc.getConv(j)+z_up(j)+reg;
	if( m > m_new ){
	  m = m_new;
	}
      }
    }
    return m;
  }

  template<typename REPAM_ARRAY>
  void DiscreteTomographyFactorCounting::MaximizePotentialAndComputePrimal(const REPAM_ARRAY& repam, PrimalSolutionStorage::Element primal) const{
    assert(repam.size() == (upSize_ + leftSize_ + rightSize_ + regSize_));
    for(INDEX i=0; i<repam.size(); ++i) { 
      assert(primal[i] == true);
      primal[i]=false;
    }
    
    REAL m = std::numeric_limits<REAL>::max();
    INDEX up,left,right,reg;
    
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

      auto z_up = [&](INDEX k){ return repam[a + d*numberOfLabels_ + k*pow(numberOfLabels_,2)];  };
      auto z_left = [&](INDEX k){ return repam[upSize_ + a + b*numberOfLabels_ + k*pow(numberOfLabels_,2)];  };
      auto z_right = [&](INDEX k){ return repam[upSize_ + leftSize_ + c + d*numberOfLabels_ + k*pow(numberOfLabels_,2)];  };
      REAL reg = repam[upSize_ + leftSize_ + rightSize_ + b + c*numberOfLabels_];

      MinConv mc(z_left,z_right,z_left_size,z_right_size,z_up_size);
      mc.CalcConv(op,z_left,z_right);

      REAL m_new = 0;
      for( INDEX j=0;j<z_up_size;j++ ){
	assert(j == op(mc.getIdxA(j),mc.getIdxB(j)));
	m_new = mc.getConv(j)+z_up(j)+reg;
	if( m > m_new ){
	  m = m_new;
	  up = a + d*numberOfLabels_ + j*pow(numberOfLabels_,2);
	  left = upSize_ + a + b*numberOfLabels_ + mc.getIdxA(j)*pow(numberOfLabels_,2);
	  right = upSize_ + leftSize_ + c + d*numberOfLabels_ + mc.getIdxB(j)*pow(numberOfLabels_,2);
	  reg = upSize_ + leftSize_ + rightSize_ + b + c*numberOfLabels_;
	}
      }
    }

    primal[up]=true;
    primal[left]=true;
    primal[right]=true;
    primal[reg]=true;
    
  }
  
  /*
    template<typename REPAM_ARRAY>
    INDEX DiscreteTomographyFactorCounting::ComputeOptimalLabeling(const REPAM_ARRAY& repam) const{
    assert(repam.size() == (upSize_ + leftSize_ + rightSize_ + regSize_));

    return 0;
    }
  */

  template<typename REPAM_ARRAY>
  REAL DiscreteTomographyFactorCounting::EvaluatePrimal(const REPAM_ARRAY& repam, const PrimalSolutionStorage::Element primal) const {
    assert(repam.size() == (upSize_ + leftSize_ + rightSize_ + regSize_));

    
    
    return 0;
  }

  void DiscreteTomographyFactorCounting::WritePrimal(const PrimalSolutionStorage::Element primal, std::ofstream& fs) const {
    
  }
  
} // end namespace LP_MP

#endif // LPMP_DTOMOGRAPHY_FACTOR_COUNTING_HXX

