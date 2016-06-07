#ifndef LPMP_DTOMOGRAPHY_FACTOR_COUNTING_HXX
#define LPMP_DTOMOGRAPHY_FACTOR_COUNTING_HXX

#include "LP_MP.h"
#include "minConv.hxx"
#include <math.h> 

namespace LP_MP{

  using MinConv = discrete_tomo::MinConv<std::function<REAL(INDEX)>,REAL,INDEX>;
  
  class DiscreteTomographyFactorCounting{

  public:

    enum class NODE {left,right,up,reg};
    
    //DiscreteTomographyFactorCounting(); //--required
    DiscreteTomographyFactorCounting(INDEX numberOfLabels,INDEX numberOfVarsLeft,INDEX numberOfVarsRight);
    
    template<typename REPAM_ARRAY>
    INDEX ComputeOptimalLabeling(const REPAM_ARRAY& repam) const; //--required

    template<typename REPAM_ARRAY>
    REAL LowerBound(const REPAM_ARRAY& repam) const; //--required

    template<typename REPAM_ARRAY>
    REAL EvaluatePrimal(const REPAM_ARRAY& repam, const PrimalSolutionStorage::Element primal) const; //--required

    void WritePrimal(const PrimalSolutionStorage::Element, std::ofstream& fs) const; //--required

    // --------------------
    // custom public methods 
    
    INDEX getSize(NODE);
    
  private:

    const INDEX numberOfLabels_,numberOfVarsLeft_,numberOfVarsRight_;
    INDEX upSize_,leftSize_,rightSize_,regSize_;
        
  };

  INDEX DiscreteTomographyFactorCounting::getSize(NODE n){
    if( n == NODE::left ){ return leftSize_;  }
    if( n == NODE::right ){ return rightSize_;  }
    if( n == NODE::up ){ return upSize_;  }
    if( n == NODE::reg ){ return regSize_;  }
  }
   
  DiscreteTomographyFactorCounting::DiscreteTomographyFactorCounting(INDEX numberOfLabels,INDEX numberOfVarsLeft,INDEX numberOfVarsRight)
    : numberOfLabels_(numberOfLabels),numberOfVarsLeft_(numberOfVarsLeft),numberOfVarsRight_(numberOfVarsRight){
    assert(numberOfLabels_ > 1);
    assert(numberOfVarsLeft_ > 1);
    assert(numberOfVarsRight_ > 1);
    
    upSize_ = pow(numberOfLabels_,2)*((numberOfVarsLeft_ + numberOfVarsRight_)*(numberOfLabels_-1)+1);
    leftSize_ = pow(numberOfLabels_,2)*(numberOfVarsLeft_*(numberOfLabels_-1)+1);
    rightSize_ = pow(numberOfLabels_,2)*(numberOfVarsRight_*(numberOfLabels_-1)+1);
    
    regSize_ = pow(numberOfLabels_,2);
  }

  template<typename REPAM_ARRAY>
  REAL DiscreteTomographyFactorCounting::LowerBound(const REPAM_ARRAY& repam) const{
    assert(repam.size() == (upSize_ + leftSize_ + rightSize_ + regSize_));
    REAL m = std::numeric_limits<REAL>::max();

    auto op = [](INDEX i,INDEX j){ return i+j;  };
    for( INDEX i=0;i<pow(numberOfLabels_,4);i++ ){
      INDEX idx = i;
      INDEX a = idx % numberOfLabels_;
      idx = ( idx - a )/numberOfLabels_;
      INDEX b = idx % numberOfLabels_;
      idx = ( idx - b )/numberOfLabels_;
      INDEX c = idx % numberOfLabels_;
      idx = ( idx - c )/numberOfLabels_;
      INDEX d = idx % numberOfLabels_;

      auto z_up = [&](INDEX k){ return repam[c + d*numberOfLabels_ + k*pow(numberOfLabels_,2)];  };
      auto z_left = [&](INDEX k){ return repam[upSize_ + a + b*numberOfLabels_ + k*pow(numberOfLabels_,2)];  };
      auto z_right = [&](INDEX k){ return repam[upSize_ + leftSize_ + c + d*numberOfLabels_ + k*pow(numberOfLabels_,2)];  };
      REAL reg = repam[upSize_ + leftSize_ + rightSize_ + b + c*numberOfLabels_];

      MinConv mc(z_left,z_right,leftSize_,rightSize_);
      mc.CalcConv(op);

      REAL m_new = 0;
      for( INDEX j=0;j<upSize_;j++ ){
	assert(j == mc.getIdxA(j) + mc.getIdxB(j));
	m_new = mc.getConv(j)+z_up(j)+reg;
	if( m > m_new ){
	  m = m_new;
	}
      }
    }
    return m;
  }
    
  template<typename REPAM_ARRAY>
  INDEX DiscreteTomographyFactorCounting::ComputeOptimalLabeling(const REPAM_ARRAY& repam) const{
    return 0;
  }

  template<typename REPAM_ARRAY>
  REAL DiscreteTomographyFactorCounting::EvaluatePrimal(const REPAM_ARRAY& repam, const PrimalSolutionStorage::Element primal) const {
    return 0;
  }

  void DiscreteTomographyFactorCounting::WritePrimal(const PrimalSolutionStorage::Element, std::ofstream& fs) const {
    
  }
  
} // end namespace LP_MP

#endif // LPMP_DTOMOGRAPHY_FACTOR_COUNTING_HXX

