#ifndef LPMP_DTOMOGRAPHY_FACTOR_COUNTING_HXX
#define LPMP_DTOMOGRAPHY_FACTOR_COUNTING_HXX

#include "LP_MP.h"
#include "minConv.hxx"
#include <math.h> 

namespace LP_MP{

  using MinConv = discrete_tomo::MinConv<std::function<REAL(INDEX)>,REAL,INDEX>;
  
  class DiscreteTomographyFactorCounting{

  public:
    
    //DiscreteTomographyFactorCounting(); //--required
    DiscreteTomographyFactorCounting(INDEX numberOfLabels,INDEX,INDEX,INDEX,INDEX);
    
    template<typename REPAM_ARRAY>
    INDEX ComputeOptimalLabeling(const REPAM_ARRAY& repam) const; //--required

    template<typename REPAM_ARRAY>
    REAL LowerBound(const REPAM_ARRAY& repam) const; //--required

    template<typename REPAM_ARRAY>
    REAL EvaluatePrimal(const REPAM_ARRAY& repam, const PrimalSolutionStorage::Element primal) const; //--required

    void WritePrimal(const PrimalSolutionStorage::Element, std::ofstream& fs) const; //--required

    // --------------------
    // custom public methods 
    
    void setRHS(INDEX b){ root_ = true; rhs_ = b; }
    void setNOISE(std::function<REAL (REAL)> n){ noise_ = true; weights_ = n;  }
    REAL eval(INDEX,INDEX,INDEX);
    
  private:

    const INDEX a_,b_,c_,d_;
    const INDEX numberOfLabels_;
    INDEX upSize_,leftSize_,rightSize_,regSize_;
    INDEX rhs_;
    bool root_ = false;
    bool noise_ = false;
    std::function<REAL (REAL)> weights_ = NULL;
        
  };

  REAL DiscreteTomographyFactorCounting::eval(INDEX left,INDEX right,INDEX up){
    REAL val = std::numeric_limits<REAL>::max();
    
  }
  
  DiscreteTomographyFactorCounting::DiscreteTomographyFactorCounting(INDEX numberOfLabels,INDEX a,INDEX b, INDEX c, INDEX d)
    : numberOfLabels_(numberOfLabels),a_(a),b_(b),c_(c),d_(d){
    assert(0 <= a_); assert(a_ <= b_); assert(b_ <= c_); assert(c_ <= d_);
    assert(numberOfLabels > 1);
    
    upSize_ = pow(numberOfLabels_,2)*((d_-a_+1)*(numberOfLabels_-1)+1);
    leftSize_ = pow(numberOfLabels_,2)*((b_-a_+1)*(numberOfLabels_-1)+1);
    rightSize_ = pow(numberOfLabels_,2)*((d_-c_+1)*(numberOfLabels_-1)+1);
    
    regSize_ = pow(numberOfLabels_,2);
  }

  template<typename REPAM_ARRAY>
  REAL DiscreteTomographyFactorCounting::LowerBound(const REPAM_ARRAY& repam) const{
    assert(repam.size() == (upSize_ + leftSize_ + rightSize_ + regSize_));
    assert( !noise_ || root_ );
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
      auto z_reg = repam[upSize_ + leftSize_ + rightSize_ + b + c*numberOfLabels_];

      MinConv mc(z_left,z_right,
		 ((b_-a_+1)*(numberOfLabels_-1)+1),((c_-d_+1)*(numberOfLabels_-1)+1),((d_-a_+1)*(numberOfLabels_-1)+1));
      mc.CalcConv(op);

      REAL m_new = 0;
      if( noise_ && root_ ){
	for( INDEX j=0;j<((d_-a_+1)*(numberOfLabels_-1)+1);j++ ){
	  m_new = mc.getConv(j)+z_up(j)+z_reg+weights_(std::abs(rhs_-j));
	  if( m > m_new ){
	    m = m_new;
	  }
	}
      }
      else if( !noise_ && root_ ){
	m_new = mc.getConv(rhs_)+z_up(rhs_)+z_reg;
	if( m > m_new ){
	  m = m_new;
	}	
      }
      else{
	for( INDEX j=0;j<((d_-a_+1)*(numberOfLabels_-1)+1);j++ ){
	  m_new = mc.getConv(j)+z_up(j)+z_reg;
	  if( m > m_new ){
	    m = m_new;
	  }
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

