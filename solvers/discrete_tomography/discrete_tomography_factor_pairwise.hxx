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

    void setREG(std::function<REAL (INDEX,INDEX)> reg){ reg_ = reg;  }
  private:
    std::function<REAL (INDEX,INDEX)> reg_;
    INDEX a_,b_,c_,d_;
  };

  DiscreteTomographyFactorCounting::DiscreteTomographyFactorCounting(INDEX numberOfLabels,INDEX a,INDEX b, INDEX c, INDEX d)
    : numberOfLabels_(numberOfLabels){
    assert(numberOfLabels > 1);
    assert(a_ <= b_); assert(c_ <= d_);
    leftSize_ = pow(numberOfLabels_,2)*((b_-a_+1)*numberOfLabels_+1);
    rightSize_ = pow(numberOfLabels_,2)*((c_-b_+1)*numberOfLabels_+1);
  }

  template<typename REPAM_ARRAY>
  REAL DiscreteTomographyFactorCounting::LowerBound(const REPAM_ARRAY& repam) const{
    assert(repam.size()==leftSize_ + rightSize_ + pow(numberOfLabels_,2);
    REAL m = std::numeric_limits<REAL>::max();
    REAL m_new = 0;
    for( INDEX k=0;k<pow(numberOfLabels_,2);k++ ){
      INDEX i = k % numberOfLabels_;
      INDEX j = ((k-i)/numberOfLabels_) % numberOfLabels_;
      m_new = reg_(i,j)
	+ repam[i + i*numberOfLabels_ + i*pow(numberOfLabels_,2)]
	+ repam[leftSize_ + i + i*numberOfLabels_ + i*pow(numberOfLabels_,2)]
	+ repam[LeftSize_ + rightSize_ + i + j*numberOfLabels_];
      if( m > m_new ){ m = m_new;  }
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

