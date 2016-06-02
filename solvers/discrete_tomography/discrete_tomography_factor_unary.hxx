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
    DiscreteTomographyFactorCounting(INDEX numberOfLabels,INDEX a,INDEX b);
    
    template<typename REPAM_ARRAY>
    INDEX ComputeOptimalLabeling(const REPAM_ARRAY& repam) const; //--required

    template<typename REPAM_ARRAY>
    REAL LowerBound(const REPAM_ARRAY& repam) const; //--required

    template<typename REPAM_ARRAY>
    REAL EvaluatePrimal(const REPAM_ARRAY& repam, const PrimalSolutionStorage::Element primal) const; //--required

    void WritePrimal(const PrimalSolutionStorage::Element, std::ofstream& fs) const; //--required

    void setFUNC(std::function<REAL (INDEX)> unary){ unary_ = unary; }
  private:
    std::function<REAL (INDEX)> unary_;
    INDEX a_,b_;
  };

  DiscreteTomographyFactorCounting::DiscreteTomographyFactorCounting(INDEX numberOfLabels,INDEX a,INDEX b)
    : numberOfLabels_(numberOfLabels){
    assert(numberOfLabels > 1);
    assert(a_ <= b_);
    nodeSize_ = pow(numberOfLabels_,2)*((b_-a_+1)*(numberOfLabels_-1)+1);
  }

  template<typename REPAM_ARRAY>
  REAL DiscreteTomographyFactorCounting::LowerBound(const REPAM_ARRAY& repam) const{
    assert(repam.size()==nodeSize_);
    REAL m = std::numeric_limits<REAL>::max();
    REAL m_new = 0;
    if( a_ == b_ ){
      for( INDEX i=0;i<numberOfLabels_;i++ ){
	m_new = repam[i + i*numberOfLabels_ + i*pow(numberOfLabels_,2)] + unary_(i);
	if( m > m_new ){ m = m_new;  }
      }
    }
    else{
      for( INDEX i=0;i<nodeSize_;i++ ){
	if( m > repam[i] ){ m = repam[i]; }
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

