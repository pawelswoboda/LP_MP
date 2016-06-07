#ifndef LPMP_DTOMOGRAPHY_FACTOR_COUNTING_HXX
#define LPMP_DTOMOGRAPHY_FACTOR_COUNTING_HXX

#include "LP_MP.h"
#include "minConv.hxx"
#include <math.h> 

namespace LP_MP{

  using MinConv = discrete_tomo::MinConv<std::function<REAL(INDEX)>,REAL,INDEX>;

  class DiscreteTomographyFactorUnary{

  public:
    
    //DiscreteTomographyFactorCounting(); //--required
    DiscreteTomographyFactorUnary(INDEX numberOfLabels,INDEX numberOfVars);
    
    template<typename REPAM_ARRAY>
    INDEX ComputeOptimalLabeling(const REPAM_ARRAY& repam) const; //--required

    template<typename REPAM_ARRAY>
    REAL LowerBound(const REPAM_ARRAY& repam) const; //--required

    template<typename REPAM_ARRAY>
    REAL EvaluatePrimal(const REPAM_ARRAY& repam, const PrimalSolutionStorage::Element primal) const; //--required

    void WritePrimal(const PrimalSolutionStorage::Element, std::ofstream& fs) const; //--required

    INDEX getSize(){ return nodeSize_; }
    
  private:
    const INDEX numberOfVars_,numebrOfLabels_,nodeSize_;
  };

  DiscreteTomographyFactorUnary::DiscreteTomographyFactorUnary(INDEX numberOfLabels,INDEX numberOfVars)
    : numberOfLabels_(numberOfLabels),numberOfVars_(numberOfVars),nodeSize_(pow(numberOfLabels_,2)*(numberOfVars_*(numberOfLabels_-1)+1)){
    assert(numberOfLabels > 1);
    assert(numberOfVars >= 1);
  }

  template<typename REPAM_ARRAY>
  REAL DiscreteTomographyFactorUnary::LowerBound(const REPAM_ARRAY& repam) const{
    assert(repam.size()==nodeSize_);
    REAL m = std::numeric_limits<REAL>::max();
    REAL m_new = 0;
    if( numberOfVars_ == 1 ){
      for( INDEX i=0;i<numberOfLabels_;i++ ){
	m_new = repam[i + i*numberOfLabels_ + i*pow(numberOfLabels_,2)] + unary(i);
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
  INDEX DiscreteTomographyFactorUnary::ComputeOptimalLabeling(const REPAM_ARRAY& repam) const{
    return 0;
  }

  template<typename REPAM_ARRAY>
  REAL DiscreteTomographyFactorUnary::EvaluatePrimal(const REPAM_ARRAY& repam, const PrimalSolutionStorage::Element primal) const {
    return 0;
  }

  void DiscreteTomographyFactorUnary::WritePrimal(const PrimalSolutionStorage::Element, std::ofstream& fs) const {
    
  }
  
} // end namespace LP_MP

#endif // LPMP_DTOMOGRAPHY_FACTOR_COUNTING_HXX

