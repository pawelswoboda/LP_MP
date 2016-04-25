#ifndef LPMP_DTOMOGRAPHY_FACTOR_COUNTING_HXX
#define LPMP_DTOMOGRAPHY_FACTOR_COUNTING_HXX

#include "LP_MP.h"
#include <math.h> 

namespace LP_MP{

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
    
  private:

    const INDEX a_,b_,c_,d_;
    const INDEX numberOfLabels_;
    
  };

  DiscreteTomographyFactorCounting::DiscreteTomographyFactorCounting(INDEX numberOfLabels,INDEX a,INDEX b, INDEX c, INDEX d)
    : numberOfLabels_(numberOfLabels),a_(a),b_(b),c_(c),d_(d){
    assert(0 <= a_);
    assert(a_ <= b_);
    assert(b_ <= c_);
    assert(c_ <= d_);
    assert(numberOfLabels > 1);
  }
  
  template<typename REPAM_ARRAY>
  REAL DiscreteTomographyFactorCounting::LowerBound(const REPAM_ARRAY& repam) const{
    REAL m = std::numeric_limits<REAL>::max();
    for(INDEX i=0;i<repam.size();i++){
      for(INDEX j=0;j<repam[i].size();j++){
	m = std::min(m,repam[i][j]);
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

