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

    // --------------------
    // custom public methods 
    
    void setRHS(INDEX b){ root_ = true; rhs_ = b; }
    void setNOISE(std::function<REAL (REAL)> n){ weights_ = n;  }
    
  private:

    const INDEX a_,b_,c_,d_;
    const INDEX numberOfLabels_;
    INDEX upSize_,leftSize_,rightSize_,regSize_;
    INDEX rhs_;
    bool root_ = false;
    std::function<REAL (REAL)> weights_ = NULL;
    
    
    bool checkConstraint(INDEX,INDEX,INDEX) const; // z_{a:b} + z_{c:d} = z_{a:d}
    INDEX getLabelPair(INDEX,INDEX) const;
    REAL getNoise(INDEX,INDEX,INDEX) const;
    
  };

  DiscreteTomographyFactorCounting::DiscreteTomographyFactorCounting(INDEX numberOfLabels,INDEX a,INDEX b, INDEX c, INDEX d)
    : numberOfLabels_(numberOfLabels),a_(a),b_(b),c_(c),d_(d){
    assert(0 <= a_); assert(a_ <= b_); assert(b_ <= c_); assert(c_ <= d_);
    assert(numberOfLabels > 1);
    
    upSize_ = pow(numberOfLabels_,3)*(d_-a_+1);
    leftSize_ = pow(numberOfLabels_,3)*(b_-a_+1);
    rightSize_ = pow(numberOfLabels_,3)*(c_-b_+1);
    regSize_ = pow(numberOfLabels_,2);
  }

  template<typename REPAM_ARRAY>
  REAL DiscreteTomographyFactorCounting::LowerBound(const REPAM_ARRAY& repam) const{
    assert(repam.size() == (upSize_ + leftSize_ + rightSize_ + regSize_));
    REAL m = std::numeric_limits<REAL>::max();
    for( INDEX u=0;u < upSize_;u++ ){
      for( INDEX l=0;l < leftSize_;l++ ){
	for( INDEX r=0;r < rightSize_;r++ ){
	  if( !checkConstraint(u,l,r) && !root_ ){ continue; }
	  INDEX noise = 0;
	  if( root_ ){
	    noise = getNoise(u,l,r);
	  }
	  m = std::min(m,noise
		       + repam[u]
		       + repam[upSize_ + l] 
		       + repam[upSize_ + leftSize_ + r]
		       + repam[upSize_ + leftSize_ + rightSize_ + getLabelPair(l,r)] );
	}
      }
    }
    return m;
  }

  REAL DiscreteTomographyFactorCounting::getNoise(INDEX u,INDEX l,INDEX r) const{
    if( weights_ == NULL ){ return 0; }
    else{
      auto exF = [&](INDEX i){
	i = (i-(i % numberOfLabels_))/numberOfLabels_;
	return (i-(i % numberOfLabels_))/numberOfLabels_;
      };
      INDEX uz = exF(u );
      INDEX lz = exF(l);
      INDEX rz = exF(r);
      REAL n = std::abs((REAL)uz - (REAL)lz - (REAL)rz);
      return weights_(n);
    }
  }
  
  INDEX DiscreteTomographyFactorCounting::getLabelPair(INDEX l,INDEX r) const{
    INDEX lb = l % numberOfLabels_;
    INDEX ra = ((r - (r % numberOfLabels_))/numberOfLabels_) % numberOfLabels_;
    return lb + ra;
  }
  
  bool DiscreteTomographyFactorCounting::checkConstraint(INDEX up,INDEX left,INDEX right) const{
    auto exF = [&](INDEX i){
      i = (i-(i % numberOfLabels_))/numberOfLabels_;
      return (i-(i % numberOfLabels_))/numberOfLabels_;
    };
    INDEX uz = exF(up);
    INDEX lz = exF(left);
    INDEX rz = exF(right);
    return (lz + rz == uz);
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

