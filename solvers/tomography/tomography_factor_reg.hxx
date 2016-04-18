#ifndef LP_MP_TOMOGRAPHY_REG_FACTOR_HXX
#define LP_MP_TOMOGRAPHY_REG_FACTOR_HXX

#include "LP_MP.h"

namespace LP_MP{

  class TomographyRegFactor{

  public:
    
    TomographyRegFactor(); //--required

    template<typename REPAM_ARRAY>
    void MaximizePotential(const REPAM_ARRAY& repam); //--required

    template<typename REPAM_ARRAY>
    INDEX ComputeOptimalLabeling(const REPAM_ARRAY& repam) const; //--required

    template<typename REPAM_ARRAY>
    REAL LowerBound(const REPAM_ARRAY& repam) const; //--required

    template<typename REPAM_ARRAY>
    REAL EvaluatePrimal(const REPAM_ARRAY& repam, const PrimalSolutionStorage::Element primal) const; //--required

    void WritePrimal(const PrimalSolutionStorage::Element, std::ofstream& fs) const; //--required
    
  private:
    
  }

} // end namespace LP_MP

#endif // LP_MP_TOMOGRAPHY_REG_FACTOR_HXX

