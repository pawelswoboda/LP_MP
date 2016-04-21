#ifndef LP_MP_SIMPLEX_FACTOR_HXX
#define LP_MP_SIMPLEX_FACTOR_HXX

#include "LP_MP.h"
#include "factors_messages.hxx"
#include "const_array_types.h"
#include <tuple>

// rename this class to SimplexFactor
// Investigate contingency tables (Knuth) for more general tabular structures.

namespace LP_MP {

// the polytope {x >= 0 : x_1 + ... + x_n = 1 }
// dual problem:
//      max_{z} z 
//      s.t.    z <= repamPot[i]

// this class does nothing and is used, when reparametrization is held explicitly by RepamStorage. Specialized classes for holding actual information can be supplied, though.
struct SimplexEmptyVector {
template<typename ARRAY> SimplexEmptyVector(const ARRAY& a) {}
};
template<typename COST_ARRAY = SimplexEmptyVector>
class SimplexFactor : public COST_ARRAY
{
public:
   SimplexFactor(const std::vector<REAL>& cost)
      : COST_ARRAY(cost)
   {}
   SimplexFactor() {}

   template<typename REPAM_ARRAY>
   void MaximizePotentialAndComputePrimal(const REPAM_ARRAY& repam, typename PrimalSolutionStorage::Element primal)
   {
      // note: currently possibly also pairwise factors are called here, although this should not be made for SRMP style rounding
      INDEX min_element = 0;
      REAL min_value = repam[0];
      for(INDEX i=1; i<repam.size(); ++i) {
         primal[i] = false;
         if(min_value > repam[i]) {
            min_value = repam[i];
            min_element = i;
         }
      }
      primal[min_element] = true;
   };

   template<typename REPAM_ARRAY>
   REAL LowerBound(const REPAM_ARRAY& repamPot) const {
      REAL min_val = repamPot[0];
      for(INDEX i=1; i<repamPot.size(); ++i) {
         min_val = std::min(min_val, repamPot[i]);
      }
      return min_val;
   }

   template<typename REPAM_ARRAY>
   REAL EvaluatePrimal(const REPAM_ARRAY& repam, const PrimalSolutionStorage::Element primal) const
   {
      REAL cost = std::numeric_limits<REAL>::max();
      INDEX noActive = 0;
      for(INDEX i=0; i<repam.size(); ++i) {
         if(primal[i]) {
            ++noActive;
            cost = repam[i];
         }
      }
      if(noActive == 1) {
         return cost;
      } else {
         return std::numeric_limits<REAL>::max();
      }
   }

   void WritePrimal(const INDEX primal, std::ofstream& fs) const
   {
      fs << primal;
   }
};

} // end namespace LP_MP

#endif // LP_MP_SIMPLEX_FACTOR_HXX

