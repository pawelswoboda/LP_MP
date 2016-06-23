#ifndef LP_MP_SIMPLEX_FACTOR_HXX
#define LP_MP_SIMPLEX_FACTOR_HXX

#include "LP_MP.h"
#include "factors_messages.hxx"
#include "const_array_types.h"

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
   static void MaximizePotentialAndComputePrimal(const REPAM_ARRAY& repam, typename PrimalSolutionStorage::Element primal)
   {
      for(INDEX i=0; i<repam.size(); ++i) { // ensure that primal has been initialized correctly to true
         assert(primal[i] == unknownState); // note: more general procedure is possible, but not implemented yet. In general, if already one state is primal, set all unknown states to false. Otherwise choose minimum over unknown states and set primal to it. use randomization?
      }
      // note: currently possibly also pairwise factors are called here, although this should not be made for SRMP style rounding
      INDEX min_element;
      REAL min_value = std::numeric_limits<REAL>::max();
      for(INDEX i=0; i<repam.size(); ++i) {
         primal[i] = false;
         if(min_value > repam[i]) {
            min_value = repam[i];
            min_element = i;
         }
      }
      assert(min_element < repam.size());
      primal[min_element] = true;
   };

   template<typename REPAM_ARRAY>
   static REAL LowerBound(const REPAM_ARRAY& repamPot) {
      REAL min_val = repamPot[0];
      for(INDEX i=1; i<repamPot.size(); ++i) {
         min_val = std::min(min_val, repamPot[i]);
      }
      return min_val;
   }

   template<typename REPAM_ARRAY>
   static REAL EvaluatePrimal(const REPAM_ARRAY& repam, const PrimalSolutionStorage::Element primal) 
   {
      REAL cost = std::numeric_limits<REAL>::max();
      INDEX noActive = 0;
      INDEX noUnknown = 0;
      INDEX unknownIndex;
      for(INDEX i=0; i<repam.size(); ++i) {
         if(primal[i] == true) { 
            ++noActive;
            cost = repam[i];
         }
         if(primal[i] == unknownState) {
            ++noUnknown;
            unknownIndex = i;
         }
      }
      if(noActive == 0 && noUnknown == 1) {
         // set 
         primal[unknownIndex] = true;
         return repam[unknownIndex];
      } else if(noActive == 1 && noUnknown == 0) {
         return cost;
      } else if(noActive == 1) {
         for(INDEX i=0; i<repam.size(); ++i) {
            if(primal[i] == unknownState) {
               primal[i] = false;
            }
         }
         return cost;
      } else {
         // most probably something is not working in primal propagation
         return std::numeric_limits<REAL>::infinity();
      }
   }

   void WritePrimal(const INDEX primal, std::ofstream& fs) const
   {
      fs << primal;
   }
};

} // end namespace LP_MP

#endif // LP_MP_SIMPLEX_FACTOR_HXX

