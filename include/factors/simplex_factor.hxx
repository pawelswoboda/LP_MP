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
   static void MaximizePotentialAndComputePrimal(const REPAM_ARRAY& repam, typename PrimalSolutionStorage::Element& primal)
   {
      //std::cout << "Address of first element of primal in simplex factor: " << &(*primal) << "\n";
      for(INDEX i=0; i<repam.size(); ++i) { // ensure that primal has been initialized correctly to true
         assert(primal[i] == true);
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
      for(INDEX i=0; i<repam.size(); ++i) {
         if(primal[i]) {
            ++noActive;
            cost = repam[i];
         }
      }
      if(noActive == 1) {
         return cost;
      } else {
         if(repam.size() < 100) {
            for(INDEX i=0; i<repam.size(); ++i) {
               std::cout << primal[i] << ",";
            }
            std::cout << "\n";
            assert(false); // this should not happen. If yes, primal propagation has done something wrong
         }
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

