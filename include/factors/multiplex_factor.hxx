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
template<typename COST_ARRAY = std::vector<REAL>, typename VAR_CAPACITY_ARRAY = std::vector<INDEX>, typename SUM = INDEX>
class MultiplexFactor
{
public:
   MultiplexFactor(const std::vector<REAL>& cost, const std::vector<INDEX>& varCapacity, INDEX sum)
      : 
      c_(cost),
      varCapacity_(varCapacity),
      sum_(sum)
   {
      assert(sum_ == 1); 
   }


   template<typename COST_ARRAY_TMP = COST_ARRAY>
   MultiplexFactor(const std::vector<REAL>& cost) // fallback to simplex
      : c_(cost),
      varCapacity_(std::vector<INDEX>(cost.size(),1)),
      sum_(1)
   {}

   template<typename REPAM_ARRAY>
   void MaximizePotential(const REPAM_ARRAY& repam) {};
   template<typename REPAM_ARRAY>
   std::pair<INDEX,REAL> MaximizePotentialAndComputePrimal(const REPAM_ARRAY& repam) 
   {
      // note: currently possibly also pairwise factors are called here, although this should not be made for SRMP style rounding
      INDEX min_element = 0;
      REAL min_value = repam[0];
      for(INDEX i=1; i<repam.size(); ++i) {
         if(min_value > repam[i]) {
            min_value = repam[i];
            min_element = i;
         }
      }
      return std::make_pair(min_element, min_value);
   };

   template<typename REPAM_ARRAY>
   REAL LowerBound(const REPAM_ARRAY& repamPot) const {
      REAL min_val = repamPot[0];
      for(INDEX i=1; i<repamPot.size(); ++i) {
         min_val = std::min(min_val, repamPot[i]);
      }
      return min_val;
   }

   // only provide those if cost needs to be stored
   const REAL operator[](const INDEX i) const { return c_[i]; }
   const INDEX size() const { return c_.size(); }

   template<typename REPAM_ARRAY>
   REAL EvaluatePrimal(const REPAM_ARRAY& repam, const INDEX primal) const
   {
      assert(primal<repam.size());
      return repam[primal];
   }
   void WritePrimal(const INDEX primal, std::ofstream& fs) const
   {
      fs << primal;
   }

private:
   // do zrobienia: privately inherit COST_ARRAY c_ to enable empty base class optimization
   const COST_ARRAY c_; // const not  doable in array constructor
   const VAR_CAPACITY_ARRAY varCapacity_; // maximum of each variable
   const SUM sum_; // sum of variables must be equal to
};

} // end namespace LP_MP

#endif // LP_MP_SIMPLEX_FACTOR_HXX

