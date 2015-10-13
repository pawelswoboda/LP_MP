#ifndef LP_MP_SIMPLEX_FACTOR_SIMD
#define LP_MP_SIMPLEX_FACTOR_SIMD

#include "LP_MP.h"
#include "factors_messages.hxx"
#include "const_array_types.h"
#include "factor_storage.hxx"

namespace LP_MP {

// expects Vc::Memory as reparametrization
// the polytope {x >= 0 : x_1 + ... + x_n = 1 }
// do zrobienia: possibly better derive from FactorStorage<STORE_FACTOR>
template<bool STORE_FACTOR>
class SimplexFactorSIMD
{
public:
   SimplexFactorSIMD(const std::vector<REAL>& cost, const std::vector<INDEX>& varCapacity, INDEX sum)
      : SimplexFactorSIMD(cost) {} // for compatibility with Multiplex
   SimplexFactorSIMD(const std::vector<REAL>& cost)
      : c_(cost)
   {
   }

   std::vector<REAL> GetPotential() {
      return c_;
   }

   void MaximizePotential() {};
   template<typename G>
   REAL LowerBound(const G& repamPot) const { 
      REAL_SIMD l(repamPot.vector(1));
      // do zrobienia: check if it automatically pads the rest correctly
      for(INDEX i=1; i<repamPot.vectorsCount(); ++i) {
         l = Vc::min(l, repamPot.vector(i));
      }
      REAL ls = l.min();
      return ls;
   }

   template<typename REPAM_ARRAY>
   REAL GetBreakpointCost(const REPAM_ARRAY& repamPot) const { 
      return LowerBound(repamPot); 
   }

   INDEX GetSum() const { return 1; }

   REAL operator[](const INDEX i) const { return c_[i]; }
   INDEX size() const { return c_.size(); }
private:
   FactorStorage<STORE_FACTOR> c_;
};


} // end namespace LP_MP

#endif // LP_MP_SIMPLEX_FACTOR_SIMD


