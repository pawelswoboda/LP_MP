#ifndef LP_MP_EMPTY_FACTOR_HXX
#define LP_MP_EMPTY_FACTOR_HXX

#include "instances.inc"

namespace LP_MP {

class EmptyFactor
{
public:
   EmptyFactor(const std::vector<REAL>& cost) {}

   // do zrobienia: is this really needed?
   std::vector<REAL> GetPotential() {
      return std::vector<REAL>(0);
   }

   void MaximizePotential() {}
   template<typename G>
   REAL LowerBound(const G& repamPot) const { 
      return 0.0;
   }

   REAL operator[](const INDEX i) const { return 0; }
   INDEX size() const { return 0; }
};


} // end namespace LP_MP

#endif // LP_MP_EMPTY_FACTOR_HXX


