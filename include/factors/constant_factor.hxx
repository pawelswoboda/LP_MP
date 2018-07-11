#ifndef LP_MP_CONSTANT_FACTOR_HXX
#define LP_MP_CONSTANT_FACTOR_HXX

#include "LP_MP.h"
#include "lp_interface/lp_interface.h"
#include <random>

namespace LP_MP {

class ConstantFactor {
public:
   ConstantFactor(const REAL offset = 0) : offset_(offset) {}

   constexpr static INDEX size() { return 0; }

   REAL LowerBound() const { return offset_; }

   REAL EvaluatePrimal() const
   {
      return offset_;
   }

   void AddToOffset(const REAL delta) { offset_ += delta; }

   void init_primal() {}
   template<typename ARCHIVE> void serialize_dual(ARCHIVE& ar) { ar(offset_); }
   template<typename ARCHIVE> void serialize_primal(ARCHIVE& ar) {}
private:
   REAL offset_;
};

} // end namespace LP_MP

#endif // LP_MP_CONSTANT_FACTOR_HXX

