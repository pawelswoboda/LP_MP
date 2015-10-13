#ifndef LP_MP_PAIRWISE_SIMPLEX_FACTOR_SIMD
#define LP_MP_PAIRWISE_SIMPLEX_FACTOR_SIMD

#include "LP_MP.h"
#include "factors_messages.hxx"
#include "const_array_types.h"
namespace LP_MP {


// expects Vc::Memory compatible as reparametrization
// the polytope {x >= 0 : x_11 + ... + x_nk = 1 }
// inherit from SimplexFactorSIMD and add GetDim function
class PairwiseSimplexFactorSIMD
{
public:
   //PairwiseSimplexFactorSIMD(const std::vector<REAL>& cost, const std::vector<INDEX>& varCapacity, INDEX sum)
   //   : PairwiseSimplexFactorSIMD(cost) { throw std::runtime_error("not dimensions specified"); } // for compatibility with Multiplex
   PairwiseSimplexFactorSIMD(const INDEX dim1, const INDEX dim2, const std::vector<REAL>& cost, const std::vector<INDEX>& varCapacity, INDEX sum)
      : dim1_(dim1),
      dim2_(dim2),
      c_(cost)
   { }

   std::vector<REAL> GetPotential() {
      return c_;
   }

   void MaximizePotential() {};
   template<typename G>
   REAL LowerBound(const G& repamPot) const { 
      REAL_SIMD l(repamPot.vector(0));
      // do zrobienia: check if it automatically pads the rest correctly
      for(INDEX i=1; i<repamPot.vectorsCount(); ++i) {
         l = Vc::min(l, repamPot.vector(i));
      }
      assert(l.min() < 1000000.0 && l.min() > -10000000.0);
      return l.min();
   }

   template<typename REPAM_ARRAY>
   REAL GetBreakpointCost(const REPAM_ARRAY& repamPot) const { 
      return LowerBound(repamPot); 
   }

   INDEX GetSum() const { return 1; }

   INDEX GetDim(const INDEX i) const { 
      if(i==0) return dim1_;
      if(i==1) return dim2_;
      throw std::runtime_error("two dimensions present");
   }

   REAL operator[](const INDEX i) const { return c_[i]; }
   INDEX size() const { return c_.size(); }
private:
   const INDEX dim1_, dim2_;
   const std::vector<REAL> c_;
};


} // end namespace LP_MP

#endif // LP_MP_PAIRWISE_SIMPLEX_FACTOR_SIMD



