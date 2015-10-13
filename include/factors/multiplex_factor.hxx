#ifndef LP_MP_MULTIPLEX_FACTOR
#define LP_MP_MULTIPLEX_FACTOR

#include "LP_MP.h"
#include "factors_messages.hxx"
#include "const_array_types.h"
//#include "multiplex_marg_message.hxx"

namespace LP_MP {


// the polytope {x >= 0 : x_i <= varCapacity_[i], x_1 + ... + x_n = sum_ }
// dual problem:
//      max_{z,y} sum_*z + varCapacity_[0]*y[0] + ... + varCapacity_[n-1]*y[n-1]
//      s.t.      z + y[i] <= \theta^{\phi}(i)
//                    y[i] <= 0
// template arguments give cost type, varCapacity type and sum type
template<typename COST_ARRAY = std::vector<REAL>, typename VAR_CAPACITY_ARRAY = std::vector<INDEX>, typename SUM = INDEX>
//template<typename COST_ARRAY = Vc::memory<REAL_SIMD>, typename VAR_CAPACITY_ARRAY = Vc::Memory<INDEX_SIMD>, typename SUM = INDEX>
class MultiplexFactor
{
public:
   //MultiplexFactor(const COST& cost, const VAR_CAPACITY& varCapacity, SUM sum)
   MultiplexFactor(const std::vector<REAL>& cost, const std::vector<INDEX>& varCapacity, INDEX sum)
      : //FactorBase(cost),
      c_(cost),
      varCapacity_(varCapacity),
      sum_(sum)
   {
      assert(sum_ == 1); // for now, later >= 1
   }

   std::vector<REAL> GetPotential() { return c_; }
   void MaximizePotential() {};
   template<typename G>
   REAL LowerBound(const G& repamPot) const { return LowerBoundImpl(repamPot, varCapacity_, sum_); }

   template<typename REPAM_ARRAY>
   REAL GetBreakpointCost(const REPAM_ARRAY& repamPot) const { return GetBreakpointCostImpl(repamPot, varCapacity_, sum_); }

   INDEX GetSum() const { return sum_; }

   REAL operator[](const INDEX i) const { return c_[i]; }
   INDEX size() const { return c_.size(); }
private:
   const COST_ARRAY c_;
   const VAR_CAPACITY_ARRAY varCapacity_; // maximum of each variable
   const SUM sum_; // sum of variables must be equal to
};

// do zrobienia: write specialization for sum_ = 1 and sum_ = 2;

template<typename REPAM_ARRAY, typename VAR_CAPACITY_ARRAY, typename SUM>
REAL LowerBoundImpl(const REPAM_ARRAY& repamPot, const VAR_CAPACITY_ARRAY& varCapacity, const SUM sum)
{
   std::cout << "no Specialization\n"; 
   std::vector<INDEX> indices(repamPot.size()); // do zrobienia: preallocate memory for that, or find way to push on stack
   for(INDEX i=0; i<indices.size(); ++i) indices[i] = i;
   // do not sort all, but only the first sum_ ones
   std::partial_sort(indices.begin(), indices.begin() + std::min(sum, indices.size()), indices.end(), [&](const INDEX i, const INDEX j) { return repamPot[i] < repamPot[j]; });
   //std::sort(indices.begin(), indices.end(), [&](const INDEX i, const INDEX j) { return pot[i] < pot[j]; });
   REAL dualCost = 0.0;
   INDEX remainingCapacity = sum;
   for(INDEX i=0; i<indices.size(); ++i) {
      dualCost += repamPot[indices[i]] * std::min(varCapacity[indices[i]], remainingCapacity);
      remainingCapacity -= varCapacity[indices[i]];
      if(remainingCapacity <= 0) return dualCost;
   }
   throw std::runtime_error("variable capacities smaller than sum");
}

template<typename REPAM_ARRAY, typename VAR_CAPACITY_ARRAY>
REAL LowerBoundImpl(const REPAM_ARRAY& repamPot, const VAR_CAPACITY_ARRAY& varCapacity, const const_one sum)
{
   REAL min_val = std::numeric_limits<REAL>::max();
   for(INDEX i=0; i<repamPot.size(); ++i) {
      min_val = std::min(min_val, repamPot[i]);
   }
   return min_val;
}

/*
REAL LowerBoundImpl(const std::valarray<REAL>& repamPot, const const_ones_array varCapacity, const const_one sum)
{
   std::cout << "kwas\n";
   return repamPot.min();
}
*/

// get reparametrized cost for which the breakpoint is active in the sense that we can reduce all cost higher than breakpoint to breakpoint and dual cost does not change
// do zrobienia: preallocate indices
template<typename REPAM_ARRAY, typename VAR_CAPACITY_ARRAY, typename SUM>
REAL GetBreakpoint(const REPAM_ARRAY& repamPot, const VAR_CAPACITY_ARRAY& varCapacity, const SUM sum)
{
   std::vector<INDEX> indices(repamPot.size());
   for(INDEX i=0; i<indices.size(); ++i) indices[i] = i;
   // do zrobienia: possibly do not sort all, but only the first sum_ ones
   std::partial_sort(indices.begin(), indices.begin() + std::min(sum, indices.size()), indices.end(), [&](const INDEX i, const INDEX j) { return repamPot[i] < repamPot[j]; });
   //std::sort(indices.begin(), indices.end(), [&](const INDEX i, const INDEX j) { return repamPot[i] < repamPot[j]; });
   INDEX remainingCapacity = sum;
   for(INDEX i=0; i<indices.size(); ++i) {
      remainingCapacity -= varCapacity[indices[i]];
      if(remainingCapacity <= 0) return repamPot[indices[i]];;
   }
   throw std::runtime_error("variable capacities smaller than sum");
}

template<typename REPAM_ARRAY, typename VAR_CAPACITY_ARRAY>
REAL GetBreakpoint(const REPAM_ARRAY& repamPot, const VAR_CAPACITY_ARRAY& varCapacity, const const_one sum)
{
   REAL min_val = std::numeric_limits<REAL>::max();
   for(INDEX i=0; i<repamPot.size(); ++i) {
      min_val = std::min(min_val, repamPot[i]);
   }
   return min_val;

}

} // end namespace LP_MP

#endif // LP_MP_MULTIPLEX_FACTOR

