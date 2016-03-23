#ifndef LP_MP_MULTIPLEX_SUM_FACTOR_HXX
#define LP_MP_MULTIPLEX_SUM_FACTOR_HXX

#include "factors/multiplex_factor.hxx"
#include "messages/message_loops.hxx"

namespace LP_MP {

// implements the factor A_1*(mu_1 + ... + mu_n) = A_2*nu,
// where mu_1,...,mu_n are multiplexes as is nu, and A is a 0/1 matrix.
// templates LEFT_MULTIPLEX_FACTOR_TYPE and RIGHT_MULTIPLEX_FACTOR_TYPE should not hold cost explicitly and cost should be zero.
// for now we assume that LEFT_MULTIPLEX_FACTOR is actually a simplex
//
// we assume that this factor is joined via multiplex_marg_message with UnaryLoops to the mu-simplices and with MarginalSummationMessage to multplex nu.
// message operations then simulate the above ones but on modified costs.
//
// do zrobienia: rename left... into summand... and right... into sum...
// do zrobienia: possibly do not hold the sum explicitly and compute.
// Then the sum factor is a view on the sum of all summands
template<typename LEFT_MULTIPLEX_FACTOR_TYPE, typename RIGHT_MULTIPLEX_FACTOR_TYPE>
class MultiplexSumFactor {
public:
   using LeftMultiplexFactorType = LEFT_MULTIPLEX_FACTOR_TYPE;
   using RightMultiplexFactorType = RIGHT_MULTIPLEX_FACTOR_TYPE;

   // possibly incorporate also capacities and individual sums. For now assume them fixed.
   MultiplexSumFactor(const INDEX leftMultiplexFactorNo, const INDEX multiplexDim)
      : lm_(LeftMultiplexFactorType(std::vector<REAL>(multiplexDim,0.0))),
      rm_(RightMultiplexFactorType(std::vector<REAL>(multiplexDim,0.0), std::vector<INDEX>(multiplexDim,leftMultiplexFactorNo), leftMultiplexFactorNo)),
      leftMultiplexFactorNo_(leftMultiplexFactorNo)
   {}

   template<typename REPAM_ARRAY>
   const REAL LowerBound(const REPAM_ARRAY& repamPot) const {
      // do zrobienia: not correct yet, take into account equality constraints
      REAL lb = 0.0;
      // summands
      for(INDEX i=0; i<leftMultiplexFactorNo_; ++i) {
         VecSlice<const REPAM_ARRAY> repamSlice (repamPot, i*lm_.size(), (i+1)*lm_.size());
         VecSlice<const REPAM_ARRAY> sum (repamPot, leftMultiplexFactorNo_*lm_.size(), (1+leftMultiplexFactorNo_)*lm_.size());
         PlusExprVec<decltype(repamSlice), decltype(sum)> r(repamSlice,sum);
         lb += lm_.LowerBound(r);
         //lb += lm_.LowerBound(repamSlice);
      } 
      // the sum
      //VecSlice<const REPAM_ARRAY> repamSlice (repamPot, leftMultiplexFactorNo_*lm_.size(), (leftMultiplexFactorNo_+1)*lm_.size());
      //lb += rm_.LowerBound(repamSlice);
      return lb;
   }

   void MaximizePotential() {}

   const REAL operator[](const INDEX i) const { return 0.0; }
   const INDEX size() const { return (1+leftMultiplexFactorNo_)*rm_.size(); }
   const INDEX GetMultiplexDim() const { return lm_.size(); }
   const INDEX GetMultiplexNo() const { return leftMultiplexFactorNo_; }
   LeftMultiplexFactorType* GetSummandFactor() { return &lm_; }
   RightMultiplexFactorType* GetSumFactor() { return &rm_; }

private:
   // do zrobienia: holding them statically would also be possible for every fixed dimension multiplexDim. Also not holding them at all and only instantiating them whenever needed is also possible. For this more efficient initialization of constant multiplex factors is needed, though.
   LeftMultiplexFactorType lm_; // the summand
   RightMultiplexFactorType rm_; // the sum
   const INDEX leftMultiplexFactorNo_; // number of multiplex factors on the left
};

} // end namespace LP_MP

#endif //  LP_MP_MULTIPLEX_SUM_FACTOR_HXX
