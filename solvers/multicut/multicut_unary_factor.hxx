#ifndef LP_MP_MULTICUT_UNARY_FACTOR_HXX
#define LP_MP_MULTICUT_UNARY_FACTOR_HXX

#include "LP_MP.h"
#include "lp_interface/lp_interface.h"
#include <random>

namespace LP_MP {

class MulticutUnaryFactor 
{
public:
   MulticutUnaryFactor(const double cost) : pot_(cost) {};
   void MaximizePotentialAndComputePrimal(typename PrimalSolutionStorage::Element primal)
   {
      if(pot_ <= 0) { 
         primal[0] = true; 
      } else { 
         primal[0] = false; 
      }
      /*
      if(std::abs(repam[0]) <= eps) { // round solution, flip coin. it does not seem to help, but I will double check this
         if(r() == 1) {
            primal[0] = true;
         } else {
            primal[0] = false;
         }
      }
      */
   }
   REAL LowerBound() const {
      return std::min(pot_,0.0);
   }

   constexpr static INDEX size() { return 1; }

   REAL EvaluatePrimal() const
   {
      return primal_*pot_;
   }

   operator REAL() const { return pot_; }
   operator REAL&() { return pot_; }
   //void CreateConstraints(LpInterfaceAdapter* lp) const {} // we do not have to do anything

   void init_primal() {}
   template<typename ARCHIVE> void serialize_dual(ARCHIVE& ar) { ar(pot_); }
   template<typename ARCHIVE> void serialize_primal(ARCHIVE& ar) { ar( primal_ ); }

private:
   REAL pot_;
   bool primal_;

   static std::uniform_int_distribution<>::param_type p;
   static decltype(std::bind(std::uniform_int_distribution<>{p}, std::default_random_engine{})) r;
};
// very ugly, make nicer
std::uniform_int_distribution<>::param_type MulticutUnaryFactor::p = decltype(MulticutUnaryFactor::p){0,1};
decltype(std::bind(std::uniform_int_distribution<>{MulticutUnaryFactor::p}, std::default_random_engine{})) MulticutUnaryFactor::r = std::bind(std::uniform_int_distribution<>{MulticutUnaryFactor::p}, std::default_random_engine{});

} // end namespace LP_MP

#endif // LP_MP_MULTICUT_UNARY_FACTOR



