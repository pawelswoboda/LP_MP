#ifndef LP_MP_MULTICUT_UNARY_FACTOR_HXX
#define LP_MP_MULTICUT_UNARY_FACTOR_HXX

#include "LP_MP.h"
#include <random>

namespace LP_MP {

// possibly better inherit from simplex factor with two labels
class MulticutUnaryFactor 
{
public:
   using PrimalType = bool;
   MulticutUnaryFactor(const double cost) {};
   template<typename REPAM_ARRAY>
   void MaximizePotential(const REPAM_ARRAY& repam) {};
   template<typename REPAM_ARRAY>
   void MaximizePotentialAndComputePrimal(const REPAM_ARRAY& repam, typename PrimalSolutionStorage::Element primal)
   {
      assert(repam.size() == 1);
      /*
      if(repam[0] == 0.0) { // round solution, flip coin. it does not seem to help, but I will double check this
         if(r() == 1) {
            return std::make_pair(true,0.0);
         } else {
            return std::make_pair(false,0.0);
         }
      }
      */
      if(repam[0] <= 0) { primal[0] = true; }
      else { primal[1] = false; }
   }
   template<typename REPAM_ARRAY>
   REAL LowerBound(const REPAM_ARRAY& repamPot) const {
      assert(repamPot.size() == 1);
      return std::min(repamPot[0],0.0);
   }

   const INDEX size() const { return 1; }

   template<typename REPAM_ARRAY>
   REAL EvaluatePrimal(const REPAM_ARRAY& repam, const PrimalSolutionStorage::Element primal) const
   {
      assert(repam.size() == 1);
      return primal[0]*repam[0];
   }
   void WritePrimal(const PrimalSolutionStorage::Element primal, std::ofstream& fs) const
   {
      assert(false);
      //fs << primal;
   }

private:

   static std::uniform_int_distribution<>::param_type p;
   static decltype(std::bind(std::uniform_int_distribution<>{p}, std::default_random_engine{})) r;
};
// very ugly, make nicer
std::uniform_int_distribution<>::param_type MulticutUnaryFactor::p = decltype(MulticutUnaryFactor::p){0,1};
decltype(std::bind(std::uniform_int_distribution<>{MulticutUnaryFactor::p}, std::default_random_engine{})) MulticutUnaryFactor::r = std::bind(std::uniform_int_distribution<>{MulticutUnaryFactor::p}, std::default_random_engine{});

} // end namespace LP_MP

#endif // LP_MP_MULTICUT_UNARY_FACTOR



