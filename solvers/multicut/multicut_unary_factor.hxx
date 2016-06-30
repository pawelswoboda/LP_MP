#ifndef LP_MP_MULTICUT_UNARY_FACTOR_HXX
#define LP_MP_MULTICUT_UNARY_FACTOR_HXX

#include "LP_MP.h"
#include <random>

namespace LP_MP {

// possibly better inherit from simplex factor with two labels
class MulticutUnaryFactor 
{
public:
   MulticutUnaryFactor(const double cost) {};
   template<typename REPAM_ARRAY>
   static void MaximizePotentialAndComputePrimal(const REPAM_ARRAY& repam, typename PrimalSolutionStorage::Element primal)
   {
      assert(repam.size() == 1);
      if(repam[0] <= 0) { 
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
   template<typename REPAM_ARRAY>
   static REAL LowerBound(const REPAM_ARRAY& repamPot) {
      assert(repamPot.size() == 1);
      return std::min(repamPot[0],0.0);
   }

   constexpr static INDEX size() { return 1; }

   template<typename REPAM_ARRAY>
   static REAL EvaluatePrimal(const REPAM_ARRAY& repam, const PrimalSolutionStorage::Element primal)
   {
      assert(repam.size() == 1);
      assert(primal[0] == false || primal[0] == true);
      return primal[0]*repam[0];
   }
   void WritePrimal(const PrimalSolutionStorage::Element primal, std::ofstream& fs) const
   {
      //fs << primal[0];
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



