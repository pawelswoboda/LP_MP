#ifndef LP_MP_SIMPLEX_FACTOR_HXX
#define LP_MP_SIMPLEX_FACTOR_HXX

#include "LP_MP.h"
#include "factors_messages.hxx"
#include "const_array_types.h"

// rename this class to SimplexFactor
// Investigate contingency tables (Knuth) for more general tabular structures.

namespace LP_MP {

// the polytope {x >= 0 : x_1 + ... + x_n = 1 }
// dual problem:
//      max_{z} z 
//      s.t.    z <= repamPot[i]

// this class does nothing and is used, when reparametrization is held explicitly by RepamStorage. Specialized classes for holding actual information can be supplied, though.
struct SimplexEmptyVector {
template<typename ARRAY> SimplexEmptyVector(const ARRAY& a) : size_(a.size()) {}
const INDEX size() const { return size_; }
private:
const INDEX size_;
};

template<typename COST_ARRAY = SimplexEmptyVector>
class SimplexFactor : public COST_ARRAY
{
public:
   SimplexFactor(const std::vector<REAL>& cost)
      : COST_ARRAY(cost)
   {}
   //SimplexFactor() {}

   template<typename REPAM_ARRAY>
   static void MaximizePotentialAndComputePrimal(const REPAM_ARRAY& repam, typename PrimalSolutionStorage::Element primal)
   {
      // note: currently possibly also pairwise factors are called here, although this should not be made for SRMP style rounding
      INDEX min_element;
      REAL min_value = std::numeric_limits<REAL>::infinity();
      for(INDEX i=0; i<repam.size(); ++i) {
         primal[i] = false;
         if(min_value >= repam[i]) {
            min_value = repam[i];
            min_element = i;
         }
      }
      assert(min_element < repam.size());
      primal[min_element] = true;
      //std::cout << ";    " << min_element;
      //std::cout << "\n";
   };

   // set to false all elements that are larger by epsilon than the minimal element
   template<typename REPAM_ARRAY>
   static void NarrowPrimal(const REPAM_ARRAY& repam, const REAL epsilon, typename PrimalSolutionStorage::Element primal)
   {
      REAL min_val = std::numeric_limits<REAL>::infinity();
      for(INDEX i=0; i<repam.size(); ++i) {
         if(primal[i] != false) {
            std::min(repam[i], min_val);
         }
      }
      for(INDEX i=0; i<repam.size(); ++i) {
         if(repam[i] >= min_val + epsilon) {
            primal[i] = false;
         }
      }
   }

   template<typename REPAM_ARRAY>
   static REAL LowerBound(const REPAM_ARRAY& repamPot) {
      REAL min_val = repamPot[0];
      for(INDEX i=1; i<repamPot.size(); ++i) {
         min_val = std::min(min_val, repamPot[i]);
      }
      return min_val;
   }

   // if there is exactly one unknownIndex (rest false), make it true
   // if there is already one true index, make all other false
   void PropagatePrimal(PrimalSolutionStorage::Element primal)
   {
      INDEX noTrue = 0;
      INDEX noUnknown = 0;
      INDEX unknownIndex;
      for(INDEX i=0; i<this->size(); ++i) {
         if(primal[i] == true) { 
            ++noTrue;
         }
         if(primal[i] == unknownState) {
            ++noUnknown;
            unknownIndex = i;
         }
      }
      if(noTrue == 0 && noUnknown == 1) {
         primal[unknownIndex] = true;
      } else if(noTrue == 1 && noUnknown > 0) {
         for(INDEX i=0; i<this->size(); ++i) {
            if(primal[i] == unknownState) {
               primal[i] = false;
            }
         }
      }
   }

   template<typename REPAM_ARRAY>
   static REAL EvaluatePrimal(const REPAM_ARRAY& repam, const PrimalSolutionStorage::Element primal) 
   {
      REAL cost = std::numeric_limits<REAL>::max();
      INDEX primalSum = 0;
      for(INDEX i=0; i<repam.size(); ++i) {
         if(primal[i] == true) { 
            cost = repam[i];
         }
         primalSum += primal[i];
      }
      if(primalSum == 1) {
         return cost;
      } else {
         //std::cout << "primal not inferred correctly: " << primalSum << "\n";
         return std::numeric_limits<REAL>::infinity();
      }
   }

  INDEX GetNumberOfAuxVariables() const { return 0; }
   
  void CreateConstraints(LpInterfaceAdapter* lp) const { 
    LinExpr lhs = lp->CreateLinExpr();
    for(INDEX i=0;i<lp->GetFactorSize();i++){
      lhs += lp->GetVariable(i);
    }
    LinExpr rhs = lp->CreateLinExpr();
    rhs += 1;
    lp->addLinearEquality(lhs,rhs);
  }
  
   void WritePrimal(const INDEX primal, std::ofstream& fs) const
   {
      fs << primal;
   }
};

} // end namespace LP_MP

#endif // LP_MP_SIMPLEX_FACTOR_HXX

