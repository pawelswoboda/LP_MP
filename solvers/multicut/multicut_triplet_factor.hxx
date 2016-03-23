#ifndef LP_MP_MULTICUT_TRIPLET_FACTOR_HXX
#define LP_MP_MULTICUT_TRIPLET_FACTOR_HXX

#include "LP_MP.h"

namespace LP_MP {

// possibly better inherit from simplex factor with two labels
class MulticutTripletFactor 
{
public:
   struct LabelingType : public std::array<bool,3> // try INDEX:1;
   {
      LabelingType(const bool b = false) : std::array<bool,3>({b,b,b}) {}
   };
   using PrimalType = LabelingType;
                         
   using IndicesType = std::array<INDEX,3>; // try INDEX:2;

   MulticutTripletFactor() {}; 
   template<typename REPAM_ARRAY>
   void MaximizePotential(const REPAM_ARRAY& repam) {};
   // do zrobienia: remove this function
   /*
   template<typename REPAM_ARRAY>
   std::pair<LabelingType,REAL> MaximizePotentialAndComputePrimal(const REPAM_ARRAY& repam) 
   {
   }
   */
   template<typename REPAM_ARRAY>
   IndicesType SortIndices(const REPAM_ARRAY& repamPot) const
   {
      //std::cout << "kwas: do zrobienia sorting\n";
      IndicesType i{0,1,2};
      std::sort(i.begin(), i.end(), [&repamPot](const INDEX i1, const INDEX i2) { return repamPot[i1] < repamPot[i2]; });
      return i;
   }
   template<typename REPAM_ARRAY>
   LabelingType ComputeOptimalLabeling(const REPAM_ARRAY& repamPot) const
   {
      assert(repamPot.size() == 3);
      // sort repamPot. Question: How to do this fastest?
      std::array<INDEX,3> sortIndices = SortIndices(repamPot);
      if(repamPot[sortIndices[0]] + repamPot[sortIndices[1]] > 0.0) {
         return LabelingType(false);//{false,false,false};
      } else if(repamPot[sortIndices[2]] > 0) {
         LabelingType labeling;
         labeling[sortIndices[0]] = true;
         labeling[sortIndices[1]] = true;
         labeling[sortIndices[2]] = false;
         return labeling;
      } else {
         return LabelingType(true);//{true,true,true};
      }
   }
   template<typename REPAM_ARRAY>
   REAL LowerBound(const REPAM_ARRAY& repamPot) const {
      assert(repamPot.size() == 3);
      IndicesType sortIndices = SortIndices(repamPot);
      if(repamPot[sortIndices[0]] + repamPot[sortIndices[1]] > 0.0) {
         return 0.0;
      } else if(repamPot[sortIndices[2]] > 0) {
         return repamPot[sortIndices[0]] + repamPot[sortIndices[1]];
      } else {
         return repamPot[0] + repamPot[1] + repamPot[2];
      }
   }

   const REAL operator[](const INDEX i) const { assert(i<3); return 0.0; }
   const INDEX size() const { return 3; }

   template<typename REPAM_ARRAY>
   REAL EvaluatePrimal(const REPAM_ARRAY& repam, const LabelingType primal) const
   {
      assert(repam.size() == 3);
      const INDEX sum = primal[0] + primal[1] + primal[2];
      if(sum == 1) return 1e11;
      else return repam[0]*primal[0] + repam[1]*primal[1] + repam[2]*primal[2];
   }
   void WritePrimal(const LabelingType, std::ofstream& fs) const
   {
   }

private:
};

} // end namespace LP_MP

#endif // LP_MP_MULTICUT_TRIPLET_FACTOR



