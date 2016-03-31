#ifndef LP_MP_MULTICUT_TRIPLET_FACTOR_HXX
#define LP_MP_MULTICUT_TRIPLET_FACTOR_HXX

#include "LP_MP.h"
#include "cpp-sort/sort.h"
#include "cpp-sort/sorters.h"
#include "cpp-sort/fixed_sorters.h"

namespace LP_MP {

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

      // default
      //std::sort(i.begin(), i.end(), [&repamPot](const INDEX i1, const INDEX i2) { return repamPot[i1] < repamPot[i2]; });
      
      // possibly faster sorting network.
      // seems marginally faster
      using network_sorter = cppsort::sorting_network_sorter<3>; //cppsort::small_array_adapter< cppsort::sorting_network_sorter<3> >;
      cppsort::sort(i, network_sorter{}, [&repamPot](const INDEX i1, const INDEX i2) { return repamPot[i1] < repamPot[i2]; }); 
      // possibly even faster: define struct {REAL val, INDEX key}; and sort by val and then use this directly in subsequent optimization inside message. Or somehow embed in low order bits the indices and sort directly
      // possibly can be sorted with SIMD. Vc has sorting routines

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
      std::array<REAL,3> sortRepam{repamPot[0], repamPot[1], repamPot[2]};
      using network_sorter = cppsort::sorting_network_sorter<3>;
      cppsort::sort(sortRepam, network_sorter{}); 
      const REAL lb =  std::min(sortRepam[0] + sortRepam[1],0.0) + std::min(sortRepam[2],0.0);
      return lb;
      //return lb;
      /*
      std::array<REAL,3> sortRepam{repamPot[0], repamPot[1], repamPot[2]};
      using network_sorter = cppsort::sorting_network_sorter<3>;
      cppsort::sort(sortRepam, network_sorter{}); 
      if(sortRepam[0] + sortRepam[1] > 0) {
         return 0;
      } else if(sortRepam[2] > 0) {
         return sortRepam[0] + sortRepam[1];
      } else {
         return sortRepam[0] + sortRepam[1] + sortRepam[2];
      }
      */

      // strangely, this seems to be marginally faster than the approach above
      /*
      IndicesType sortIndices = SortIndices(repamPot);
      if(repamPot[sortIndices[0]] + repamPot[sortIndices[1]] > 0.0) {
         assert(lb == 0.0);
         return 0.0;
      } else if(repamPot[sortIndices[2]] > 0) {
         assert(std::abs(lb - (repamPot[sortIndices[0]] + repamPot[sortIndices[1]])) < eps);
         return repamPot[sortIndices[0]] + repamPot[sortIndices[1]];
      } else {
         assert(std::abs(lb - (repamPot[0] + repamPot[1] + repamPot[2])) <= eps);
         return repamPot[0] + repamPot[1] + repamPot[2];
      }
      */
   }

   const REAL operator[](const INDEX i) const { assert(i<3); return 0.0; }
   const INDEX size() const { return 3; }

   template<typename REPAM_ARRAY>
   REAL EvaluatePrimal(const REPAM_ARRAY& repam, const PrimalSolutionStorage::Element primal) const
   {
      assert(repam.size() == 3);
      const INDEX sum = INDEX(primal[0]) + INDEX(primal[1]) + INDEX(primal[2]);
      if(sum == 1) return std::numeric_limits<REAL>::max();
      else return repam[0]*primal[0] + repam[1]*primal[1] + repam[2]*primal[2];
   }
   void WritePrimal(const PrimalSolutionStorage::Element, std::ofstream& fs) const
   {
   }

   template<typename REPAM_POT, typename MSG_ARRAY>
   void MakeFactorUniform(const REPAM_POT& repam, MSG_ARRAY& msgs, const double omega) const
   {
      //const REAL omega = std::accumulate(omegaIt, omegaIt+3,0.0);
      //const REAL omega = 1.0; // do zrobienia: for now, in general this will not converge
      assert(omega <= 1.0 + eps);
      //std::cout << "(" << *omegaIt << "," << *(omegaIt+1) << "," << *(omegaIt+2) << ")\n";
      const auto sortIndices = SortIndices(repam);
      const REAL x = repam[sortIndices[0]] + repam[sortIndices[1]];

      assert(msgs.size() == 3);
      // do zrobienia: stupid interface. Make references. See also in factors_messages.hxx

      if(x > 0.0) { // labeling 000
         if(repam[sortIndices[0]] < 0.0) {
            //msgs[sortIndices[2]]->operator[](0) += omega*(-repam[sortIndices[2]] + repam[sortIndices[1]] - x);
            //msgs[sortIndices[1]]->operator[](0) -= omega*x;
            const REAL delta = repam[sortIndices[1]] - 0.5*x;
            //msgs[sortIndices[0]]->operator[](0) -= 0.5*omega*x;
            //msgs[sortIndices[1]]->operator[](0) -= 0.5*omega*x;
            //msgs[sortIndices[2]]->operator[](0) += omega*(-repam[sortIndices[2]] + delta);
            msgs[sortIndices[0]] -= 0.5*omega*x;
            msgs[sortIndices[1]] -= 0.5*omega*x;
            msgs[sortIndices[2]] += omega*(-repam[sortIndices[2]] + delta);
         } else {
            for(INDEX i=0; i<3; ++i) {
               //msgs[i]->operator[](0) -= omega*repam[i];
               msgs[i] -= omega*repam[i];
            }
         }
      } else if(repam[sortIndices[2]] > 0.0) { //labeling 110
         if(repam[sortIndices[1]] > 0.0) {
            //msgs[sortIndices[0]]->operator[](0) -= omega*x;
            //msgs[sortIndices[2]]->operator[](0) += omega*(-repam[sortIndices[2]] + repam[sortIndices[1]]);
            // do zrobiebia: make decrease of repam[sortIndices[2]] and increase of repam[sortIndices[1]] as equal as pssible.
            const REAL b = std::min(-0.5*x,0.5*(repam[sortIndices[2]] - repam[sortIndices[1]]));
            const REAL delta = repam[sortIndices[1]] + b;
            //msgs[sortIndices[1]]->operator[](0) += omega*b;
            //msgs[sortIndices[2]]->operator[](0) += omega*(-repam[sortIndices[2]] + delta);
            //msgs[sortIndices[0]]->operator[](0) -= omega*(0.5*x + b); // do zrobienia: can possibly be heightened more
            msgs[sortIndices[1]] += omega*b;
            msgs[sortIndices[2]] += omega*(-repam[sortIndices[2]] + delta);
            msgs[sortIndices[0]] -= omega*(0.5*x + b); // do zrobienia: can possibly be heightened more
         } else {
            for(INDEX i=0; i<3; ++i) {
               //msgs[i]->operator[](0) -= omega*repam[i];
               msgs[i] -= omega*repam[i];
            }
         }
      } else { // labeling 111
         for(INDEX i=0; i<3; ++i) {
            //msgs[i]->operator[](0) -= omega*repam[i];
            msgs[i] -= omega*repam[i];
         }
      }
   }


private:
};

} // end namespace LP_MP

#endif // LP_MP_MULTICUT_TRIPLET_FACTOR



