/*
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
   template<typename REPAM_ARRAY>
   IndicesType SortIndices(const REPAM_ARRAY& repamPot) const
   {
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
   }

   constexpr static INDEX size() { return 3; }

   template<typename REPAM_ARRAY, typename PRIMAL>
   REAL EvaluatePrimal(const REPAM_ARRAY& repam, const PRIMAL primal) const
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

      
      //
      if(x > 0.0) { // labeling 000
         if(repam[sortIndices[0]] < 0.0) {
            const REAL delta = repam[sortIndices[1]] - 0.5*x;
            msgs[sortIndices[0]] -= 0.5*omega*x;
            msgs[sortIndices[1]] -= 0.5*omega*x;
            msgs[sortIndices[2]] += omega*(-repam[sortIndices[2]] + delta);
         } else {
            for(INDEX i=0; i<3; ++i) {
               msgs[i] -= omega*repam[i];
            }
         }
      } else if(repam[sortIndices[2]] > 0.0) { //labeling 110
         if(repam[sortIndices[1]] > 0.0) {
            const REAL b = std::min(-0.5*x,0.5*(repam[sortIndices[2]] - repam[sortIndices[1]]));
            const REAL delta = repam[sortIndices[1]] + b;
            msgs[sortIndices[1]] += omega*b;
            msgs[sortIndices[2]] += omega*(-repam[sortIndices[2]] + delta);
            msgs[sortIndices[0]] -= omega*(0.5*x + b); // do zrobienia: can possibly be heightened more
         } else {
            for(INDEX i=0; i<3; ++i) {
               msgs[i] -= omega*repam[i];
            }
         }
      } else { // labeling 111
         for(INDEX i=0; i<3; ++i) {
            msgs[i] -= omega*repam[i];
         }
      }
   }
};


// left factor must be MulticutUnaryFactor and right factor must be MulticutTripletFactor
// possibly templatize for index i_
// templatize for either SRMP- or MPLP-type message passing, i.e. uanry factors are active or triplet factors are active
enum class MessageSending { SRMP, MPLP }; // do zrobienia: place this possibly more global, also applies to pairwise factors in MRFs
template<MessageSending MST = MessageSending::SRMP>
class MulticutUnaryTripletMessage
{
public:
   MulticutUnaryTripletMessage(const INDEX i) : i_(i) { assert(i < 3); }; // i is the index in the triplet factor

   template<typename RIGHT_FACTOR, typename G1, typename G2, MessageSending MST_TMP = MST>
   typename std::enable_if<MST_TMP == MessageSending::SRMP,void>::type
   ReceiveMessageFromRight(RIGHT_FACTOR* const r, const G1& rightPot, G2& msg) 
   {
      static_assert(MST_TMP == MST,"");
      const auto sortIndices = r->SortIndices(rightPot);
      // now adjust costs such that labeling just stays optimal for modified rightPot
      const REAL x = rightPot[sortIndices[0]] + rightPot[sortIndices[1]];

      // investigate all possiblities: labelings + position in labeling -> 3x3 decisions
      // do zrobienia: devise one formula which takes into account all those decisions. faster?
      if(x > 0.0) { // labeling 000
         if(sortIndices[0] == i_) {
            msg[0] += -x;
         } else if(sortIndices[1] == i_) {
            msg[0] += -x;
         } else {
            msg[0] += -rightPot[i_] - rightPot[sortIndices[0]];
         }
      } else if(rightPot[sortIndices[2]] > 0) { // labeling 110
         if(sortIndices[0] == i_) {
            msg[0] += -rightPot[i_] + std::min(-rightPot[sortIndices[1]], rightPot[sortIndices[2]]);
         } else if(sortIndices[1] == i_) {
            msg[0] += -rightPot[i_] + std::min(-rightPot[sortIndices[0]], rightPot[sortIndices[2]]);
         } else {
            msg[0] += -rightPot[i_] + std::max(rightPot[sortIndices[1]], 0.0);
         }
      } else { // labeling 111 -> all reparametrized values <= 0
         if(sortIndices[0] == i_) {
            msg[0] -= rightPot[i_];
         } else if(sortIndices[1] == i_) {
            msg[0] -= rightPot[i_];
         } else {
            msg[0] -= rightPot[i_];
         }
      }
   }

   template<typename LEFT_FACTOR, typename RIGHT_FACTOR, typename G1, typename G2, typename G3, MessageSending MST_TMP = MST>
   typename std::enable_if<MST_TMP == MessageSending::SRMP,void>::type
   SendMessageToRight(LEFT_FACTOR* const l, RIGHT_FACTOR* const r, const G1& leftPot, const G2& rightPot, G3& msg, const REAL omega)
   {
      static_assert(MST_TMP == MST,"");
      msg[0] += omega*leftPot[0];
   }

   template<typename LEFT_FACTOR, typename G1, typename G2, MessageSending MST_TMP = MST>
   typename std::enable_if<MST_TMP == MessageSending::MPLP,void>::type
   ReceiveMessageFromLeft(LEFT_FACTOR* const l, const G1& leftPot, G2& msg) 
   {
      static_assert(MST_TMP == MST,"");
      msg[0] += leftPot[0];
   }

   template<typename RIGHT_FACTOR, typename MSG_ARRAY, typename RIGHT_REPAM, typename ITERATOR, MessageSending MST_TMP = MST>
   static typename std::enable_if<MST_TMP == MessageSending::MPLP,void>::type
   SendMessagesToLeft(const RIGHT_FACTOR& rightFactor, const RIGHT_REPAM& rightRepam, const MSG_ARRAY& msgs, ITERATOR omegaIt)
   {
      static_assert(MST_TMP == MST,"");
      assert(msgs[0]->GetMessageOp().i_ == 0);
      assert(msgs[1]->GetMessageOp().i_ == 1);
      assert(msgs[2]->GetMessageOp().i_ == 2);
      const REAL omega = std::accumulate(omegaIt, omegaIt+3,0.0);
      //const REAL omega = 1.0; // do zrobienia: for now, in general this will not converge
      rightFactor.MakeFactorUniform(rightRepam, msgs, omega);
    }

   template<typename G>
   void RepamLeft(G& repamPot, const REAL msg, const INDEX msg_dim)
   {
      assert(msg_dim == 0);
      repamPot[0] -= msg;
   }
   template<typename G>
   void RepamRight(G& repamPot, const REAL msg, const INDEX msg_dim)
   {
      assert(msg_dim == 0);
      repamPot[i_] += msg;
   }

   // compute primal functions, how to do it? look into ultra-metric rounding
   void ComputeRightFromLeftPrimal(const bool leftPrimal, MulticutTripletFactor::LabelingType& rightPrimal)
   {
      rightPrimal[i_] = leftPrimal;
   }


private:
   const INDEX i_; // index of the affected variable in the cycle factor, do zrobienia: needs only two bits
};





} // end namespace LP_MP

#endif // LP_MP_MULTICUT_TRIPLET_FACTOR
*/



#ifndef LP_MP_MULTICUT_TRIPLET_FACTOR_HXX
#define LP_MP_MULTICUT_TRIPLET_FACTOR_HXX

#include "LP_MP.h"

namespace LP_MP {

// encoding with 4 entries corresponding to the four possible labelings in this order:
// 011 101 110 111. Labeling 000 is always assigned cost 0. Labelings 100 010 001 are forbidden.
class MulticutTripletFactor
{
public:
   MulticutTripletFactor() {};

   using TripletEdges = std::array<std::array<INDEX,2>,3>;
   static TripletEdges SortEdges(const INDEX i1, const INDEX i2, const INDEX i3)
   {
      std::array<INDEX,3> ti{i1,i2,i3}; // the node indices of the triplet
      std::sort(ti.begin(), ti.end()); // do zrobienia: use faster sorting
      assert(ti[0] < ti[1] < ti[2]);
      TripletEdges te{{{ti[0],ti[1]},{ti[0],ti[2]},{ti[1],ti[2]}}}; 
      return te;
   }

   template<typename REPAM_ARRAY>
   void MaximizePotential(const REPAM_ARRAY& repam) {};

   constexpr static INDEX size() { return 4; }

   template<typename REPAM_ARRAY>
   REAL LowerBound(const REPAM_ARRAY& repamPot) const 
   {
      assert(repamPot.size() == 4);
      return std::min(std::min(std::min( repamPot[0], repamPot[1]), std::min(repamPot[2], repamPot[3])),0.0); // possibly, if we have a SIMD factor, use a simdized minimum
   }

   template<typename REPAM_ARRAY, typename PRIMAL>
   REAL EvaluatePrimal(const REPAM_ARRAY& repam, const PRIMAL primal) const
   {
      assert(repam.size() == 4);
      return 3e111;
   }

};


// message between unary factor and triplet factor
//enum class MessageSending { SRMP, MPLP }; // do zrobienia: place this possibly more global, also applies to pairwise factors in MRFs
template<MessageSendingType MST = MessageSendingType::SRMP>
class MulticutUnaryTripletMessage
{
public:
   // i is the index in the triplet factor
   MulticutUnaryTripletMessage(const INDEX i) : i_(i) 
   { 
      assert(i < 3); 
   }; 
   ~MulticutUnaryTripletMessage()
   {
      static_assert(MST == MessageSendingType::SRMP,"");
   }

   constexpr static INDEX size() { return 1; }

   template<typename RIGHT_FACTOR, typename G1, typename G2, MessageSendingType MST_TMP = MST>
   //typename std::enable_if<MST_TMP == MessageSending::SRMP,void>::type
   void
   ReceiveMessageFromRight(RIGHT_FACTOR* const r, const G1& rightPot, G2& msg) 
   {
      assert(msg.size() == 1);
      msg[0] -= std::min(std::min(rightPot[(i_+1)%3], rightPot[(i_+2)%3]), rightPot[3]) - std::min(rightPot[i_],0.0);
   }

   template<typename LEFT_FACTOR, typename G1, typename G3, MessageSendingType MST_TMP = MST>
   //typename std::enable_if<MST_TMP == MessageSending::SRMP,void>::type
   void
   SendMessageToRight(LEFT_FACTOR* const l, const G1& leftPot, G3& msg, const REAL omega)
   {
      assert(msg.size() == 1);
      static_assert(MST_TMP == MST,"");
      msg[0] -= omega*leftPot[0];
   }

   template<typename G>
   void RepamLeft(G& repamPot, const REAL msg, const INDEX msg_dim)
   {
      assert(msg_dim == 0);
      repamPot[0] += msg;
   }
   template<typename G>
   void RepamRight(G& repamPot, const REAL msg, const INDEX msg_dim)
   {
      assert(msg_dim == 0);
      // the labeling with two 1s in them
      repamPot[(i_+1)%3] += msg;
      repamPot[(i_+2)%3] += msg;
      // for the 111 labeling which is always affected
      repamPot[3] += msg; 
   }

   // compute primal functions, how to do it? look into ultra-metric rounding
   void ComputeRightFromLeftPrimal(const bool leftPrimal, typename PrimalSolutionStorage::Element rightPrimal)
   {
      // do zrobienia
   }

private:
   const INDEX i_;
};

} // end namespace LP_MP

#endif // LP_MP_MULTICUT_TRIPLET_FACTOR

