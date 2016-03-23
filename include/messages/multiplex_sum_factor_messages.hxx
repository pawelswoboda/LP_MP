#ifndef LP_MP_MULTIPLEX_SUM_FACTOR_MESSAGES_HXX
#define LP_MP_MULTIPLEX_SUM_FACTOR_MESSAGES_HXX

#include "template_utilities.hxx"

namespace LP_MP {

//////////////////////////////////////////////////////////////////////////////////////////////////
// kwaskwas: do zrobienia: this class is not working. Same for MultiplexSumFactor. Remove them. //
//////////////////////////////////////////////////////////////////////////////////////////////////

// two message classes need to be defined: One from the simplices to MultiplexSumFactor
// the other from MultiplexSumFactor to summated MultiplexFactor

// for both message classes: check if left and right factors are of correct type, e.g. Multiplex and MultiplexSum

// left factor is multiplex factor, right one is MultiplexSumFactor, connects to i-th individual multiplex
template<class MESSAGE_TYPE> // the underlying message, e.g. equality message between two multiplex factors, to be transformed into a message between MultiplexFactor and MultiplexSumFactor
class MultiplexToMultiplexSummandMessage {
public:
   using BaseMessageType = MESSAGE_TYPE;

   MultiplexToMultiplexSummandMessage(const BaseMessageType& baseMessage, const INDEX multiplexSumFactorIndex)
      : baseMessage_(baseMessage),
      multiplexSumFactorIndex_(multiplexSumFactorIndex)
   {}

   template<typename RIGHT_FACTOR, typename G1, typename G2>
   void ReceiveMessageFromRight(RIGHT_FACTOR* const r, const G1& rightPot, G2& msg)
   {
      baseMessage_.ReceiveMessageFromRight(r->GetSummandFactor(), GetSumPlusSummand(rightPot,r), msg);
      //baseMessage_.ReceiveMessageFromRight(r->GetSummandFactor(), GetSummand(rightPot,r), msg);
   }

   template<typename LEFT_FACTOR, typename G1, typename G2>
   void ReceiveMessageFromLeft(LEFT_FACTOR* l, const G1& leftPot, G2& msg)
   {
      baseMessage_.ReceiveMessageFromLeft(l, leftPot, msg);
   }

   template<typename LEFT_FACTOR, typename RIGHT_FACTOR, typename G1, typename G2, typename G3>
   void SendMessageToRight(LEFT_FACTOR* const l, RIGHT_FACTOR* const r, const G1& leftPot, const G2& rightPot, G3& msg, const REAL omega)
   {
      baseMessage_.SendMessageToRight(l, r->GetSummandFactor(), leftPot, GetSumPlusSummand(rightPot,r), msg, omega);
      //baseMessage_.SendMessageToRight(l, r->GetSummandFactor(), leftPot, GetSummand(rightPot,r), msg, omega);
   }

   template<typename LEFT_FACTOR, typename RIGHT_FACTOR, typename G1, typename G2, typename G3>
   void SendMessageToLeft(LEFT_FACTOR* const l, RIGHT_FACTOR* const r, const G1& leftPot, const G2& rightPot, G3& msg, const REAL omega)
   {
      baseMessage_.SendMessageToRight(l, r->GetSummandFactor(), leftPot, GetSumPlusSummand(rightPot,r), msg, omega);
      //baseMessage_.SendMessageToLeft(l, r->GetSummandFactor(), leftPot, GetSummand(rightPot,r), msg, omega);
   }

   template<typename REPAM_ARRAY>
   void RepamLeft(REPAM_ARRAY& leftRepamPot, const REAL msg, const INDEX dim)
   {  
      baseMessage_.RepamLeft(leftRepamPot, msg, dim);
   }
   template<typename REPAM_ARRAY>
   void RepamRight(REPAM_ARRAY& rightRepamPot, const REAL msg, const INDEX dim)
   {
      auto r = GetSummand(rightRepamPot, rightRepamPot.GetFactor());
      baseMessage_.RepamRight(r, msg, dim);
   }
   
private:

   template<typename ARRAY_TYPE, typename FACTOR_TYPE>
   VecSlice<ARRAY_TYPE> GetSummand(ARRAY_TYPE& v, FACTOR_TYPE* f) 
   {
      const INDEX stride = f->GetMultiplexDim();
      return VecSlice<ARRAY_TYPE>(v, multiplexSumFactorIndex_*stride, (multiplexSumFactorIndex_+1)*stride);
   }

   /*
   template<typename ARRAY_TYPE, typename MULTIPLEX_SUM_FACTOR>
   VecSlice<ARRAY_TYPE> GetSum(ARRAY_TYPE& v, const MULTIPLEX_SUM_FACTOR& f)
   {
      const INDEX multiplexNo = f->GetMultiplexNo();
      const INDEX stride = f->GetMultiplexDim();
      std::cout << multiplexNo << ", " << stride << "\n";
      return VecSlice<ARRAY_TYPE>(v, multiplexNo*stride, (multiplexNo+1)*stride);
   }
   */

   template<typename ARRAY_TYPE, typename MULTIPLEX_SUM_FACTOR>
   //PlusExprVec<VecSlice<ARRAY_TYPE>,VecSlice<ARRAY_TYPE>> 
   std::vector<REAL>
   GetSumPlusSummand(ARRAY_TYPE& v, const MULTIPLEX_SUM_FACTOR& f)
   {
      const INDEX multiplexNo = f->GetMultiplexNo();
      const INDEX stride = f->GetMultiplexDim();

      //VecSlice<ARRAY_TYPE> summand = GetSummand(v,f);
      //VecSlice<ARRAY_TYPE> sum = GetSum(v,f);
      std::vector<REAL> s(stride,0.0);
      for(INDEX i=0; i<s.size(); ++i) {
         s[i] = v[multiplexSumFactorIndex_*stride + i] + v[multiplexNo*stride+i];
      }
      return s;
      //return PlusExprVec<decltype(summand), decltype(sum)> (summand,sum);
   }

   BaseMessageType baseMessage_; // possibly share? Should be the same for all summands. Could be stored in MultiplexSumFactor.
   const INDEX multiplexSumFactorIndex_;
};

// left factor is MultiplexSumFactor, right one is Multiplex, connects from the sum multiplex in MultiplexSumFactor
// message is of length (1 + number of multiplexes which are summed) x (size of multiplexes)
template<class MESSAGE_TYPE> // the underlying message, e.g. equality message between two multiplex factors, to be transformed into a message between MultiplexFactor and MultiplexSumFactor
class MultiplexSumToMultiplexMessage {
   using BaseMessageType = MESSAGE_TYPE;
public:
   MultiplexSumToMultiplexMessage(const BaseMessageType& baseMessage)
      : baseMessage_(baseMessage)
   {}

   template<typename RIGHT_FACTOR, typename G1, typename G2>
   void ReceiveMessageFromRight(RIGHT_FACTOR* const r, const G1& rightPot, G2& msg)
   {
      baseMessage_.ReceiveMessageFromRight(r, rightPot, msg);
      return;
      // do zrobienia: this is all dirty tricks here
      const INDEX multiplexDim = 2;
      const INDEX multiplexNo = msg.size()/2;
      const REAL frac = 1.0/REAL(multiplexNo);
      std::vector<REAL> scaledRightPot(rightPot.size()); // do zrobienia: must be temporarily stored, as rightPot is modified all the time
      for(INDEX i=0; i<scaledRightPot.size(); ++i) {
         scaledRightPot[i] = frac*rightPot[i];
      }
      for(INDEX i=0; i<multiplexNo; ++i) {
         baseMessage_.ReceiveMessageFromRight(r, scaledRightPot, msg);
      }
   }

   template<typename LEFT_FACTOR, typename G1, typename G2>
   void ReceiveMessageFromLeft(LEFT_FACTOR* l, const G1& leftPot, G2& msg)
   {
      const INDEX multiplexDim = l->GetMultiplexDim();
      const INDEX multiplexNo = l->GetMultiplexNo();
      assert(msg.size() == multiplexDim*multiplexNo);
      baseMessage_.ReceiveMessageFromLeft(l, , msg);
      return;
      for(INDEX i=0; i<multiplexNo; ++i) {
         auto curLeftPot = GetSummand(leftPot, leftPot.GetFactor(), i);
         baseMessage_.ReceiveMessageFromLeft(l, curLeftPot, msg);
      }
      
      //baseMessage_.ReceiveMessageFromLeft(l, GetSum(leftPot, l), msg);
   }

   template<typename LEFT_FACTOR, typename RIGHT_FACTOR, typename G1, typename G2, typename G3>
   void SendMessageToRight(LEFT_FACTOR* const l, RIGHT_FACTOR* const r, const G1& leftPot, const G2& rightPot, G3& msg, const REAL omega)
   {
      return;
      const INDEX multiplexDim = l->GetMultiplexDim();
      const INDEX multiplexNo = l->GetMultiplexNo();
      for(INDEX i=0; i<multiplexNo; ++i) {
         auto curLeftPot = GetSummand(leftPot, l, i);
         baseMessage_.SendMessageToRight(l, r, curLeftPot, rightPot, msg, omega); // do zrobienia: this will not work in general, just for MultiplexMargMessage, which does not depend on the right factor
      }
      //baseMessage_.SendMessageToRight(l->GetSumFactor(), r, GetSum(leftPot, l), rightPot, msg, omega);
   }

   template<typename LEFT_FACTOR, typename RIGHT_FACTOR, typename G1, typename G2, typename G3>
   void SendMessageToLeft(LEFT_FACTOR* const l, RIGHT_FACTOR* const r, const G1& leftPot, const G2& rightPot, G3& msg, const REAL omega)
   {
      return;
      const INDEX multiplexDim = l->GetMultiplexDim();
      const INDEX multiplexNo = l->GetMultiplexNo();
      const REAL omega_frac = omega/REAL(multiplexNo);
      const REAL frac = 1/REAL(multiplexNo);

      std::vector<REAL> freezedRightPot(rightPot.size()); // do zrobienia: must be temporarily stored, as rightPot is modified all the time
      for(INDEX i=0; i<freezedRightPot.size(); ++i) {
         freezedRightPot[i] = rightPot[i];
      }

      for(INDEX i=0; i<multiplexNo; ++i) {
         auto curLeftPot = GetSummand(leftPot, l, i);
         baseMessage_.SendMessageToLeft(l->GetSummandFactor(), r, curLeftPot, freezedRightPot, msg, omega_frac);
      }
      //baseMessage_.SendMessageToLeft(l->GetSumFactor(), r, GetSum(leftPot, l), rightPot, msg, omega);
   }

   template<typename REPAM_ARRAY>
   void RepamLeft(REPAM_ARRAY& leftRepamPot, const REAL msg, const INDEX dim)
   {
      //std::cout << "kwas before: " << "(" << leftRepamPot[leftRepamPot.size()-2] << ", " << leftRepamPot[leftRepamPot.size()-1] << ")" << ", " << msg << "\n";
      //std::vector<REAL> kwas(2,0.0);
      VecSlice<REPAM_ARRAY> s = GetSum(leftRepamPot, leftRepamPot.GetFactor());
      baseMessage_.RepamLeft(s, msg, dim);
      //baseMessage_.RepamLeft(kwas, msg, dim);
      //std::cout << "kwas  after: " << "(" << leftRepamPot[leftRepamPot.size()-2] << ", " << leftRepamPot[leftRepamPot.size()-1] << ")" << ", " << msg << "\n";
      //std::cout << "kwas  test1: " << "(" << kwas[0] << ", " << kwas[1] << ")" << ", " << msg << "\n";
      //std::cout << "kwas  test2: " << "(" << s[0] << ", " << s[1] << ")" << ", " << msg << "\n";
      //std::cout << "kwas test full: ";
      //for(INDEX i=0; i<leftRepamPot.size(); ++i) {
      //   std::cout << leftRepamPot[i] << ",";
      //}
      //std::cout << "\n";
      return;
      // do zrobienia: this could be made faster by reparametrizing everything at once
      const INDEX multiplexDim = leftRepamPot.GetFactor()->GetMultiplexDim();
      const INDEX multiplexNo = dim/multiplexDim;
      auto l = GetSummand(leftRepamPot, leftRepamPot.GetFactor(), multiplexNo);
      baseMessage_.RepamLeft(l, msg, dim%multiplexDim);
      /*
      for(INDEX i=0; i<leftRepamPot.GetFactor()->GetMultiplexNo(); ++i) {
         auto l = GetSummand(leftRepamPot,leftRepamPot.GetFactor(),i);
         baseMessage_.RepamLeft(l, msg, dim);
      }
      */
   }
   template<typename REPAM_ARRAY>
   void RepamRight(REPAM_ARRAY& rightRepamPot, const REAL msg, const INDEX dim)
   {
      baseMessage_.RepamRight(rightRepamPot, msg, dim);
      return;
      const INDEX multiplexDim = 2; // do zrobienia: this need not hold in general!
      baseMessage_.RepamRight(rightRepamPot, msg, dim%multiplexDim);
   }
   
private:

   template<typename ARRAY_TYPE, typename FACTOR_TYPE>
   VecSlice<ARRAY_TYPE> GetSummand(ARRAY_TYPE& v, FACTOR_TYPE* f, const INDEX i) 
   {
      assert(i < f->GetMultiplexNo());
      const INDEX stride = f->GetMultiplexDim();
      return VecSlice<ARRAY_TYPE>(v, i*stride, (i+1)*stride);
   }

   template<typename ARRAY_TYPE, typename MULTIPLEX_SUM_FACTOR>
   VecSlice<ARRAY_TYPE> GetSum(ARRAY_TYPE& v, const MULTIPLEX_SUM_FACTOR& f)
   {
      const INDEX multiplexNo = f->GetMultiplexNo();
      const INDEX stride = f->GetMultiplexDim();
      return VecSlice<ARRAY_TYPE>(v, multiplexNo*stride, (multiplexNo+1)*stride);
   }
   /*
   template<typename ARRAY_TYPE, typename MULTIPLEX_SUM_FACTOR>
   std::vector<REAL> GetSum(ARRAY_TYPE& v, const MULTIPLEX_SUM_FACTOR& f)
   {
      //std::cout << f->GetMultiplexDim() << "\n";
      std::vector<REAL> sum(f->GetMultiplexDim(), 0.0);
      INDEX stride = f->GetMultiplexDim();
      std::cout << f->size()/sum.size() << "=" << f->GetMultiplexNo() << "\n";
      const INDEX no_multiplex = f->GetMultiplexNo();
      for(INDEX multiplex_idx=0; multiplex_idx<no_multiplex; ++multiplex_idx) {
         for(INDEX i=0; i<sum.size(); ++i) { // entry in multiplex
            sum[i] += v[multiplex_idx*stride + i];
         }
      }
      return sum;
   }
   */

   BaseMessageType baseMessage_;
};

}

#endif // LP_MP_MULTIPLEX_SUM_FACTOR_MESSAGES_HXX
