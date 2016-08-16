#ifndef LP_MP_CONFIG_HXX
#define LP_MP_CONFIG_HXX

#include "MinCost/MinCost.h"
#include "Vc/Vc"
//#include "Vc/Memory"

// type definitions for LP_MP

#ifdef _MSC_VER
#pragma warning(disable: 4661)
#endif



namespace LP_MP {

   // data types for all floating point/integer operations 
   using REAL = double;
   using INDEX = unsigned int;
   using UNSIGNED_INDEX = INDEX;
   using SIGNED_INDEX = int; // note: must be the same as flow type in MinCost
   using SHORT_INDEX = unsigned char;
   using LONG_SIGNED_INDEX = long int;
   using LONG_INDEX = long unsigned int;

   // data types for all floating point/integer operations performed with SIMD
   using REAL_SIMD = Vc::double_v;
   using REAL_MASK_SIMD = Vc::double_m;
   using INDEX_SIMD = Vc::int_v;
   using INDEX_MASK_SIMD = Vc::int_m;

   enum class Chirality {left,right};
   enum class MessageSendingType {SRMP,MPLP};
   enum class ReparametrizationMode {Explicit,Implicit};

   constexpr REAL eps = 1e-8;
   
   // shortcuts to indicate how many messages a factor holds
   constexpr SIGNED_INDEX variableMessageNumber = 0;
   constexpr SIGNED_INDEX atMostOneMessage = -1;
   constexpr SIGNED_INDEX atMostTwoMessages = -2;
   constexpr SIGNED_INDEX atMostThreeMessages = -3;
   constexpr SIGNED_INDEX atMostFourMessages = -4;
   constexpr SIGNED_INDEX atMostFiveMessages = -5;
   constexpr SIGNED_INDEX atMostSixMessages = -6;

   // shortcut to indicate how big the message is: here it is determined only at runtime
   constexpr SIGNED_INDEX variableMessageSize = -1;


   // possibly make better class out of this having multiple fields|flags
   /*
   enum class LPVisitorReturnType {
      ReparametrizeUniform,ReparametrizeLowerBoundUniform,ReparametrizeLowerBoundPrimalUniform,ReparametrizePrimalUniform,
      ReparametrizeAnisotropic,ReparametrizeLowerBoundAnisotropic,ReparametrizeLowerBoundPrimalAnisotropic,ReparametrizePrimalAnisotropic,
      Break,Error
   };
   bool Round(const LpVisitorReturnType t) {
      if(t == ReparametrizeLowerBoundPrimalUniform ||
         t == ReparametrizePrimalUniform || 
         t == ReparametrizeLowerBoundPrimalAnisotropic || 
         t == ReparametrizePrimalAnisotropic) {
         return true;
      } else {
         return false;
      }
   }
   LpVisitorReturnType RemoveRounding(const LpVisitorReturnType t) {
      switch(t) {
         case ReparametrizeUniform:
            return ReparametrizeUniform;
            break;
         case ReparametrizeLowerBoundUniform:
            return ReparametrizeLowerBoundUniform;
            break;
         case ReparametrizeLowerBoundPrimalUniform:
            return ReparametrizeLowerBoundUniform
            break;
         case ReparametrizePrimalUniform:
            return ReparametrizeUniform
            break;
         case ReparametrizeAnisotropic:
            return ReparametrizeAnisotropic;
            break;
         case ReparametrizeLowerBoundAnisotropic:
            return ReparametrizeLowerBoundAnisotropic;
            break;
         case ReparametrizeLowerBoundPrimalAnisotropic:
            return ReparametrizeLowerBoundAnisotropic;
            break;
         case ReparametrizePrimalAnisotropic:
            return ReparametrizeAnisotropic;
            break;
         case Break:
            return Break;
            break;
         case Error:
            return Error;
            break;
         default:
            assert(false); // unknown case
      }
   }
   */
   // do zrobienia: maybe put this into LP_MP.h
   enum class LPReparametrizationMode {Anisotropic, Uniform, Undefined};

   // steers optimization of LP solver. Is returned by visitor and processed by solver.
   // do zrobienia: nicer design by implementing it via callbacks to solver.
   // also put this into solver.hxx
   //template<typename SOLVER>
   class LpControl {
   public:
      //void ScheduleTighten(const INDEX noConstraints, const REAL minDualIncrease) 
      //{ 
      //   tighten = true; 
      //   tightenConstraints = noConstraints;
      //   minDualIncrease = minDualIncrease; 
      //}
      //void SchedulePrimal() { computePrimal = true; }
      //void SchedulaLowerBound() { computeLowerBound = true; }
      //void SetReparametrization(LPReparametrizationMode r) { repam_ = r; }
      //void End() { end = true; }
   //private:
      LPReparametrizationMode repam = LPReparametrizationMode::Undefined;
      bool computePrimal = false;
      bool computeLowerBound = false;
      bool tighten = false;
      bool end = false; // terminate optimization
      bool error = false;
      INDEX tightenConstraints = 0; // when given as return type, indicates how many constraints are to be added. When given as parameter to visitor, indicates how many were added.
      REAL tightenMinDualIncrease = 0.0;
   };


   // hash function for various types
   namespace hash {
      static auto array2 = [](const std::array<INDEX,2> x) { return std::hash<INDEX>()(x[0])^std::hash<INDEX>()(x[1]); };
      static auto array3 = [](const std::array<INDEX,3> x) { return std::hash<INDEX>()(x[0])^std::hash<INDEX>()(x[1])^std::hash<INDEX>()(x[2]); };
      static auto array4 = [](const std::array<INDEX,4> x) { return std::hash<INDEX>()(x[0])^std::hash<INDEX>()(x[1])^std::hash<INDEX>()(x[2])^std::hash<INDEX>()(x[3]); };
   }
}

template class MinCost<LP_MP::SIGNED_INDEX,LP_MP::REAL>;

//template class MinCost<int,int>;
//template class MinCost<int,size_t>;
//template class MinCost<int,double>;
//template class MinCost<int,float>;

#endif // LP_MP_CONFIG_HXX

