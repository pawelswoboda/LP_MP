#ifndef LP_MP_CONFIG_HXX
#define LP_MP_CONFIG_HXX

//#include "MinCost/MinCost.h"
//#include "Vc/Vc"
//#include "Vc/Memory"

//#include <fenv.h>
#include <stdexcept>
#include <string>
#include <array>
#include <cmath>
#include <cassert>
#include <limits>

// type definitions for LP_MP

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
//   using REAL_SIMD = Vc::double_v;
//   using REAL_MASK_SIMD = Vc::double_m;
//   using INDEX_SIMD = Vc::int_v;
//   using INDEX_MASK_SIMD = Vc::int_m;

   enum class Chirality {left,right};
   enum class MessageSendingType {SRMP,MPLP}; // also add full, for always sending and receiving messages

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

   // do zrobienia: maybe put this into LP_MP.h
   struct LpReparametrizationModetest {
      enum class LPReparametrizationMode {Anisotropic, Uniform, Undefined} mode;
      REAL weight;
   };
   enum class LPReparametrizationMode {Anisotropic, Uniform, Undefined};

   inline LPReparametrizationMode LPReparametrizationModeConvert(const std::string& s)
   {
      //feenableexcept(FE_INVALID | FE_OVERFLOW);
      const std::string uniform = "uniform";
      if(s == "anisotropic") {
         //return LpReparametrizationMode({mode::anisotropic,0.0});
         return LPReparametrizationMode::Anisotropic;
      } else if(s.length() >= uniform.length() && std::equal(uniform.begin(), uniform.end(), s.begin())) {
         //if(s.length() == uniform.length()) {
         //   return LpReparametrizationMode({mode::uniform,0.0});
         //} else {
         //   const REAL weight = std::stod( s.substr(uniform.length()) );
         //   return LpReparametrizationMode({mode::uniform,weight});
         //}
         return LPReparametrizationMode::Uniform;
      } else {
         throw std::runtime_error("reparametrization mode " + s + " unknown");
      }
   }

   // steers optimization of LP solver. Is returned by visitor and processed by solver.
   // also put this into solver.hxx
   class LpControl {
   public:
      LPReparametrizationMode repam = LPReparametrizationMode::Undefined;
      bool computePrimal = false;
      bool computeLowerBound = false;
      bool tighten = false;
      bool end = false; // terminate optimization
      bool error = false;
      INDEX tightenConstraints = 0; // when given as return type, indicates how many constraints are to be added. When given as parameter to visitor, indicates how many were added.
      REAL tightenMinDualIncrease = 0.0; // do zrobienia: obsolete
   };


   // hash function for various types
   namespace hash {
      static auto array2 = [](const std::array<INDEX,2> x) { return std::hash<INDEX>()(x[0])^std::hash<INDEX>()(x[1]); };
      static auto array3 = [](const std::array<INDEX,3> x) { return std::hash<INDEX>()(x[0])^std::hash<INDEX>()(x[1])^std::hash<INDEX>()(x[2]); };
      static auto array4 = [](const std::array<INDEX,4> x) { return std::hash<INDEX>()(x[0])^std::hash<INDEX>()(x[1])^std::hash<INDEX>()(x[2])^std::hash<INDEX>()(x[3]); };
   }

   REAL normalize(const REAL x) {
      assert(!std::isnan(x));
      if(std::isfinite(x)) {
         return x;
      } else {
         return std::numeric_limits<REAL>::infinity();
      }
   }

}

//template class MinCost<LP_MP::SIGNED_INDEX,LP_MP::REAL>;

//template class MinCost<int,int>;
//template class MinCost<int,size_t>;
//template class MinCost<int,double>;
//template class MinCost<int,float>;

#endif // LP_MP_CONFIG_HXX

