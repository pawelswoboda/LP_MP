#ifndef LP_MP_CONFIG_HXX
#define LP_MP_CONFIG_HXX

//#include "MinCost/MinCost.h"

//#include <fenv.h>
#include <stdexcept>
#include <string>
#include <array>
#include <cmath>
#include <cassert>
#include <limits>
#include "tclap/CmdLine.h"

#define SIMDPP_ARCH_X86_AVX2
#include "simdpp/simd.h"

// type definitions for LP_MP

namespace LP_MP {

   // data types for all floating point/integer operations 
   // float is inaccurate for large problems and I observed oscillation. Possibly, some sort of numerical stabilization needs to be employed
   //using REAL = float;
   //constexpr std::size_t REAL_ALIGNMENT = 8;
   //using REAL_VECTOR = simdpp::float32<REAL_ALIGNMENT>;

   using REAL = double;
   constexpr std::size_t REAL_ALIGNMENT = 4;
   using REAL_VECTOR = simdpp::float64<REAL_ALIGNMENT>;

   using INDEX = std::size_t;
   using UNSIGNED_INDEX = INDEX;
   using SIGNED_INDEX = long int; // note: must be the same as flow type in MinCost
   using SHORT_INDEX = unsigned char;
   using LONG_SIGNED_INDEX = long int;
   using LONG_INDEX = long unsigned int;

   enum class Chirality {left,right};
   enum class MessageSendingType {SRMP,MPLP}; // TODO: remove
   enum class Direction {forward, backward};

   enum class message_passing_schedule {
     left, // messages are received from left and sent by left
     right, // messages are received from right and sent by right
     full, // messages are received and send in both directions
     only_send, // messages are only sent
     none // message is not called during message passing
   }; 

   constexpr REAL eps = std::is_same<REAL,float>::value ? 1e-6 : 1e-8;
   // verbosity levels: 0: silent
   //                   1: important diagnostics, e.g. lower bound, upper bound, runtimes
   //                   2: debug informations
   static INDEX verbosity = 0; 
   static bool diagnostics() { return verbosity >= 1; }
   static bool debug() { return verbosity >= 2; }
   
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
   enum class LPReparametrizationMode {Anisotropic, Anisotropic2, Uniform, DampedUniform, Mixed, Undefined};

   inline LPReparametrizationMode LPReparametrizationModeConvert(const std::string& s)
   {
      //feenableexcept(FE_INVALID | FE_OVERFLOW);
      const std::string uniform = "uniform";
      if(s == "anisotropic") {
         return LPReparametrizationMode::Anisotropic;
      } else if(s == "anisotropic2") {
         return LPReparametrizationMode::Anisotropic2;
      } else if(s == "uniform") {
         return LPReparametrizationMode::Uniform;
      } else if(s == "damped_uniform") {
         return LPReparametrizationMode::DampedUniform;
      } else if(s == "mixed") {
         return LPReparametrizationMode::Mixed;
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
      // equivalent of boost hash combine
      inline size_t hash_combine( size_t lhs, size_t rhs ) {
         lhs^= rhs + 0x9e3779b9 + (lhs << 6) + (lhs >> 2);
         return lhs;
      }

      template<typename T, size_t N>
      size_t hash_array(const std::array<T,N>& x)
      {
         size_t hash = std::hash<T>()(x[0]);
         for(INDEX i=0; i<N; ++i) {
            hash = hash_combine(hash, std::hash<T>()(x[i]));
         }
         return hash; 
      }
   }

   template<typename T>
   T normalize(const T x) {
      static_assert(std::is_same<T,double>::value || std::is_same<T,float>::value,"");
      assert(!std::isnan(x));
      if(std::isfinite(x)) {
         return x;
      } else {
         return std::numeric_limits<REAL>::infinity();
      }
   }

   // TCLAP constraints
   class PositiveRealConstraint : public TCLAP::Constraint<REAL>
   {
      public:
         std::string description() const { return "positive real constraint"; };
         std::string shortID() const { return "positive real number"; };
         bool check(const REAL& value) const { return value >= 0.0; };
   };
   class OpenUnitIntervalConstraint: public TCLAP::Constraint<REAL>
   {
      public:
         std::string description() const { return "0<x<1 real constraint"; };
         std::string shortID() const { return "positive real number smaller 1"; };
         bool check(const REAL& value) const { return value > 0.0 && value < 1.0; };
   };
   class PositiveIntegerConstraint : public TCLAP::Constraint<INDEX>
   {
      public:
         std::string description() const { return "strictly positive integer constraint"; };
         std::string shortID() const { return "strictly positive integer"; };
         bool check(const INDEX& value) const { return value > 0; };
   };
   static PositiveIntegerConstraint positiveIntegerConstraint;
}

// insert hash functions from above into standard namespace
namespace std
{
    template<size_t N> struct hash<std::array<LP_MP::INDEX,N>>
    {
        typedef std::array<LP_MP::INDEX,N> argument_type;
        typedef std::size_t result_type;
        result_type operator()(argument_type const& s) const
        {
            return LP_MP::hash::hash_array(s);
        }
    };
}

#endif // LP_MP_CONFIG_HXX

