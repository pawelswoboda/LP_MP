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

// type definitions for LP_MP

namespace LP_MP {

   // data types for all floating point/integer operations 
   using REAL = float;
   using INDEX = unsigned int;
   using UNSIGNED_INDEX = INDEX;
   using SIGNED_INDEX = int; // note: must be the same as flow type in MinCost
   using SHORT_INDEX = unsigned char;
   using LONG_SIGNED_INDEX = long int;
   using LONG_INDEX = long unsigned int;

   enum class Chirality {left,right};
   enum class MessageSendingType {SRMP,MPLP}; // also add full, for always sending and receiving messages
   enum class Direction {forward, backward};

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
      size_t hash_combine( size_t lhs, size_t rhs ) {
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

   REAL normalize(const REAL x) {
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

//template class MinCost<LP_MP::SIGNED_INDEX,LP_MP::REAL>;

//template class MinCost<int,int>;
//template class MinCost<int,size_t>;
//template class MinCost<int,double>;
//template class MinCost<int,float>;

#endif // LP_MP_CONFIG_HXX

