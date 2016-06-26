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

   // shortcut to indicate how big the message is
   constexpr SIGNED_INDEX variableMessageSize = -1;
}

template class MinCost<LP_MP::SIGNED_INDEX,LP_MP::REAL>;

//template class MinCost<int,int>;
//template class MinCost<int,size_t>;
//template class MinCost<int,double>;
//template class MinCost<int,float>;

#endif // LP_MP_CONFIG_HXX

