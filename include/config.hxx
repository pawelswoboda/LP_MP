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

   constexpr REAL eps = 1e-8;
}

template class MinCost<LP_MP::SIGNED_INDEX,LP_MP::REAL>;

//template class MinCost<int,int>;
//template class MinCost<int,size_t>;
//template class MinCost<int,double>;
//template class MinCost<int,float>;

#endif // LP_MP_CONFIG_HXX

