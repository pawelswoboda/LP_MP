#ifndef LP_MP_INSTANCES_INC
#define LP_MP_INSTANCES_INC

#include "MinCost/MinCost.h"
#include "Vc/Vc"
#include "Vc/Memory"

// type definitions for LP_MP

#ifdef _MSC_VER
#pragma warning(disable: 4661)
#endif


template class MinCost<int,int>;
template class MinCost<int,size_t>;
template class MinCost<int,double>;
template class MinCost<int,float>;

namespace LP_MP {

   // data type for all floating poINDEX/INDEXeger operations 
   using REAL = double;
   using INDEX = size_t;
   using SHORT_INDEX = unsigned char;
   // data type for all floating poINDEX/INDEXeger operations performed with SIMD
   using REAL_SIMD = Vc::double_v; // albo Vc::float_v
   using REAL_MASK_SIMD = Vc::double_m;
   using INDEX_SIMD = Vc::int_v;
   using INDEX_MASK_SIMD = Vc::int_m;

   enum class Chirality {left,right};

}

#endif // LP_MP_INSTANCES_INC
