#ifndef LP_MP_BINARY_TOMOGRAPHY
#define LP_MP_BINARY_TOMOGRAPHY

#include "factors_messages.hxx"
#include "LP_MP.h"
#include "problem_decomposition.hxx"
#include "factors/multiplex_factor.hxx"
#include "const_array_types.h"
#include "messages/multiplex_marg_message.hxx"

namespace LP_MP {

// file format
// problem 0 # mrf
// unary v 0 cost [0 0]
// pairwise v 0 1 cost [0 1; 1 0]
// ...
// problem 1 # counting constraints
// v 0 1 2 sum [${cost0} ${cost1} ${cost2}]
// ...

class BinaryTomography {
public:
  
   // factor-message-network
   struct FMC; // forward declaration
   
   struct FMC {
     using factor_list = meta::list<>;
     using msg_list = meta::list<>;
     using problem_decomposition = meta::list<>;
   };
};

}

#endif // LP_MP_BINARY_TOMOGRAPHY
