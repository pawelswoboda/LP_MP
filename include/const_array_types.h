#ifndef LP_MP_CONST_ARRAY_TYPES
#define LP_MP_CONST_ARRAY_TYPES

#include "instances.inc"

namespace LP_MP {
// used for multiplex_factor to efficientyl hold constant var_capacity and sum for known values

class const_ones_array {
public:
   constexpr const_ones_array(const std::vector<INDEX>& c) {} // assert that c = {1,...,1}
   constexpr const INDEX operator[](const INDEX i) const { return 1; }
};

// instroduce other constant classes. std::integral_type better?
class const_one {
public:
   constexpr const_one(const INDEX i) {} // assert(i == 1)
   constexpr operator INDEX() const { return 1; }
};

}

#endif // LP_MP_CONST_ARRAY_TYPES
