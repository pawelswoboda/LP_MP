#ifndef LP_MP_CONST_ARRAY_TYPES
#define LP_MP_CONST_ARRAY_TYPES

#include "instances.inc"

namespace LP_MP {
// used for multiplex_factor to efficiently hold constant var_capacity and sum for known values
// do zrobienia: rename constant_types.hxx

// runtime constant array
class constant_array {
public:
   constant_array(const std::vector<INDEX>& c) { 
      for(INDEX i=0; i<c.size()-1; ++i) {
         assert(c[i] == c[i+1]); 
      } 
      if(c.size() > 0) { 
         v_ = c[0]; 
      }
   }
   const INDEX operator[](const INDEX i) const { return v_; }
private:
   INDEX v_;
};

// compile time constant array
template<INDEX N>
class constexpr_array {
public:
   constexpr constexpr_array(const std::vector<INDEX>& c) {}
   constexpr const INDEX operator[](const INDEX i) const { return N; }
};

// do zrobienia: remove this class
class const_ones_array {
public:
   constexpr const_ones_array(const std::vector<INDEX>& c) {}
   constexpr const INDEX operator[](const INDEX i) const { return 1; }
};


template<INDEX N>
class constant {
public:
   constexpr constant(const INDEX i) {}// assert(i == N); }
   constexpr operator INDEX() const { return N; }
};

// do zrobienia: remove this class
class const_one {
public:
   constexpr const_one(const INDEX i) {}
   constexpr operator INDEX() const { return 1; }
};

} // end namespace LP_MP

#endif // LP_MP_CONST_ARRAY_TYPES
