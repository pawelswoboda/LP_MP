#ifndef LP_MP_FACTOR_STORAGE_CLASSES
#define LP_MP_FACTOR_STORAGE_CLASSES

#include "LP_MP.h"

// classes for holding/not holding data for factors: Sometimes we do not want to hold the original factor, as it will not be used in the algorithm

namespace LP_MP {

template<bool STORE_FACTOR>
class FactorStorage {};

template<>
class FactorStorage<true>
{
public: 
   FactorStorage(const std::vector<REAL>& cost) : cost_(cost) {}
   const INDEX size() const { return cost_.size(); } 
   const REAL operator[](const INDEX i) { return cost_[i]; }
private:
   const std::vector<double> cost_; // possibly use some allocator here
};

template<>
class FactorStorage<false>
{
public:
   FactorStorage(const std::vector<REAL>& cost) {}
   const INDEX size() const { assert(false); throw std::runtime_error("factor not held"); return 0; }
   const REAL operator[](const INDEX i) const { assert(false); throw std::runtime_error("factor not held"); return 0; }
};

} // end namespace LP_MP

#endif // LP_MP_FACTOR_STORAGE_CLASSES
