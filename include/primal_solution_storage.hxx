#ifndef LP_MP_PRIMAL_SOLUTION_STORAGE_HXX
#define LP_MP_PRIMAL_SOLUTION_STORAGE_HXX

#include "config.hxx"
#include <vector>

namespace LP_MP {

// do zrobienia: better use some custom vector with packing. Essentially only 3 values are needed: true, false, unknown, which can be stored with 2 bits. Also using 3 bits for two primal values is possible, but runtime cost may be too high with such a difficult construction. Alternatively, pack 2 entries into 3 bits and then use chunk of, say, 64 bits, such that every k entries we have alignment.
static constexpr unsigned char unknownState = 3;
class PrimalSolutionStorage : public std::vector<char> {
public:
   using Element = typename std::vector<char>::iterator;

   PrimalSolutionStorage() {}

   template<typename FACTOR_ITERATOR>
   PrimalSolutionStorage(FACTOR_ITERATOR factorIt, FACTOR_ITERATOR factorItEnd)
   {
      INDEX size = 0;
      for(;factorIt != factorItEnd; ++factorIt) {
         size += (*factorIt)->PrimalSize();
      }
      this->resize(size);
   }

   void Initialize()
   {
      //std::fill(&p_[0], &p_[0] + size_, 0);
   }

private:
};  
using bit_vector = typename std::vector<bool>::iterator;

// holds a vector of primal classes for each factor type
// what do to with factors that contain primal types with size only known at runtime -> possibly let primal indicate vector of elements, then allocate as many primal elements as needed
/*
template<typename FACTOR_LIST>
class PrimalSolutionStorageTest {
public:
   // iterate over all factors of given type and record how many
   PrimalSolutionStorageTest() {}

   template<typename FACTOR>
   typename std::vector<typename FACTOR::primal>& get()
   {
      constexpr INDEX n = meta::find_index<FACTOR_LIST, FACTOR>::value;
      return std::get<n>(p_);
   }

private:
   struct get_primal_type {
      template<class FACTOR>
         using invoke = std::vector<typename FACTOR::FactorType::primal>;
   };

   using primal_type_list = meta::transform< FACTOR_LIST, get_primal_type >;

   tuple_from_list<primal_type_list> p_;

};
*/

} // end namespace LP_MP

#endif // LP_MP_PRIMAL_SOLUTION_STORAGE_HXX
