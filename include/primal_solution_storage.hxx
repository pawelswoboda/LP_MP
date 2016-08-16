#ifndef LP_MP_PRIMAL_SOLUTION_STORAGE_HXX
#define LP_MP_PRIMAL_SOLUTION_STORAGE_HXX

#include <vector>

namespace LP_MP {

// do zrobienia: better use some custom vector with packing. Essentially only 3 values are needed: true, false, unknown, which can be stored with 2 bits. Also using 3 bits for two primal values is possible, but runtime cost may be too high with such a difficult construction. Alternatively, pack 2 entries into 3 bits and then use chunk of, say, 64 bits, such that every k entries we have alignment.
static constexpr unsigned char unknownState = 3;
class PrimalSolutionStorage : public std::vector<unsigned char> {
public:
   using Element = std::vector<unsigned char>::iterator;

   PrimalSolutionStorage() {}
   template<typename FACTOR_ITERATOR>
   PrimalSolutionStorage(FACTOR_ITERATOR factorIt, FACTOR_ITERATOR factorItEnd)
   {
      INDEX size = 0;
      for(;factorIt != factorItEnd; ++factorIt) {
         size += (*factorIt)->PrimalSize();
      }
      this->resize(size,unknownState);
      //this->shrink_to_fit();
      //std::cout << "primal size = " << this->size() << "\n";
   }

   void Initialize()
   {
      std::fill(this->begin(), this->end(), unknownState);
   }

};  


} // end namespace LP_MP

#endif // LP_MP_PRIMAL_SOLUTION_STORAGE_HXX
