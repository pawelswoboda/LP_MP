#ifndef LP_MP_TWO_DIMENSIONAL_VARIABLE_ARRAY_HXX
#define LP_MP_TWO_DIMENSIONAL_VARIABLE_ARRAY_HXX

#include "instances.inc"
#include <vector>

namespace LP_MP {

// general two-dimensional array with variable first and second dimension sizes, i.e. like vector<vector<T>>. Holds all data contiguously and therefore may be more efficient than vector<vector<T>>

// do zrobienia: - alignment?
//               - custom memory allocator?
//               - iterators
template<typename T>
class TwoDimVariableArray
{
public:
   TwoDimVariableArray(const std::vector<std::vector<T>>& data) : dim1_(data.size())
   {
      class ArraySizeObject {
         public:
            ArraySizeObject(const std::vector<std::vector<T>>& data) : data_(data){}
            const INDEX size() const { return data_.size(); }
            struct SizeObject {
               SizeObject(const INDEX size) : size_(size) {}
               const INDEX size() const { return size_; }
               const INDEX size_;
            };
            SizeObject operator[](const INDEX i) const { return SizeObject(data_[i].size()); }
         private:
            const std::vector<std::vector<T>>& data_;
      };
      AllocateAndInitializeMemory( ArraySizeObject(data) );
      //std::cout << data.size() << "=" << this->size() << "\n";
      for(INDEX i=0; i<data.size(); ++i) {
         //std::cout << data[i].size() << "=" << this->operator[](i).size() << "\n";
         for(INDEX j=0; j<data[i].size(); ++j) {
            (*this)[i][j] = data[i][j];
         };
      }
   }
   TwoDimVariableArray(const std::vector<INDEX>& dimensions) : dim1_(dimensions.size())
   {
      AllocateAndInitializeMemory(dimensions);
   }
   ~TwoDimVariableArray()
   {
      assert(p_ != nullptr);
      free(p_);
   }

   struct ArrayAccessObject
   {
      //ArrayAccessObject(T** p) : p_(p) {}
      ArrayAccessObject(T* begin, T* end) : begin_(begin), end_(end) {}
      //const T operator[](const INDEX i) const { return p_[0][i]; } // possibly do not implement this function but return reference, for copying elements may be costly for sizeof(T) big. Look up std::vector
      const T operator[](const INDEX i) const { return begin_[i]; } // possibly do not implement this function but return reference, for copying elements may be costly for sizeof(T) big. Look up std::vector
      //T& operator[](const INDEX i) { return p_[0][i]; }
      T& operator[](const INDEX i) { return begin_[i]; }
      //const INDEX size() const {  return (p_[1] - p_[0]); }
      const INDEX size() const {  return (begin_ - end_); }
      private:
         //T** p_;
         T* begin_;
         T* end_;
   };
   ArrayAccessObject operator[](const INDEX i) {
      assert(i<dim1_);
      return ArrayAccessObject( &p_[i], &p_[i+1] ); 
   }
   const INDEX size() const { return dim1_; }

private:
   template<typename ARRAY_SIZE_OBJECT>
   void AllocateAndInitializeMemory(const ARRAY_SIZE_OBJECT& s)
   {
      // first calculate amount of memory needed in bytes
      INDEX neededMem = (dim1_+1)*sizeof(T*);
      assert(s.size() == dim1_);
      for(INDEX i=0; i<s.size(); ++i) {
         neededMem += s[i].size()*sizeof(T);
      }
      p_ = static_cast<T**>(malloc(neededMem)); // what about new ... ?
      if(p_ == nullptr) throw std::runtime_error("Not enough memory");
      T* data_pointer = reinterpret_cast<T*>(&p_[s.size()+1]); // pointer to where data starts
      for(INDEX i=0; i<s.size(); ++i) {
         p_[i] = data_pointer;
         data_pointer += s[i].size();
      }
      p_[s.size()] = data_pointer;
      static_assert(sizeof(char) == 1,"");
      assert(reinterpret_cast<char*>(data_pointer) == reinterpret_cast<char*>(p_) + neededMem); // char is supposed to be one byte
   }

   const INDEX dim1_; // how many arrays are there?
   T** p_; // pointer to continuous memory for two-dimensional array
   // memory is laid out like this:
   // pointer[1], pointer[2], ... , pointer[dim1], pointer[end], data[1], data[2], ..., data[dim1] .
   // first pointers to the arrays in the second dimension are stored, then contiguously the data itself.
};

} // namespace LP_MP

#endif // LP_MP_TWO_DIMENSIONAL_VARIABLE_ARRAY_HXX
