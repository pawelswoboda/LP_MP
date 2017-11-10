#ifndef LP_MP_TWO_DIMENSIONAL_VARIABLE_ARRAY_HXX
#define LP_MP_TWO_DIMENSIONAL_VARIABLE_ARRAY_HXX

#include "config.hxx"
#include <vector>

namespace LP_MP {

// general two-dimensional array with variable first and second dimension sizes, i.e. like vector<vector<T>>. Holds all data contiguously and therefore may be more efficient than vector<vector<T>>

// do zrobienia: - alignment?
//               - custom memory allocator?
//               - iterators
template<typename T>
class two_dim_variable_array
{
public:
  template<typename TWO_DIM_ARRAY>
   two_dim_variable_array(const TWO_DIM_ARRAY& data)
   {
     InitializeFromArray(data);
   }
  template<typename TWO_DIM_ARRAY>
  void InitializeFromArray(const TWO_DIM_ARRAY& data) {

      class ArraySizeObject {
         public:
            ArraySizeObject(const TWO_DIM_ARRAY& data) : data_(data){}
            INDEX size() const { return data_.size(); }
            INDEX operator[](const INDEX i) const { return data_[i].size(); }
         private:
            const TWO_DIM_ARRAY& data_;
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
   two_dim_variable_array(const std::vector<INDEX>& dimensions) 
   {
      AllocateAndInitializeMemory(dimensions);
   }
   two_dim_variable_array() 
     :dim1_(0),
     p_(nullptr)
  {}
   two_dim_variable_array<T>& operator=(const two_dim_variable_array& o)
   {
     if(p_ != nullptr) { free(p_); }
     InitializeFromArray(o);
     return *this;
   }
   two_dim_variable_array(const two_dim_variable_array& o)
   {
     if(p_ != nullptr) {
       free(p_);
     }
     InitializeFromArray(o);
   }
   two_dim_variable_array(two_dim_variable_array<T>&& o)
   {
      dim1_ = o.dim1_;
      p_ = o.p_;
      o.dim1_ = 0;
      o.p_ = nullptr;
   }
   ~two_dim_variable_array()
   {
      if(p_ != nullptr) {
        free(p_);
      }
   }


   struct ArrayAccessObject
   {
      //ArrayAccessObject(T** p) : p_(p) {}
      ArrayAccessObject(T* begin, T* end) : begin_(begin), end_(end) { assert(begin <= end); }
      template<typename VEC>
      void operator=(const VEC& o) 
      { 
        assert(o.size() == this->size());
        for(INDEX i=0; i<this->size(); ++i) {
          (*this)[i] = o[i];
        }
      }
      T operator[](const INDEX i) const { assert(i < size()); return begin_[i]; } // possibly do not implement this function but return reference, for copying elements may be costly for sizeof(T) big. Look up std::vector
      T& operator[](const INDEX i) { assert(i < size()); return begin_[i]; }
      INDEX size() const {  return (end_ - begin_); }
      T* begin() { return begin_; }
      T* end() { return end_; }
      T const* begin() const { return begin_; }
      T const* end() const { return end_; }
      private:
         T* begin_;
         T* end_;
   };
   ArrayAccessObject operator[](const INDEX i) {
      assert(i<dim1_);
      return ArrayAccessObject( p_[i], p_[i+1] ); 
   }
   ArrayAccessObject operator[](const INDEX i) const {
      assert(i<dim1_);
      return ArrayAccessObject( p_[i], p_[i+1] ); 
   }
   INDEX size() const { return dim1_; }

   struct iterator : public std::iterator< std::random_access_iterator_tag, T* > {
     iterator(T** x) : x_(x) {}
     void operator++() { ++x_; }
     iterator& operator+=(const INDEX i) { x_+=i; return *this; }
     iterator operator+(const INDEX i) { iterator it({x_ + i}); return it; }
     iterator operator-(const INDEX i) { iterator it({x_ - i}); return it; }
     const INDEX operator-(const iterator it) { return x_ - it.x_; }
     ArrayAccessObject operator*() { return ArrayAccessObject(*x_,*(x_+1)); }
     T** x_; // pointer to current
   };

   iterator begin() { return iterator(p_); }
   iterator end() { return iterator(p_+dim1_); }



private:
   template<typename VEC_SIZE>
   void AllocateAndInitializeMemory(const VEC_SIZE& s)
   {
      // first calculate amount of memory needed in bytes
      dim1_ = s.size();
      INDEX neededMem = (dim1_+1)*sizeof(T*);
      for(INDEX i=0; i<s.size(); ++i) {
        assert(s[i] >= 0);
        neededMem += s[i]*sizeof(T);
      }
      p_ = static_cast<T**>(malloc(neededMem)); // what about new ... ?
      if(p_ == nullptr) throw std::runtime_error("Not enough memory");
      T* data_pointer = reinterpret_cast<T*>(&p_[s.size()+1]); // pointer to where data starts
      for(INDEX i=0; i<s.size(); ++i) {
         p_[i] = data_pointer;
         data_pointer += s[i];
      }
      p_[s.size()] = data_pointer;
      static_assert(sizeof(char) == 1,"");
      assert(reinterpret_cast<char*>(data_pointer) == reinterpret_cast<char*>(p_) + neededMem); // char is supposed to be one byte
   }

   INDEX dim1_; // how many arrays are there?
   T** p_; // pointer to continuous memory for two-dimensional array
   // memory is laid out like this:
   // pointer[1], pointer[2], ... , pointer[dim1], pointer[end], data[1], data[2], ..., data[dim1] .
   // first pointers to the arrays in the second dimension are stored, then contiguously the data itself.
};

} // namespace LP_MP

#endif // LP_MP_TWO_DIMENSIONAL_VARIABLE_ARRAY_HXX
