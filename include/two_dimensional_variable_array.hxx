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
   template<typename T2>
   two_dim_variable_array(const two_dim_variable_array<T2>& data)
   {
     allocate_memory(data.size_begin(), data.size_end());
     initialize(data);
   }
   template<typename T2>
   void initialize(const two_dim_variable_array<T2>& data) {
      for(INDEX i=0; i<data.size(); ++i) {
         for(INDEX j=0; j<data[i].size(); ++j) {
            (*this)[i][j] = data[i][j];
         };
      }
   }
   // do zrobienia: remove this constructor
   two_dim_variable_array(const std::vector<INDEX>& dimensions) 
   {
      allocate_memory(dimensions.begin(), dimensions.end());
   }
   // iterator holds size of each dimension of the two dimensional array
   template<typename ITERATOR>
   two_dim_variable_array(ITERATOR begin, ITERATOR end)
   {
      allocate_memory(begin,end);
   }
   two_dim_variable_array() 
     :dim1_(0),
     p_(nullptr)
   {}
   two_dim_variable_array<T>& operator=(const two_dim_variable_array& o)
   {
     if(p_ != nullptr) { free(p_); }
     allocate_memory(o.size_begin(), o.size_end());
     initialize(o);
     return *this;
   }
   two_dim_variable_array(const two_dim_variable_array& o)
   {
     if(p_ != nullptr) {
       free(p_);
     }
     allocate_memory(o.size_begin(), o.size_end());
     initialize(o);
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

   friend std::ostream& operator<<(std::ostream& os, const two_dim_variable_array<T>& a) {
     for(INDEX i=0; i<a.size(); ++i) {
       for(INDEX j=0; j<a[i].size(); ++j) {
         os << a(i,j) << " ";
       }
       os << "\n";
     }
   }

   template<typename ITERATOR>
   void resize(ITERATOR begin, ITERATOR end)
   {
      if(p_ != nullptr) { free(p_); }
      allocate_memory(begin,end);
   }

   struct ArrayAccessObject
   {
      //ArrayAccessObject(T** p) : p_(p) {}
      ArrayAccessObject(T* begin, T* end) : begin_(begin), end_(end) { assert(begin <= end); }
      template<typename VEC>
      void operator=(const VEC& o) 
      { 
        assert(o.size() == this->size());
        const auto s = this->size();
        for(INDEX i=0; i<s; ++i) {
          (*this)[i] = o[i];
        }
      }
      T operator[](const INDEX i) const { assert(i < size()); return begin_[i]; } // possibly do not implement this function but return reference, for copying elements may be costly for sizeof(T) big. Look up std::vector
      T& operator[](const INDEX i) { assert(i < size()); return begin_[i]; }
      INDEX size() const {  return (end_ - begin_); }

      T* begin() { return begin_; }
      T* end() { return end_; }
      auto rbegin() { return std::make_reverse_iterator(end()); }
      auto rend() { return std::make_reverse_iterator(begin()); }

      T const* begin() const { return begin_; }
      T const* end() const { return end_; }
      auto rbegin() const { return std::make_reverse_iterator(end()); }
      auto rend() const { return std::make_reverse_iterator(begin()); }

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
   const T& operator()(const INDEX i, const INDEX j) const
   {
      assert(i < size() && j< (*this)[i].size());
      return *(p_[i] + j);
   }
   T& operator()(const INDEX i, const INDEX j)
   {
      assert(i < size() && j< (*this)[i].size());
      return *(p_[i] + j);
   }

   INDEX size() const { return dim1_; }

   struct iterator : public std::iterator< std::random_access_iterator_tag, T* > {
     iterator(T** x) : x_(x) {}
     void operator++() { ++x_; }
     void operator--() { --x_; }
     iterator& operator+=(const INDEX i) { x_+=i; return *this; }
     iterator& operator-=(const INDEX i) { x_-=i; return *this; }
     iterator operator+(const INDEX i) { iterator it({x_ + i}); return it; }
     iterator operator-(const INDEX i) { iterator it({x_ - i}); return it; }
     auto operator-(const iterator it) const { return x_ - it.x_; }
     ArrayAccessObject operator*() { return ArrayAccessObject(*x_,*(x_+1)); }
     const ArrayAccessObject operator*() const { return ArrayAccessObject(*x_,*(x_+1)); }
     bool operator==(const iterator it) const { return x_ == it.x_; }
     bool operator!=(const iterator it) const { return !(*this == it); }
     T** x_; // pointer to current
   };

   struct reverse_iterator : public std::iterator_traits< T* > {
     reverse_iterator(T** x) : x_(x) {}
     void operator++() { --x_; }
     void operator--() { ++x_; }
     reverse_iterator& operator+=(const INDEX i) { x_-=i; return *this; }
     reverse_iterator& operator-=(const INDEX i) { x_+=i; return *this; }
     reverse_iterator operator+(const INDEX i) { reverse_iterator it({x_ - i}); return it; }
     reverse_iterator operator-(const INDEX i) { reverse_iterator it({x_ + i}); return it; }
     auto operator-(const reverse_iterator it) const { return it.x_ - x_; }
     ArrayAccessObject operator*() { return ArrayAccessObject(*(x_-1),*x_); }
     const ArrayAccessObject operator*() const { return ArrayAccessObject(*(x_-1),*x_); }
     bool operator==(const reverse_iterator it) const { return x_ == it.x_; }
     bool operator!=(const reverse_iterator it) const { return !(*this == it); }
     T** x_; // pointer to current
   };

   iterator begin() { return iterator(p_); }
   iterator end() { return iterator(p_+dim1_); }

   reverse_iterator rbegin() { return reverse_iterator(p_+dim1_); }
   reverse_iterator rend() { return reverse_iterator(p_); }

   struct size_iterator : public std::iterator< std::random_access_iterator_tag, std::size_t > {
     size_iterator(T** x) : x_(x) {}
     void operator++() { ++x_; }
     size_iterator& operator+=(const INDEX i) { x_+=i; return *this; }
     size_iterator operator+(const INDEX i) { size_iterator it({x_ + i}); return it; }
     size_iterator operator-(const INDEX i) { size_iterator it({x_ - i}); return it; }
     const INDEX operator-(const size_iterator it) { return x_ - it.x_; }
     std::size_t operator*() { return *(x_+1) - *x_; }
     T** x_; // pointer to current
   };

   size_iterator size_begin() const { return size_iterator(p_); }
   size_iterator size_end() const { return size_iterator(p_+dim1_); }

private:
   template<typename ITERATOR>
   void allocate_memory(ITERATOR begin, ITERATOR end)
   {
      // first calculate amount of memory needed in bytes
      dim1_ = std::distance(begin, end);
      INDEX neededMem = (dim1_+1)*sizeof(T*);
      for(INDEX i=0; i<dim1_; ++i) {
        assert(*(begin+i) >= 0);
        neededMem += *(begin+i)*sizeof(T);
      }
      p_ = static_cast<T**>(malloc(neededMem)); // what about new ... ?
      if(p_ == nullptr) throw std::runtime_error("Not enough memory");
      T* data_pointer = reinterpret_cast<T*>(&p_[dim1_+1]); // pointer to where data starts
      for(INDEX i=0; i<dim1_; ++i) {
         p_[i] = data_pointer;
         data_pointer += *(begin+i);
      }
      p_[dim1_] = data_pointer;
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
