#ifndef LP_MP_TWO_DIMENSIONAL_VARIABLE_ARRAY_HXX
#define LP_MP_TWO_DIMENSIONAL_VARIABLE_ARRAY_HXX

#include <vector>
#include <cassert>

namespace LP_MP {

// general two-dimensional array with variable first and second dimension sizes, i.e. like vector<vector<T>>. Holds all data contiguously and therefore may be more efficient than vector<vector<T>>

// do zrobienia: - alignment?
//               - custom memory allocator?
//               - iterators
template<typename T>
class two_dim_variable_array
{
public:
   two_dim_variable_array() {}

/*
   template<typename T2>
   two_dim_variable_array(const two_dim_variable_array<T2>& data)
   {
     allocate_memory(data.size_begin(), data.size_end());
     if constexpr(std::is_same<T,T2>::value) {
	     initialize(data);
     }
   }
   */

   template<typename T2>
   void initialize(const two_dim_variable_array<T2>& data) {
      for(std::size_t i=0; i<data.size(); ++i) {
         for(std::size_t j=0; j<data[i].size(); ++j) {
            (*this)[i][j] = data[i][j];
         };
      }
   }

   // iterator holds size of each dimension of the two dimensional array
   template<typename I>
   two_dim_variable_array(const std::vector<I>& size)
   {
      allocate_memory(size.begin(), size.end());
   }
   // iterator holds size of each dimension of the two dimensional array
   template<typename ITERATOR>
   two_dim_variable_array(ITERATOR begin, ITERATOR end)
   {
      allocate_memory(begin,end);
   }

   /*
   friend std::ostream& operator<<(std::ostream& os, const two_dim_variable_array<T>& a) {
     for(std::size_t i=0; i<a.size(); ++i) {
       for(std::size_t j=0; j<a[i].size(); ++j) {
         os << a(i,j) << " ";
       }
       os << "\n";
     }
     return os;
   }
   */

   template<typename ITERATOR>
   void resize(ITERATOR begin, ITERATOR end)
   {
      allocate_memory(begin,end);
   }

   struct ConstArrayAccessObject
   {
      ConstArrayAccessObject(const T* begin, const T* end) : begin_(begin), end_(end) { assert(begin <= end); }
      const T& operator[](const std::size_t i) const { assert(i < size()); return begin_[i]; }
      std::size_t size() const {  return (end_ - begin_); }

      const T* begin() { return begin_; }
      const T* end() { return end_; }
      auto rbegin() { return std::make_reverse_iterator(end()); }
      auto rend() { return std::make_reverse_iterator(begin()); }

      const T* begin() const { return begin_; }
      const T* end() const { return end_; }
      auto rbegin() const { return std::make_reverse_iterator(end()); }
      auto rend() const { return std::make_reverse_iterator(begin()); }

      private:
         const T* begin_;
         const T* end_;
   };
   struct ArrayAccessObject
   {
      ArrayAccessObject(T* begin, T* end) : begin_(begin), end_(end) { assert(begin <= end); }
      template<typename VEC>
      void operator=(const VEC& o) 
      { 
        assert(o.size() == this->size());
        const auto s = this->size();
        for(std::size_t i=0; i<s; ++i) {
          (*this)[i] = o[i];
        }
      }
      const T& operator[](const std::size_t i) const { assert(i < size()); return begin_[i]; }
      T& operator[](const std::size_t i) { assert(i < size()); return begin_[i]; }
      std::size_t size() const {  return (end_ - begin_); }

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
   ArrayAccessObject operator[](const std::size_t i) {
      assert(i<this->size());
      return ArrayAccessObject( &data_[offsets_[i]], &data_[offsets_[i+1]] );
   }
   ConstArrayAccessObject operator[](const std::size_t i) const {
      assert(i<this->size());
      return ConstArrayAccessObject( &data_[offsets_[i]], &data_[offsets_[i+1]] );
   }
   const T& operator()(const std::size_t i, const std::size_t j) const
   {
      assert(i<size() && j< (*this)[i].size());
      return data_[offsets_[i]+j];
   }
   T& operator()(const std::size_t i, const std::size_t j)
   {
      assert(i < size() && j< (*this)[i].size());
      return data_[offsets_[i]+j];
   }

   std::size_t size() const { assert(offsets_.size() > 0); return offsets_.size()-1; }

   struct iterator : public std::iterator< std::random_access_iterator_tag, T* > {
     iterator(T* _data, std::size_t* _offset) : data(_data), offset(_offset) {}
     void operator++() { ++offset; }
     void operator--() { --offset; }
     iterator& operator+=(const std::size_t i) { offset+=i; return *this; }
     iterator& operator-=(const std::size_t i) { offset-=i; return *this; }
     iterator operator+(const std::size_t i) { iterator it(data,offset+i); return it; }
     iterator operator-(const std::size_t i) { iterator it(data,offset-i); return it; }
     auto operator-(const iterator it) const { return offset - it.offset; }
     ArrayAccessObject operator*() { return ArrayAccessObject(data+*offset,data+*(offset+1)); }
     const ArrayAccessObject operator*() const { return ArrayAccessObject(data+*offset,data+*(offset+1)); }
     bool operator==(const iterator it) const { return data == it.data && offset == it.offset; }
     bool operator!=(const iterator it) const { return !(*this == it); }
     T* data;
     std::size_t* offset;
   };

   struct reverse_iterator : public std::iterator_traits< T* > {
     reverse_iterator(T* _data, std::size_t* _offset) : data(_data), offset(_offset) {}
     void operator++() { --offset; }
     void operator--() { ++offset; }
     reverse_iterator& operator+=(const std::size_t i) { offset-=i; return *this; }
     reverse_iterator& operator-=(const std::size_t i) { offset+=i; return *this; }
     iterator operator+(const std::size_t i) { iterator it(data,offset-i); return it; }
     iterator operator-(const std::size_t i) { iterator it(data,offset+i); return it; }
     auto operator-(const reverse_iterator it) const { return it.offset - offset; }
     ArrayAccessObject operator*() { return ArrayAccessObject(data+*(offset-1),data+*offset); }
     const ArrayAccessObject operator*() const { return ArrayAccessObject(data+*(offset-1),data+*offset); }
     bool operator==(const reverse_iterator it) const { return data == it.data && offset == it.offset; }
     bool operator!=(const reverse_iterator it) const { return !(*this == it); }
     T* data;
     std::size_t* offset;
   };

   iterator begin() { return iterator(&data_[0],&offsets_[0]); }
   iterator end() { return iterator(&data_[0],&offsets_.back()); }

   reverse_iterator rbegin() { return reverse_iterator(&data_[0],&offsets_.back()); }
   reverse_iterator rend() { return reverse_iterator(&data_[0],&offsets_[0]); }

   auto size_begin() const { assert(offsets_.size() > 0); return offsets_.begin(); }
   auto size_end() const { assert(offsets_.size() > 0); return offsets_.end()-1; }

   void set(const T& val)
   {
	   for(std::size_t i=0; i<size(); ++i) {
		   for(std::size_t j=0; j<(*this)[i].size(); ++j) {
			   (*this)(i,j) = val;
		   }
	   }
   }

private:
   template<typename ITERATOR>
   void allocate_memory(ITERATOR begin, ITERATOR end)
   {
      // first calculate amount of memory needed in bytes
      const auto s = std::distance(begin, end);
      std::size_t data_size = 0;
      offsets_.reserve(s);
      offsets_.push_back(0);
      for(auto it=begin; it!=end; ++it) {
          offsets_.push_back( offsets_.back() + *it );
          data_size += *it;
      }
      data_.resize(data_size);
   }

   std::vector<std::size_t> offsets_;
   std::vector<T> data_;
   // memory is laid out like this:
   // pointer[1], pointer[2], ... , pointer[dim1], pointer[end], data[1], data[2], ..., data[dim1] .
   // first pointers to the arrays in the second dimension are stored, then contiguously the data itself.
};

} // namespace LP_MP

#endif // LP_MP_TWO_DIMENSIONAL_VARIABLE_ARRAY_HXX
