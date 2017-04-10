#ifndef LP_MP_VECTOR_HXX
#define LP_MP_VECTOR_HXX

#include "memory_allocator.hxx"
#include "cereal/archives/binary.hpp"

namespace LP_MP {

// fixed size vector allocated from block allocator
// possibly holding size explicitly is not needed: It is held by allocator as well

// primitive expression templates for vector
template<typename T, typename E>
class vector_expression {
public:
   T operator[](const INDEX i) const { return static_cast<E const&>(*this)[i]; }
   T operator()(const INDEX i1, const INDEX i2) const { return static_cast<E const&>(*this)(i1,i2); }
   T operator()(const INDEX i1, const INDEX i2, const INDEX i3) const { return static_cast<E const&>(*this)(i1,i2,i3); }
   INDEX size() const { return static_cast<E const&>(*this).size(); }
   INDEX dim1() const { return static_cast<E const&>(*this).dim1(); }
   INDEX dim2() const { return static_cast<E const&>(*this).dim2(); }
   INDEX dim3() const { return static_cast<E const&>(*this).dim3(); }

   E& operator()() { return static_cast<E&>(*this); }
   const E& operator()() const { return static_cast<const E&>(*this); }
};

template<typename T, typename E>
class matrix_expression {
public:
   T operator[](const INDEX i) const { return static_cast<E const&>(*this)[i]; }
   INDEX size() const { return static_cast<E const&>(*this).size(); }

   T operator()(const INDEX i1, const INDEX i2) const { return static_cast<E const&>(*this)(i1,i2); }
   INDEX dim1() const { return static_cast<E const&>(*this).dim1(); }
   INDEX dim2() const { return static_cast<E const&>(*this).dim2(); }

   E& operator()() { return static_cast<E&>(*this); }
   const E& operator()() const { return static_cast<const E&>(*this); }
};

template<typename T, typename E>
class tensor3_expression {
public:
   T operator[](const INDEX i) const { return static_cast<E const&>(*this)[i]; }
   INDEX size() const { return static_cast<E const&>(*this).size(); }

   T operator()(const INDEX i1, const INDEX i2) const { return static_cast<E const&>(*this)(i1,i2); }
   T operator()(const INDEX i1, const INDEX i2, const INDEX i3) const { return static_cast<E const&>(*this)(i1,i2,i3); }
   INDEX dim1() const { return static_cast<E const&>(*this).dim1(); }
   INDEX dim2() const { return static_cast<E const&>(*this).dim2(); }
   INDEX dim3() const { return static_cast<E const&>(*this).dim3(); }

   E& operator()() { return static_cast<E&>(*this); }
   const E& operator()() const { return static_cast<const E&>(*this); }
};


// possibly also support different allocators: a pure stack allocator without block might be a good choice as well for short-lived memory
template<typename T=REAL>
class vector : public vector_expression<T,vector<T>> {
public:
  template<typename ITERATOR>
  vector(ITERATOR begin, ITERATOR end)
  {
    const INDEX size = std::distance(begin,end);
    assert(size > 0);
    begin_ = global_real_block_allocator_array[stack_allocator_index].allocate(size);
    assert(begin_ != nullptr);
    end_ = begin_ + size;
    for(auto it=this->begin(); begin!=end; ++begin, ++it) {
      (*it) = *begin;
    }
  }

  vector(const INDEX size) 
  {
    begin_ = global_real_block_allocator_array[stack_allocator_index].allocate(size);
    assert(size > 0);
    assert(begin_ != nullptr);
    end_ = begin_ + size;
  }
  vector(const INDEX size, const T value) 
     : vector(size)
  {
     assert(size > 0);
     std::fill(begin_,end_,value);
  }
  ~vector() {
     if(begin_ != nullptr) {
        global_real_block_allocator_array[stack_allocator_index].deallocate(begin_,1);
     }
  }
   vector(const vector& o)  {
      begin_ = global_real_block_allocator_array[stack_allocator_index].allocate(o.size());
      end_ = begin_ + o.size();
      assert(begin_ != nullptr);
      auto it = begin_;
      for(auto o_it = o.begin(); o_it!=o.end(); ++it, ++o_it) { *it = *o_it; }
   }
   vector(vector&& o) {
      begin_ = o.begin_;
      end_ = o.end_;
      o.begin_ = nullptr;
      o.end_ = nullptr;
   }
   template<typename E>
   void operator=(const vector_expression<T,E>& o) {
      assert(size() == o.size());
      for(INDEX i=0; i<o.size(); ++i) { 
         (*this)[i] = o[i]; }
   }
   template<typename E>
   void operator-=(const vector_expression<T,E>& o) {
      assert(size() == o.size());
      for(INDEX i=0; i<o.size(); ++i) { 
         (*this)[i] -= o[i]; } 
   }
   template<typename E>
   void operator+=(const vector_expression<T,E>& o) {
      assert(size() == o.size());
      for(INDEX i=0; i<o.size(); ++i) { 
         (*this)[i] += o[i]; } 
   }


   // force construction from expression template
   template<typename E>
   vector(vector_expression<T,E>& v) : vector(v.size())
   {
      for(INDEX i=0; v.size(); ++i) {
         (*this)[i] = v[i];
      }
   }

   INDEX size() const { return end_ - begin_; }

   T operator[](const INDEX i) const {
      assert(i<size());
      assert(!std::isnan(begin_[i]));
      return begin_[i];
   }
   T& operator[](const INDEX i) {
      assert(i<size());
      return begin_[i];
   }
   using iterator = T*;
   T* begin() const { return begin_; }
   T* end() const { return end_; }

   template<typename ARCHIVE>
   void serialize(ARCHIVE& ar)
   {
      ar( cereal::binary_data( begin_, sizeof(T)*size()) );
   }
private:
  T* begin_;
  T* end_;
};

template<typename T, INDEX N>
using array_impl = std::array<T,N>;



template<typename T, INDEX N>
class array : public vector_expression<T,array<T,N>> {
public:
   array() {}

   template<typename ITERATOR>
   array(ITERATOR begin, ITERATOR end)
   {
      const INDEX size = std::distance(begin,end);
      assert(size == N);
      for(auto it=array_.begin(); it!=array_.end(); ++it) {
         (*it) = *begin;
      }
   }

   array(const T value) 
   {
      std::fill(array_.begin(), array_.end(), value);
   }
   array(const array<T,N>& o)  {
      assert(size() == o.size());
      auto it = begin();
      for(auto o_it = o.begin(); o_it!=o.end(); ++it, ++o_it) { *it = *o_it; }
   }
   template<typename E>
   void operator=(const vector_expression<T,E>& o) {
      assert(size() == o.size());
      for(INDEX i=0; i<o.size(); ++i) { 
         (*this)[i] = o[i]; }
   }
   template<typename E>
   void operator-=(const vector_expression<T,E>& o) {
      assert(size() == o.size());
      for(INDEX i=0; i<o.size(); ++i) { 
         (*this)[i] -= o[i]; } 
   }
   template<typename E>
   void operator+=(const vector_expression<T,E>& o) {
      assert(size() == o.size());
      for(INDEX i=0; i<o.size(); ++i) { 
         (*this)[i] += o[i]; } 
   }


   // force construction from expression template
   template<typename E>
   array(vector_expression<T,E>& v) 
   {
      for(INDEX i=0; v.size(); ++i) {
         (*this)[i] = v[i];
      }
   }

   constexpr static INDEX size() { return N; }

   T operator[](const INDEX i) const {
      assert(i<size());
      assert(!std::isnan(array_[i]));
      return array_[i];
   }
   T& operator[](const INDEX i) {
      assert(i<size());
      return array_[i];
   }
   using iterator = T*;
   auto begin() const { return array_.begin(); }
   auto end() const { return array_.end(); }
   auto begin() { return array_.begin(); }
   auto end() { return array_.end(); }

   template<typename ARCHIVE>
   void serialize(ARCHIVE& ar)
   {
      ar( array_ );
   } 

private:
   std::array<T,N> array_;
};


template<typename T=REAL>
class matrix : public matrix_expression<T,matrix<T>> {
public:
   T& operator[](const INDEX i) { return vec_[i]; }
   const T operator[](const INDEX i) const { return vec_[i]; }
   INDEX size() const { return vec_.size(); }
   template<typename ARCHIVE>
   void serialize(ARCHIVE& ar)
   {
      ar( vec_ );
   }

   T* begin() const { return vec_.begin(); }
   T* end() const { return vec_.end(); }

   matrix(const INDEX d1, const INDEX d2) : vec_(d1*d2), dim2_(d2) {
      assert(d1 > 0 && d2 > 0);
   }
   matrix(const INDEX d1, const INDEX d2, const T val) : vec_(d1*d2), dim2_(d2) {
      assert(d1 > 0 && d2 > 0);
      std::fill(this->begin(), this->end(), val);
   }
   matrix(const matrix& o) 
      : vec_(o.vec_),
      dim2_(o.dim2_) 
   {}
   matrix(matrix&& o) 
      : vec_(std::move(o.vec_)),
      dim2_(o.dim2_) 
   {}
   void operator=(const matrix<T>& o) {
      assert(this->size() == o.size() && o.dim2_ == dim2_);
      for(INDEX i=0; i<o.size(); ++i) { 
         vec_[i] = o[i]; 
      }
   }
   T& operator()(const INDEX x1, const INDEX x2) {
      assert(x1<dim1() && x2<dim2());
      return vec_[x1*dim2_ + x2]; 
   }
   T operator()(const INDEX x1, const INDEX x2) const { return vec_[x1*dim2_ + x2]; }
   T operator()(const INDEX x1, const INDEX x2, const INDEX x3) const { assert(x3 == 0); return vec_[x1*dim2_ + x2]; } // sometimes we treat a matrix as a tensor with trivial last dimension
   T& operator()(const INDEX x1, const INDEX x2, const INDEX x3) { assert(x3 == 0); return vec_[x1*dim2_ + x2]; } // sometimes we treat a matrix as a tensor with trivial last dimension
   const INDEX dim1() const { return vec_.size()/dim2_; }
   const INDEX dim2() const { return dim2_; }

   void transpose() {
      assert(dim1() == dim2());
      for(INDEX x1=0; x1<dim1(); ++x1) {
         for(INDEX x2=0; x2<x1; ++x2) {
            std::swap((*this)(x1,x2), (*this)(x2,x1));
         }
      }
   }
protected:
   vector<T> vec_;
   const INDEX dim2_;
};

template<typename T=REAL>
class tensor3 : public vector<T> { // do zrobienia: remove
public:
   tensor3(const INDEX d1, const INDEX d2, const INDEX d3) : vector<T>(d1*d2*d3), dim2_(d2), dim3_(d3) {}
   tensor3(const INDEX d1, const INDEX d2, const INDEX d3, const T val) : tensor3<T>(d1,d2,d3) {
      std::fill(this->begin(), this->end(), val);
   }
   tensor3(const tensor3<T>& o) :
      vector<T>(o),
      dim2_(o.dim2_),
      dim3_(o.dim3_)
   {}
   tensor3(tensor3&& o) :
      vector<T>(std::move(o)),
      dim2_(o.dim2_),
      dim3_(o.dim3_)
   {}
   void operator=(const tensor3<T>& o) {
      assert(this->size() == o.size() && o.dim2_ == dim2_ && o.dim3_ == dim3_);
      for(INDEX i=0; i<o.size(); ++i) { 
         (*this)[i] = o[i]; 
      }
   }
   T& operator()(const INDEX x1, const INDEX x2, const INDEX x3) { 
      assert(x1<dim1() && x2<dim2() && x3<dim3());
      return (*this)[x1*dim2_*dim3_ + x2*dim3_ + x3]; 
   }
   T operator()(const INDEX x1, const INDEX x2, const INDEX x3) const { return (*this)[x1*dim2_*dim3_ + x2*dim3_ + x3]; }
   const INDEX dim1() const { return this->size()/(dim2_*dim3_); }
   const INDEX dim2() const { return dim2_; }
   const INDEX dim3() const { return dim3_; }
protected:
   const INDEX dim2_, dim3_;
};

template<INDEX FIXED_DIM, typename T=REAL>
class matrix_view_of_tensor : public vector_expression<T,matrix_view_of_tensor<FIXED_DIM,T>> {
public:
   matrix_view_of_tensor(tensor3<T>& t, const INDEX fixed_index) : fixed_index_(fixed_index), t_(t) {}
   ~matrix_view_of_tensor() {
      static_assert(FIXED_DIM < 3,"");
   }
   const INDEX size() const { 
      if(FIXED_DIM==0) {
         return t_.dim2()*t_.dim3();
      } else if(FIXED_DIM == 1) {
         return t_.dim1()*t_.dim2();
      } else {
         return t_.dim2()*t_.dim3();
      }
   }
   const INDEX dim1() const { 
      if(FIXED_DIM==0) {
         return t_.dim2();
      } else {
         return t_.dim1();
      }
   }
   const INDEX dim2() const { 
      if(FIXED_DIM==2) {
         return t_.dim2();
      } else {
         return t_.dim3();
      }
   }
   T& operator()(const INDEX x1, const INDEX x2) { 
      if(FIXED_DIM==0) {
         return t_(fixed_index_,x1,x2);
      } else if(FIXED_DIM == 1) {
         return t_(x1,fixed_index_,x2);
      } else {
         return t_(x1,x2,fixed_index_);
      }
   }
   T operator()(const INDEX x1, const INDEX x2) const { 
      if(FIXED_DIM==0) {
         return t_(fixed_index_,x1,x2);
      } else if(FIXED_DIM == 1) {
         return t_(x1,fixed_index_,x2);
      } else {
         return t_(x1,x2,fixed_index_);
      }
   }

private:
   const INDEX fixed_index_;
   tensor3<T>& t_;
};

// primitive expression templates for all the above linear algebraic classes
template<typename T, typename E>
struct scaled_vector : public vector_expression<T,scaled_vector<T,E>> {
   scaled_vector(const T& omega, const E& a) : omega_(omega), a_(a) {}
   const T operator[](const INDEX i) const {
      return omega_*a_[i];
   }
   const T operator()(const INDEX i, const INDEX j) const {
      return omega_*a_(i,j);
   }
   const T operator()(const INDEX i, const INDEX j, const INDEX k) const {
      return omega_*a_(i,j,k);
   }
   INDEX size() const { return a_.size(); }
   INDEX dim1() const { return a_.dim1(); }
   INDEX dim2() const { return a_.dim2(); }
   INDEX dim3() const { return a_.dim3(); }
   private:
   const T omega_;
   const E& a_;
};


template<typename T, typename E>
scaled_vector<T,vector_expression<T,E>> 
operator*(const T omega, const vector_expression<T,E> & v) {
   return scaled_vector<T,vector_expression<T,E>>(omega, v);
}

template<typename T, typename E>
scaled_vector<T,matrix_expression<T,E>> 
operator*(const T omega, const matrix_expression<T,E> & v) {
   return scaled_vector<T,matrix_expression<T,E>>(omega, v);
}


template<typename T, typename E>
struct minus_vector : public vector_expression<T,minus_vector<T,E>> {
   minus_vector(const E& a) : a_(a) {}
   const T operator[](const INDEX i) const {
      return -a_[i];
   }
   const T operator()(const INDEX i, const INDEX j) const {
      return -a_(i,j);
   }
   const T operator()(const INDEX i, const INDEX j, const INDEX k) const {
      return -a_(i,j,k);
   }
   INDEX size() const { return a_.size(); }
   INDEX dim1() const { return a_.dim1(); }
   INDEX dim2() const { return a_.dim2(); }
   INDEX dim3() const { return a_.dim3(); }
   private:
   const E& a_;
};

template<typename T, typename E>
struct minus_matrix : public matrix_expression<T,minus_matrix<T,E>> {
   minus_matrix(const E& a) : a_(a) {}
   const T operator[](const INDEX i) const {
      return -a_[i];
   }
   const T operator()(const INDEX i, const INDEX j) const {
      return -a_(i,j);
   }
   INDEX size() const { return a_.size(); }
   INDEX dim1() const { return a_.dim1(); }
   INDEX dim2() const { return a_.dim2(); }
   private:
   const E& a_;
};

template<typename T, typename E>
minus_vector<T,vector_expression<T,E>> 
operator-(vector_expression<T,E> const& v) {
   return minus_vector<T,vector_expression<T,E>>(v);
}

template<typename T, typename E>
minus_matrix<T,matrix_expression<T,E>> 
operator-(matrix_expression<T,E> const& v) {
   return minus_matrix<T,matrix_expression<T,E>>(v);
}


} // end namespace LP_MP

#endif // LP_MP_VECTOR_HXX
