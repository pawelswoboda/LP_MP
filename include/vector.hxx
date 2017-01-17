#ifndef LP_MP_VECTOR_HXX
#define LP_MP_VECTOR_HXX

#include "memory_allocator.hxx"

namespace LP_MP {

// fixed size vector allocated from block allocator
// possibly holding size explicitly is not needed: It is held by allocator as well

// primitive expression templates for vector
template<typename E>
class vector_expression {
public:
   double operator[](const INDEX i) const { return static_cast<E const&>(*this)[i]; }
   INDEX size() const { return static_cast<E const&>(*this).size(); }

   E& operator()() { return static_cast<E&>(*this); }
   const E& operator()() const { return static_cast<const E&>(*this); }
};

class vector : public vector_expression<vector> {
public:
  template<typename ITERATOR>
  vector(ITERATOR begin, ITERATOR end)
  {
    const INDEX size = std::distance(begin,end);
    begin_ = global_real_block_allocator.allocate(size);
    assert(begin_ != nullptr);
    end_ = begin_ + size;
    for(auto it=this->begin(); begin!=end; ++begin, ++it) {
      (*it) = *begin;
    }
  }

  vector(const INDEX size) 
  {
    begin_ = global_real_block_allocator.allocate(size);
    assert(begin_ != nullptr);
    end_ = begin_ + size;
  }
  vector(const INDEX size, const REAL value) 
     : vector(size)
  {
     std::fill(begin_,end_,value);
  }
  ~vector() {
     global_real_block_allocator.deallocate(begin_,1);
  }
   vector(const vector& o)  {
      begin_ = global_real_block_allocator.allocate(o.size());
      end_ = begin_ + o.size();
      assert(begin_ != nullptr);
      auto it = begin_;
      for(auto o_it = o.begin(); o_it!=o.end(); ++it, ++o_it) { *it = *o_it; }
   }
   template<typename E>
   void operator=(const vector_expression<E>& o) {
      assert(size() == o.size());
      for(INDEX i=0; i<o.size(); ++i) { 
         (*this)[i] = o[i]; }
   }
   template<typename E>
   void operator-=(const vector_expression<E>& o) {
      assert(size() == o.size());
      for(INDEX i=0; i<o.size(); ++i) { 
         (*this)[i] -= o[i]; } 
   }
   template<typename E>
   void operator+=(const vector_expression<E>& o) {
      assert(size() == o.size());
      for(INDEX i=0; i<o.size(); ++i) { 
         (*this)[i] += o[i]; } 
   }


   // force construction from expression template
   template<typename E>
   vector(vector_expression<E>& v) : vector(v.size())
   {
      for(INDEX i=0; v.size(); ++i) {
         (*this)[i] = v[i];
      }
   }

   INDEX size() const { return end_ - begin_; }

   REAL operator[](const INDEX i) const {
      assert(i<size());
      assert(!std::isnan(begin_[i]));
      return begin_[i];
   }
   REAL& operator[](const INDEX i) {
      assert(i<size());
      return begin_[i];
   }
   REAL* begin() const { return begin_; }
   REAL* end() const { return end_; }

   template<typename ARCHIVE>
   void serialize(ARCHIVE& ar)
   {
      ar( cereal::binary_data( begin_, sizeof(REAL)*size()) );
   }
private:
  REAL* begin_;
  REAL* end_;
};

class matrix : public vector {
public:
   matrix(const INDEX d1, const INDEX d2) : vector(d1*d2), dim2_(d2) {}
   matrix(const INDEX d1, const INDEX d2, const REAL val) : vector(d1*d2), dim2_(d2) {
      std::fill(this->begin(), this->end(), val);
   }
   //matrix(const matrix& o) vector(o) {
   //   assert(o.dim2_ == dim2_);
   //}
   void operator=(const matrix& o) {
      assert(size() == o.size() && o.dim2_ == dim2_);
      for(INDEX i=0; i<o.size(); ++i) { 
         (*this)[i] = o[i]; 
      }
   }
   REAL& operator()(const INDEX x1, const INDEX x2) {
      assert(x1<dim1() && x2<dim2());
      return (*this)[x1*dim2_ + x2]; 
   }
   REAL operator()(const INDEX x1, const INDEX x2) const { return (*this)[x1*dim2_ + x2]; }
   REAL operator()(const INDEX x1, const INDEX x2, const INDEX x3) const { assert(x3 == 0); return (*this)[x1*dim2_ + x2]; } // sometimes we treat a matrix as a tensor with trivial last dimension
   REAL& operator()(const INDEX x1, const INDEX x2, const INDEX x3) { assert(x3 == 0); return (*this)[x1*dim2_ + x2]; } // sometimes we treat a matrix as a tensor with trivial last dimension
   const INDEX dim1() const { return this->size()/dim2_; }
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
   const INDEX dim2_;
};

class tensor3 : public vector {
public:
   tensor3(const INDEX d1, const INDEX d2, const INDEX d3) : vector(d1*d2*d3), dim2_(d2), dim3_(d3) {}
   tensor3(const INDEX d1, const INDEX d2, const INDEX d3, const REAL val) : tensor3(d1,d2,d3) {
      std::fill(this->begin(), this->end(), val);
   }
   tensor3(const tensor3& o) = default;
   void operator=(const tensor3& o) {
      assert(size() == o.size() && o.dim2_ == dim2_ && o.dim3_ == dim3_);
      for(INDEX i=0; i<o.size(); ++i) { 
         (*this)[i] = o[i]; 
      }
   }
   REAL& operator()(const INDEX x1, const INDEX x2, const INDEX x3) { 
      assert(x1<dim1() && x2<dim2() && x3<dim3());
      return (*this)[x1*dim2_*dim3_ + x2*dim3_ + x3]; 
   }
   REAL operator()(const INDEX x1, const INDEX x2, const INDEX x3) const { return (*this)[x1*dim2_*dim3_ + x2*dim3_ + x3]; }
   const INDEX dim1() const { return this->size()/(dim2_*dim3_); }
   const INDEX dim2() const { return dim2_; }
   const INDEX dim3() const { return dim3_; }
protected:
   const INDEX dim2_, dim3_;
};

template<INDEX FIXED_DIM>
class matrix_view_of_tensor : public vector_expression<matrix_view_of_tensor<FIXED_DIM>> {
public:
   matrix_view_of_tensor(tensor3& t, const INDEX fixed_index) : fixed_index_(fixed_index), t_(t) {}
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
   REAL& operator()(const INDEX x1, const INDEX x2) { 
      if(FIXED_DIM==0) {
         return t_(fixed_index_,x1,x2);
      } else if(FIXED_DIM == 1) {
         return t_(x1,fixed_index_,x2);
      } else {
         return t_(x1,x2,fixed_index_);
      }
   }
   REAL operator()(const INDEX x1, const INDEX x2) const { 
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
   tensor3& t_;
};

// primitive expression templates for all the above linear algebraic classes
template<typename T>
struct scaled_vector : public vector_expression<scaled_vector<T>> {
   scaled_vector(const REAL& omega, const T& a) : omega_(omega), a_(a) {}
   const REAL operator[](const INDEX i) const {
      return omega_*a_[i];
   }
   const REAL operator()(const INDEX i, const INDEX j) const {
      return omega_*a_(i,j);
   }
   const REAL operator()(const INDEX i, const INDEX j, const INDEX k) const {
      return omega_*a_(i,j,k);
   }
   INDEX size() const { return a_.size(); }
   INDEX dim1() const { return a_.dim1(); }
   INDEX dim2() const { return a_.dim2(); }
   INDEX dim3() const { return a_.dim3(); }
   private:
   const T& a_;
   const REAL omega_;
};


template<typename T>
scaled_vector<T> const
operator*(const REAL omega, T const& v) {
   return scaled_vector<T>(omega, v);
}

template<typename T>
struct minus_vector : public vector_expression<minus_vector<T>> {
   minus_vector(const T& a) : a_(a) {}
   const REAL operator[](const INDEX i) const {
      return -a_[i];
   }
   const REAL operator()(const INDEX i, const INDEX j) const {
      return -a_(i,j);
   }
   const REAL operator()(const INDEX i, const INDEX j, const INDEX k) const {
      return -a_(i,j,k);
   }
   INDEX size() const { return a_.size(); }
   INDEX dim1() const { return a_.dim1(); }
   INDEX dim2() const { return a_.dim2(); }
   INDEX dim3() const { return a_.dim3(); }
   private:
   const T& a_;
};

template<typename T>
minus_vector<T> const
operator-(T const& v) {
   return minus_vector<T>(v);
}

} // end namespace LP_MP

#endif // LP_MP_VECTOR_HXX
