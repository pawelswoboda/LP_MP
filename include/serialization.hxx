#ifndef LP_MP_SERIALIZE_HXX
#define LP_MP_SERIALIZE_HXX

#include <array>
#include <vector>
#include "vector.hxx"
#include <cstring>

namespace LP_MP {

  // to do: use perfect forwarding where applicable

template<typename T>
class binary_data {
public:
   binary_data(T* const _p, const INDEX _no_elements) : pointer(_p), no_elements(_no_elements) {}
   T* const pointer;
   const INDEX no_elements; 
};

class allocate_archive {
public:
  // for plain data
  template<typename T>
  static typename std::enable_if<std::is_arithmetic<T>::value,INDEX>::type
  serialize(const T&)
  {
    return sizeof(T); 
  }
  
  // for arrays
  template<typename T>
  static INDEX serialize(const T*, const INDEX s)
  {
    return sizeof(T)*s;
  }
  template<typename T>
  static INDEX serialize( const binary_data<T> b )
  {
     return serialize(b.pointer, b.no_elements); 
  }

  // for std::array<T>
  template<typename T, std::size_t N>
  static INDEX serialize(const std::array<T,N>&)
  {
    return sizeof(T)*N;
  }

  // for vector<T>
  template<typename T>
  static INDEX serialize(const vector<T>& v)
  {
     return serialize(v.begin(), v.size());
  }

  template<typename T>
  static INDEX serialize(const matrix<T>& m)
  {
     return serialize(m.begin(), m.size());
  }

  // for std::vector<T>
  template<typename T>
  static INDEX serialize(const std::vector<T>& v)
  {
     return serialize(v.data(), v.size());
  }

  template<typename... T_REST>
  void operator()(T_REST&&... types)
  {}
  template<typename T, typename... T_REST>
  void operator()(T&& t, T_REST... types)
  {
     size_in_bytes_ += serialize(t);
     (*this)(types...);
  }

  INDEX size() const
  {
    return size_in_bytes_; 
  }

private:
  INDEX size_in_bytes_ = 0;
};

class serialization_archive {
public:
  template<typename ITERATOR, typename SERIALIZATION_FUN>
  serialization_archive(ITERATOR begin, ITERATOR end, SERIALIZATION_FUN serialization_fun)
  {
    allocate_archive s;
    for(auto it=begin; it!=end; ++it) {
      serialization_fun(*it,s);
    }

    const INDEX size_in_bytes = s.size();

    archive_ = new char[size_in_bytes];
    end_ = archive_ + size_in_bytes;
    assert(archive_ != nullptr);
    cur_ = archive_;
  }

  serialization_archive(allocate_archive& a)
  {
     const INDEX size_in_bytes = a.size();

     archive_ = new char[size_in_bytes];
     end_ = archive_ + size_in_bytes;
     assert(archive_ != nullptr);
     cur_ = archive_;
  }

  serialization_archive(void* mem, INDEX size_in_bytes)
  {
     assert(mem != nullptr);
     archive_ = (char*) mem;
     cur_ = archive_;
     end_ = archive_ + size_in_bytes;
  }

  ~serialization_archive()
  {
    assert(archive_ != nullptr);
    if(archive_ != nullptr) {
      delete[] archive_;
    }
  }

  void release_memory()
  { 
     archive_ = nullptr;
     cur_ = nullptr;
     end_ = nullptr;
  }

  bool operator==(serialization_archive& o) const
  {
     if(size() != o.size()) { 
        return false;
     }
     return std::memcmp(archive_, o.archive_, size()) == 0;
  }

  INDEX size() const 
  {
     return end_ - archive_;
  }

  char* cur_address() const
  {
    return cur_;
  }

  void advance(const SIGNED_INDEX bytes) 
  {
    cur_ += bytes;
    assert(cur_ >= archive_);
    assert(cur_ <= end_);
  }

  void reset_cur()
  {
     cur_ = archive_;
  } 

private:
  char* archive_ = nullptr;
  char* end_ = nullptr;
  char* cur_;

};

// write data into archive
class save_archive {
public:
   save_archive(serialization_archive& a) 
      : ar(a) 
   {
      ar.reset_cur();
   }

  // for arrays
  template<typename T>
  void serialize(const T* p, const INDEX size)
  {
      INDEX* s = (INDEX*) ar.cur_address();
      std::memcpy((void*) s, (void*) p, size*sizeof(T));
      const INDEX size_in_bytes = sizeof(T)*size;
      ar.advance(size_in_bytes);
  }
  template<typename T>
  void serialize( const binary_data<T> b )
  {
     serialize(b.pointer, b.no_elements); 
  }

  // for std::array<T,N>
  template<typename T, std::size_t N>
  void serialize(const std::array<T,N>& v)
  {
    T* s = (T*) ar.cur_address();
    for(auto x : v) {
       *s = x;
       ++s; 
    }
    const INDEX size_in_bytes = sizeof(T)*N;
    ar.advance(size_in_bytes);
  } 
 
  // for vector<T>
  template<typename T>
  void serialize(const vector<T>& v)
  {
     serialize(v.begin(), v.size());
  }

  template<typename T>
  void serialize(const matrix<T>& m)
  {
     serialize(m.begin(), m.size());
  }

  // for std::vector<T>
  template<typename T>
  void serialize(const std::vector<T>& v)
  {
     serialize(v.data(), v.size());
  }

   // for plain data
   template<typename T>
   typename std::enable_if<std::is_arithmetic<T>::value>::type
   serialize(T& t)
   {
      T* s = (T*) ar.cur_address();
      *s = t;
      ar.advance(sizeof(T));
   }


   // save multiple entries
  template<typename... T_REST>
  void operator()(T_REST&&... types)
  {}
  template<typename T, typename... T_REST>
  void operator()(T&& t, T_REST&&... types)
  {
     serialize(t);
     (*this)(types...);
  }


private:
  serialization_archive& ar; 
};

// write data from archive into objects
class load_archive
{
public:
   load_archive(serialization_archive& a) 
      : ar(a) 
   {
      ar.reset_cur(); 
   }

   // for arrays
  template<typename T>
  void serialize(T* pointer, const INDEX size)
  {
      T* s = (T*) ar.cur_address();
      std::memcpy((void*) pointer, (void*) s, size*sizeof(T));
      const INDEX size_in_bytes = sizeof(T)*size;
      ar.advance(size_in_bytes);
  }
  template<typename T>
  void serialize( binary_data<T> b )
  {
     serialize(b.pointer, b.no_elements); 
  }

  // for std::array<T,N>
  template<typename T, std::size_t N>
  void serialize(std::array<T,N>& v)
  {
    T* s = (T*) ar.cur_address();
    for(auto& x : v) {
       x = *s;
       ++s; 
    }
    const INDEX size_in_bytes = sizeof(T)*N;
    ar.advance(size_in_bytes);
  } 
 
  // for vector<T>
  template<typename T>
  void serialize(vector<T>& v)
  {
     assert(v.size() == *((INDEX*)ar.cur_address()));
     serialize(v.begin(), v.size());
  }

  template<typename T>
  void serialize(matrix<T>& m)
  {
     assert(m.size() == *((INDEX*)ar.cur_address()));
     serialize(m.begin(), m.size());
  }

  // for std::vector<T>
  template<typename T>
  void serialize(std::vector<T>& v)
  {
     assert(v.size() == *((INDEX*)ar.cur_address()));
     serialize(v.data(), v.size());
  } 

   // for plain data
   template<typename T>
   typename std::enable_if<std::is_arithmetic<T>::value>::type
   serialize(T& t)
   {
      T* s = (T*) ar.cur_address();
      t = *s;
      ar.advance(sizeof(T));
   }

   // save multiple entries
  template<typename... T_REST>
  void operator()(T_REST&&... types)
  {}
  template<typename T, typename... T_REST>
  void operator()(T&& t, T_REST&&... types)
  {
     serialize(t);
     (*this)(types...);
  }

private:
  serialization_archive& ar;
};

// add numeric values stored in archive to variables
template<SIGNED_INDEX PREFIX>
class addition_archive {
public:
   addition_archive(serialization_archive& a) 
      : ar(a) 
   {
      ar.reset_cur(); 
   }

   // for arrays
   template<typename T>
     void serialize(T* pointer, const INDEX size)
     {
       static_assert(std::is_same<T,float>::value || std::is_same<T,double>::value,"");
       T* val = (T*) ar.cur_address();
       for(INDEX i=0; i<size; ++i) {
         pointer[i] += T(PREFIX) * val[i]; 
       }
       const INDEX size_in_bytes = sizeof(T)*size;
       ar.advance(size_in_bytes);
     }
   template<typename T>
     void serialize( binary_data<T> b )
     {
       serialize(b.pointer, b.no_elements); 
     }

   // for std::array<T,N>
   template<typename T, std::size_t N>
     void serialize(std::array<T,N>& v)
     {
       static_assert(std::is_same<T,float>::value || std::is_same<T,double>::value,"");
       T* s = (T*) ar.cur_address();
       for(auto& x : v) {
         x += T(PREFIX) * (*s);
         ++s; 
       }
       const INDEX size_in_bytes = sizeof(T)*N;
       ar.advance(size_in_bytes);
     } 

   // for vector<T>
   template<typename T>
     void serialize(vector<T>& v)
     {
       assert(v.size() == *((INDEX*)ar.cur_address()));
       serialize(v.begin(), v.size());
     }

   template<typename T>
     void serialize(matrix<T>& m)
     {
       assert(m.size() == *((INDEX*)ar.cur_address()));
       serialize(m.begin(), m.size());
     }

   // for std::vector<T>
   template<typename T>
     void serialize(std::vector<T>& v)
     {
       assert(v.size() == *((INDEX*)ar.cur_address()));
       serialize(v.data(), v.size());
     } 

   // for plain data
   template<typename T>
     typename std::enable_if<std::is_arithmetic<T>::value>::type
     serialize(T& t)
     {
       static_assert(std::is_same<T,float>::value || std::is_same<T,double>::value,"");
       T* s = (T*) ar.cur_address();
       t += T(PREFIX) * (*s);
       ar.advance(sizeof(T));
     }

   // save multiple entries
   template<typename... T_REST>
     void operator()(T_REST&&... types)
     {}
   template<typename T, typename... T_REST>
     void operator()(T&& t, T_REST&&... types)
     {
       serialize(t);
       (*this)(types...);
     }

private:
  serialization_archive& ar;
};

} // end namespace LP_MP
#endif // LP_MP_SERIALIZE_HXX

