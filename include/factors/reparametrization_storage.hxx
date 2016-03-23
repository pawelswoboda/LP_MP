#ifndef LP_MP_POTENTIAL_STORAGE_CLASSES_HXX
#define LP_MP_POTENTIAL_STORAGE_CLASSES_HXX

#include "LP_MP.h"
#include "tolerance.hxx"

namespace LP_MP {

// specialized storage classes for reparametrized potentials

template<INDEX N>
class FixedSizeExplicitRepamStorage {
public:
   template<typename FACTOR_CONTAINER>
   class type : public std::array<REAL,N> {
      using ARRAY_TYPE = std::array<REAL,N>;
   public:
      template<typename FACTOR_TYPE, typename ARRAY>
      type(const FACTOR_TYPE& f, const ARRAY& cost)
      {
         for(INDEX i=0; i<cost.size(); ++i) {
            ARRAY_TYPE::operator[](i) = cost[i];
         }
      }

      const REAL operator[](const INDEX i) const {
         assert(i < N);
         return ARRAY_TYPE::operator[](i);
      }
      REAL& operator[](const INDEX i) {
         assert(i < N);
         return ARRAY_TYPE::operator[](i);
      }
      constexpr INDEX size() const { return N; }
   };
};

template<class FACTOR_CONTAINER> 
class ExplicitRepamStorage : public std::vector<REAL>
{
   //using FactorType = typename FACTOR_CONTAINER::FactorType;
   using ARRAY_TYPE = std::vector<REAL>;
public:
   template<typename FACTOR_TYPE, typename ARRAY>
   ExplicitRepamStorage(const FACTOR_TYPE& f, const ARRAY& cost)
   : ARRAY_TYPE(cost.size()) 
   {
      for(INDEX i=0; i<cost.size(); ++i) {
         ARRAY_TYPE::operator[](i) = cost[i];
      }
   }

   const REAL operator[](const INDEX i) const {
      assert(i < this->size());
      /*
      const REAL orig_pot = static_cast<const FACTOR_CONTAINER*>(this)->factor_.operator[](i);
      const REAL test_val = static_cast<const FACTOR_CONTAINER*>(this)->factor_.operator[](i) + static_cast<const FACTOR_CONTAINER*>(this)->GetMessageSum(i);
      const REAL true_val = ARRAY_TYPE::operator[](i);

      assert(std::abs( test_val - true_val) < eps);
      */
      return ARRAY_TYPE::operator[](i);
   }
   REAL& operator[](const INDEX i) {
      assert(i < this->size());
      return ARRAY_TYPE::operator[](i);
   }
   /*
   //void GetRepam(std::vector<REAL>& pot) { pot = repam_; }
   //std::vector<REAL>& GetRepam() { return repam_; }

   //template<typename G>
   //inline void UpdateRepam(const G& d) {
   //   for(INDEX i=0; i<repam_.size(); ++i) repam_[i] += d[i];
   //}
   //void UpdateRepam(const REAL diff, const INDEX dim) {
   //   repam_[dim] += diff;
   //}
   const REAL operator[](const INDEX i) const { 
      assert(i < repam_.size()); 
      return repam_[i]; 
   }
   REAL& operator[](const INDEX i) { assert(i < repam_.size()); return repam_[i]; }

   // possibly remove those, when valarray will be used (or write custom ones?)
   auto begin() -> decltype(ARRAY_TYPE().begin()) { return repam_.begin(); }
   auto end()   -> decltype(ARRAY_TYPE().end()) { return repam_.end(); }
   auto cbegin() -> decltype(ARRAY_TYPE().cbegin()) const { return repam_.cbegin(); }
   auto cend()   -> decltype(ARRAY_TYPE().cend()) const { return repam_.cend(); }
private:
   ARRAY_TYPE repam_; // do zrobienia: change to valarray
   */
};



// possibly introduce traits, so that factor knows whether to hold original potential, based on whether reparametrization storage needs access to it
// reparametrization storage
template<class FACTOR_CONTAINER>
class ExplicitRepamStorageSIMD : public Vc::Memory<REAL_SIMD>
{
public:
   using Memory = Vc::Memory<REAL_SIMD>;
   template<typename ARRAY>
   ExplicitRepamStorageSIMD(const ARRAY& pot) : Vc::Memory<REAL_SIMD>(pot.size())
   {
      for(INDEX i=0; i<Memory::entriesCount(); ++i) {
         Memory::scalar(i) = pot[i];
      }
   }
   const size_t size() const { return Memory::entriesCount(); }
   REAL& operator[](const INDEX i) { return Memory::scalar(i); }
   const REAL operator[](const INDEX i) const { return Memory::scalar(i); }
   /*
   template<typename ARRAY>
   ExplicitRepamStorageSIMD(const ARRAY& pot) :val_(pot.size()) {
      for(INDEX i=0; i<val_.entriesCount(); ++i) {
         val_[i] = pot[i];
      }
   }
   const size_t size() const { return val_.entriesCount(); }
   const REAL operator[](const INDEX i) const {
      const REAL x = val_[i];
      return x;
   }
   REAL& operator[](const INDEX i) {
      REAL& x = val_[i];
      return x;
   }

private:
   Vc::Memory<REAL_SIMD> val_;
   */
};




template<class FACTOR_TYPE_CRTP>
class PairwiseStorageSIMD {
// reparametrized potential is held in padded Vc::Memory
public:
   template<typename FACTOR_TYPE, typename ARRAY>
   PairwiseStorageSIMD(const FACTOR_TYPE& f, const ARRAY& cost)
   :  
      dim_({{f.GetDim(0), f.GetDim(1)}}),
      repam_(NumberOfRequiredVectors(dim_[0]) * REAL_SIMD::Size * dim_[1])
   {
      assert(repam_.entriesCount() >= cost.size());
      //const INDEX v1 = dim_[0]/REAL_SIMD::Size + (dim_[0]%REAL_SIMD::Size == 0 ? 0 : 1) ; // number of vectors for each dimension
      
      // pad with infinity all entries of vectors that are not occupied by actual data. Fill rest with actual potential.
      // dirty hack: first put inf everywhere. Better: only where needed.
      for(INDEX i=0; i<repam_.entriesCount(); ++i) {
         repam_[i] = std::numeric_limits<REAL>::max();
      }
      assert(this->size() == cost.size());
      for(INDEX i=0; i<cost.size(); i++) {
         this->operator[](i) = cost[i];
      }
   }

   const REAL operator[](const INDEX idx) const {
      const INDEX scalar_idx = GetScalarIndex(idx);
      return repam_[scalar_idx];
   }

   REAL& operator[](const INDEX idx) {
      const INDEX scalar_idx = GetScalarIndex(idx);
      return repam_[scalar_idx];
   }
   
   auto vector(const INDEX i) -> decltype(std::declval<Vc::Memory<REAL_SIMD>>().vector(0)) {
      return repam_.vector(i);
   }
   auto vector(const INDEX i) const -> decltype(std::declval<const Vc::Memory<REAL_SIMD>>().vector(0)) {
      return repam_.vector(i);
   }

   //template<INDEX DIM>
   const INDEX vectorsCount() const { 
      return repam_.vectorsCount();
   }
   // number of vectors in dimension 0
   const INDEX vectorsCount(const INDEX i) const {
      if(i==0) {
         return NumberOfRequiredVectors(dim_[0]);
      }
      throw std::runtime_error("not implemented for i == 1");
      //static_assert("not to be called","");
   }
   //template<INDEX DIM>
   const INDEX entriesCount(const INDEX i) const {
      if(i==1) {
         return dim_[1];
      }
      throw std::runtime_error("not implemented for i==0");
      return -1;
      //static_assert("only specialization to be called","");
   }

   const INDEX size() const {
      return dim_[0]*dim_[1];
   }

private:

   // how many vectors do we need for i scalars?
   const INDEX NumberOfRequiredVectors(const INDEX i) const {
      return i/REAL_SIMD::Size + (i%REAL_SIMD::Size == 0 ? 0 : 1);
   }

   using dim = std::array<INDEX,2>;
   dim GetDim(const INDEX i) const { 
      assert(i<dim_[0]*dim_[1]);
      return dim({{i%dim_[0], i/dim_[0]}});
   }

   INDEX GetScalarIndex(const INDEX idx) const
   {
      dim i = GetDim(idx);
      const INDEX v1 = i[0]/REAL_SIMD::Size;
      const INDEX v_offset = i[0]%REAL_SIMD::Size;
      const INDEX v1_max = NumberOfRequiredVectors(dim_[0]);
      const INDEX v = v1 + i[1]*v1_max;
      const INDEX scalar_idx = v*REAL_SIMD::Size + v_offset;
      return scalar_idx;
   }

   const dim dim_; // do zrobienia: do not store size information here: It can be accessed through the inherited factor via CRTP
   Vc::Memory<REAL_SIMD> repam_;
};

template<typename FACTOR_CONTAINER>
class ImplicitRepamStorage
{
   //using FactorType = typename FACTOR_CONTAINER::FactorType;
public:
   // do zrobienia: possibly infer FACTOR_TYPE from FACTOR_CONTAINER, gives one template less
   template<typename FACTOR_TYPE, typename ARRAY>
   ImplicitRepamStorage(const FACTOR_TYPE& f, const ARRAY& cost) {}
   const REAL operator[](const INDEX i) const
   {
      return static_cast<const FACTOR_CONTAINER*>(this)->factor_.operator[](i) + static_cast<const FACTOR_CONTAINER*>(this)->GetMessageSum(i);
   }
   // note: we do not provide write access via operator[], as this allows the message to omit all writes completely
   const INDEX size() const 
   { 
      return static_cast<const FACTOR_CONTAINER*>(this)->factor_.size(); 
   }
};





} // end namespace LP_MP

#endif // LP_MP_POTENTIAL_STORAGE_CLASSES_HXX
