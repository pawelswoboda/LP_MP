#ifndef LP_MP_MESSAGE_LOOP_HXX
#define LP_MP_MESSAGE_LOOP_HXX

#include <utility>
#include "LP_MP.h"

namespace LP_MP {
// for use in multiplex_marg_message.hxx and to templatize loops
// do zrobienia: specify loop constraints or loop concept?
// do zrobienia: check if LAMBDA1, LAMBDA2 and LAMBDA3 are lambdas with correct function specification
// do zrobienia: think about switching loop order for better aligning when COMMON_IDX is true

// loop must implement
// for(i=1:outer) {
//    f1(i);
//    for(j=1:inner) {
//       f2(j + i*outer,i);
//    }
//    f3(i);
// }

template<class ITER_LIMIT = std::array<INDEX,1> >
class UnaryLoop
{
public:
   // default constructor for derived classes with constant limits
   UnaryLoop() {}
   UnaryLoop(const INDEX dim) : dim_({dim}) {}

   template<typename Lambda1, typename Lambda2, typename Lambda3> 
   inline void loop(const Lambda1 f1, const Lambda2 f2, const Lambda3 f3) {
      for(INDEX idx=0; idx<dim_[0]; ++idx) {
         f1(idx);
         f2(idx, idx);
         f3(idx);
      }
   }

   // iterate over all dimensions where msg[msg_index] works on and apply function f on them
   template<typename Lambda>
   inline void loop(const INDEX msg_index, const Lambda f) {
      f(msg_index);
   }

   const INDEX GetDim(const INDEX i) const { assert(i==0); return dim_[0]; }

   const INDEX GetMsgIndex(const INDEX potIndex) const { return potIndex; } // for dimension i of potential, which message dimension acts on that point?
private:
   const ITER_LIMIT dim_;
};

// do zrobienia: write optimized loops for fixed ITER_LIMIT
// pairwise loops for SIMD with VC
// do zrobienia: possibly derive from PairwiseLoop<0>
template<class ITER_LIMIT = std::array<INDEX,2>>
class LeftPairwiseLoopSIMD
{
public:
   LeftPairwiseLoopSIMD(ITER_LIMIT iter_limit) : iter_limit_( iter_limit ) {}
   LeftPairwiseLoopSIMD() {}
   ~LeftPairwiseLoopSIMD() {}

   template<typename Lambda1, typename Lambda2, typename Lambda3> 
   inline void loop(const Lambda1 f1, const Lambda2 f2, const Lambda3 f3) 
   {
      /*
      for(INDEX_SIMD i0 = INDEX_SIMD::IndexesFromZero();
            !(i0<iter_limit_[0]).isEmpty();
            i0 += REAL_SIMD::Size) {
         f1(i0);

         for(INDEX_SIMD i1 = 0; i1<iter_limit_[1]; ++i1) { // do zrobienia: possibly better with INDEX?
            f2( i0 + i1*iter_limit_[0], i0 );
         }
            
         f3(i0);
      }
      */
      // do zrobienia: last i0 indices with ordinary loop
   }

   template<typename Lambda>
   inline void loop(const INDEX msg_index, const Lambda f) {
      for(INDEX i1=0; i1<iter_limit_[1]; ++i1) {
         const INDEX index = msg_index + i1*iter_limit_[0];
         f(index);
      }
   }

   const INDEX GetDim(const INDEX i) const { assert(i==0 || i==1); return iter_limit_[i]; }

private:
   ITER_LIMIT iter_limit_;
};


template<class ITER_LIMIT = std::array<INDEX,2>>
class RightPairwiseLoopSIMD
{
public:
   RightPairwiseLoopSIMD(ITER_LIMIT iter_limit) : iter_limit_( iter_limit ) {}
   RightPairwiseLoopSIMD() {}
   ~RightPairwiseLoopSIMD() {}

   template<typename Lambda1, typename Lambda2, typename Lambda3> 
   inline void loop(const Lambda1 f1, const Lambda2 f2, const Lambda3 f3) 
   {
      //for(INDEX_SIMD i1 = 0; i1<iter_limit_[1]; ++i1) {
      for(INDEX i1 = 0; i1<iter_limit_[1]; ++i1) {
         f1(i1);

         // possibly do that four times and obtain four messages
         //for(INDEX_SIMD i0 = INDEX_SIMD::IndexesFromZero();
         //      !(i0<iter_limit_[0]).isEmpty();
         //      i0 += REAL_SIMD::Size) {

         // attention: this refers to vectors!
         // first go over the first unaligned scalars
         // go over the aligned scalars
         // finally, go over the remaining unaligned scalars

         f3(i1);
      }
   }

   template<typename Lambda>
   inline void loop(const INDEX msg_index, const Lambda f) {
      for(INDEX i0=0; i0<iter_limit_[0]; ++i0) {
         const INDEX index = i0 + msg_index*iter_limit_[0];
         f(index);
      }
   }

   const INDEX GetDim(const INDEX i) const { assert(i==0 || i==1); return iter_limit_[i]; }

private:
   ITER_LIMIT iter_limit_;
};



// left loop corresponds to COMMON_IDX = 0, right one to COMMON_IDX = 1
// do zrobienia: rename common_idx INDEX message_index and full_index INDEXo potential index
template<INDEX COMMON_IDX, class ITER_LIMIT = std::array<INDEX,2> >
class PairwiseLoop
{
public:
   constexpr static INDEX commonIdx = COMMON_IDX;

   //PairwiseLoop(const INDEX first_dim, const INDEX second_dim) : iter_limit_( {{first_dim, second_dim}} ) { }
   PairwiseLoop(ITER_LIMIT iter_limit) : iter_limit_( iter_limit ) { }
   PairwiseLoop() {} // e.g. for non-type template for constant limits -> loop unrolling etc.
   ~PairwiseLoop() {
      static_assert(0 <= COMMON_IDX && COMMON_IDX <= 1, "first template argument must be 0/1");
   }

   template<typename Lambda1, typename Lambda2, typename Lambda3> 
   inline void loop(const Lambda1 f1, const Lambda2 f2, const Lambda3 f3) {
      std::array<INDEX,2> i;
      for(i[COMMON_IDX]=0; i[COMMON_IDX]<iter_limit_[COMMON_IDX]; ++i[COMMON_IDX]) {
         f1(i[COMMON_IDX]);

         for(i[1-COMMON_IDX]=0; i[1-COMMON_IDX]<iter_limit_[1-COMMON_IDX]; ++i[1-COMMON_IDX]) {
            f2( Label(i[0],i[1]) , i[COMMON_IDX]);
         }
         
         f3(i[COMMON_IDX]);
      }
   }

   // iterate over all dimensions where msg[msg_index] works on and apply function f on them
   template<typename Lambda>
   inline void loop(const INDEX msg_index, const Lambda f) {
      std::array<INDEX,2> i;
      i[COMMON_IDX] = msg_index;
      for(i[1-COMMON_IDX]=0; i[1-COMMON_IDX]<iter_limit_[1-COMMON_IDX]; ++i[1-COMMON_IDX]) {
         f(Label(i[0],i[1]));
      }
   }

   const INDEX GetDim(const INDEX i) const { assert(i==0 || i==1); return iter_limit_[i]; }

   // for dimension i of potential, which message dimension acts on that point?
   const INDEX GetMsgIndex(const INDEX potIndex) const { 
      if(COMMON_IDX == 0) {
         return potIndex%iter_limit_[0];
      }
      if(COMMON_IDX == 1) {
         return potIndex/iter_limit_[0];
      }
   } 

   const INDEX Label(const INDEX i0, const INDEX i1) const 
   { 
      assert(i0 < iter_limit_[0]);
      assert(i1 < iter_limit_[1]);
      return i0 + i1*iter_limit_[0]; 
   }

   // labelled is corresponds to the loop f1, to_label to f2
   void PropagateLabel(const typename PrimalSolutionStorage::Element labelled, typename PrimalSolutionStorage::Element to_label) const
   {
      // if labelled[i] == false, set to false all corresponding entries in to_label
      std::array<INDEX,2> i;
      for(i[COMMON_IDX]=0; i[COMMON_IDX]<iter_limit_[COMMON_IDX]; ++i[COMMON_IDX]) {
         if( labelled[ i[COMMON_IDX] ] == false) {
            for(i[1-COMMON_IDX]=0; i[1-COMMON_IDX]<iter_limit_[1-COMMON_IDX]; ++i[1-COMMON_IDX]) {
               to_label[ Label(i[0],i[1]) ] = false;
            }
         } 
      }
   }

private:
   ITER_LIMIT iter_limit_; // possibly derive from ITER_LIMIT to enable empty base class optimization
};

template<INDEX COMMON_IDX1, INDEX COMMON_IDX2, class ITER_LIMIT = std::array<INDEX,3> >
class PairwiseTripletLoop
{
public:
   constexpr static INDEX commonIdx1 = COMMON_IDX1;
   constexpr static INDEX commonIdx2 = COMMON_IDX2;
   // the complementary index
   constexpr static INDEX tripletIdx = 3 - commonIdx1 - commonIdx2;

   //PairwiseLoop(const INDEX first_dim, const INDEX second_dim) : iter_limit_( {{first_dim, second_dim}} ) { }
   PairwiseTripletLoop(ITER_LIMIT iter_limit) : iter_limit_( iter_limit ) { }
   PairwiseTripletLoop() {} // e.g. for non-type template for constant limits -> loop unrolling etc.
   ~PairwiseTripletLoop() {
      static_assert(0 <= COMMON_IDX1 && COMMON_IDX1 < COMMON_IDX2 && COMMON_IDX2 <= 2, "first two template arguments must be in {0,1,2}");
      static_assert(tripletIdx != COMMON_IDX1 && tripletIdx != COMMON_IDX2,"");
   }

   template<typename Lambda1, typename Lambda2, typename Lambda3> 
   inline void loop(const Lambda1 f1, const Lambda2 f2, const Lambda3 f3) {
      std::array<INDEX,3> i;
      for(i[COMMON_IDX1]=0; i[COMMON_IDX1]<iter_limit_[COMMON_IDX1]; ++i[COMMON_IDX1]) {
         for(i[COMMON_IDX2]=0; i[COMMON_IDX2]<iter_limit_[COMMON_IDX2]; ++i[COMMON_IDX2]) {
            f1(PairwiseLabel(i));
            for(i[tripletIdx]=0; i[tripletIdx]<iter_limit_[tripletIdx]; ++i[tripletIdx]) {
               f2( TripletLabel(i) , PairwiseLabel(i));
            }
            f3(PairwiseLabel(i));
         }
      }
   }

   // iterate over all dimensions where msg[msg_index] works on and apply function f on them
   template<typename Lambda>
   inline void loop(const INDEX msg_index, const Lambda f) {
      std::array<INDEX,3> i;
      // correct?
      i[COMMON_IDX1] = msg_index%iter_limit_[commonIdx1];
      i[COMMON_IDX2] = msg_index/iter_limit_[commonIdx1];
      for(i[tripletIdx]=0; i[tripletIdx]<iter_limit_[tripletIdx]; ++i[tripletIdx]) {
         // do zrobienia: possibly update value of triplet label with fewer instructions by reusing previous triplet label. Same for other loop type.
         f(TripletLabel(i));
      }
   }

   const INDEX GetDim(const INDEX i) const { assert(i==0 || i==1 || i==2); return iter_limit_[i]; }

   const INDEX GetMsgIndex(const INDEX potIndex) const {  // for dimension i of potential, which message dimension acts on that point?
      assert(false); // test this method!
      const INDEX i1 = potIndex%iter_limit_[0]; 
      const INDEX i2 = (potIndex/iter_limit_[0])%iter_limit_[1]; 
      const INDEX i3 = potIndex/(iter_limit_[0]*iter_limit_[1]);
      if(tripletIdx == 0) {
         return i2*i3;
      }
      if(tripletIdx == 1) {
         return i1*i3;
      }
      if(tripletIdx == 2) {
         return i1*i2;
      }
      assert(false);
   }

   const INDEX TripletLabel(const std::array<INDEX,3>& i) const 
   { 
      assert(i[0] < iter_limit_[0]);
      assert(i[1] < iter_limit_[1]);
      assert(i[1] < iter_limit_[2]);
      return i[0] + i[1]*iter_limit_[0] + i[2]*iter_limit_[0]*iter_limit_[1]; 
   }
   // ignore index not coming from pairwise potential
   const INDEX PairwiseLabel(const std::array<INDEX,3>& i) const
   { 
      assert(i[COMMON_IDX1] < iter_limit_[COMMON_IDX1]);
      assert(i[COMMON_IDX2] < iter_limit_[COMMON_IDX2]);
      return i[COMMON_IDX1] + i[COMMON_IDX2]*iter_limit_[COMMON_IDX1];
   }

   void PropagateLabel(const typename PrimalSolutionStorage::Element labelled, typename PrimalSolutionStorage::Element to_label) const
   {
      std::array<INDEX,3> i;

      for(i[COMMON_IDX1]=0; i[COMMON_IDX1]<iter_limit_[COMMON_IDX1]; ++i[COMMON_IDX1]) {
         for(i[COMMON_IDX2]=0; i[COMMON_IDX2]<iter_limit_[COMMON_IDX2]; ++i[COMMON_IDX2]) {
            if( labelled[ PairwiseLabel(i) ] == false) {
               for(i[tripletIdx]=0; i[tripletIdx]<iter_limit_[tripletIdx]; ++i[tripletIdx]) {
                  to_label[ TripletLabel(i) ] = false;
               }
            }
         }
      }
   }


private:
   ITER_LIMIT iter_limit_; // possibly derive from ITER_LIMIT to enable empty base class optimization
};
/*
template<INDEX ITER_LIMIT_1, INDEX ITER_LIMIT_2>
struct ConstantIterLimit {
   constexpr INDEX operator[](const INDEX i) const noexcept {
      return i == 0 ? ITER_LIMIT_1 : ( i == 1 ? ITER_LIMIT_2 : throw std::logic_error("i must be 0,1") );
   }
};
*/

// class for flexible number of constant iteration limits
// do zrobienia: possibly make this like tuple, i.e. make a templated getter function instead of constexpr operator[]
template <INDEX ... REST> struct ConstantIterLimit {
   constexpr INDEX operator[](const INDEX i) const noexcept {
      // do zrobienia: make this cleaner
      return i!=-1000 ? throw std::logic_error("accessing empty struct") : -1;
   }
};

template<INDEX ITER_LIMIT, INDEX ... REST>
struct ConstantIterLimit<ITER_LIMIT, REST...> : ConstantIterLimit<REST...> {
   ConstantIterLimit() : ConstantIterLimit<REST...>() {}
   constexpr INDEX operator[](const INDEX i) const noexcept {
      return i == 0 ? ITER_LIMIT : ConstantIterLimit<REST...>::operator[](i-1);
   }
};

// classes for constant iteration limits
template<INDEX ITER_LIMIT>
class UnaryLoopConstant : public UnaryLoop<ConstantIterLimit<ITER_LIMIT> >
{
public:
   UnaryLoopConstant(const INDEX dim) {}
};

template<INDEX COMMON_IDX, INDEX ITER_LIMIT_1, INDEX ITER_LIMIT_2>
class PairwiseLoopConstant : public PairwiseLoop<COMMON_IDX, ConstantIterLimit<ITER_LIMIT_1, ITER_LIMIT_2> >
{
public:
   PairwiseLoopConstant(std::array<INDEX,2>& iter_limit) {}
   //using PairwiseLoop<COMMON_IDX, ConstantIterLimit<ITER_LIMIT_1, ITER_LIMIT_2> >::PairwiseLoop;
};



/*
 * needs C++14

template<INDEX N>
struct ConstArray
{
   INDEX data[N];
   constexpr INDEX& operator[](INDEX i){return data[i];}
   constexpr const INDEX& operator[](INDEX i) const {return data[i];}
};

// transform INDEX template pack as follows:
// let a... be a subsequence of 0,1,2,...,N
// compute the complement.
template<INDEX N, INDEX...x>
constexpr auto complementIntegerPack()
{
   constexpr INDEX x_length = sizeof...(x);
   ConstArray<x_length> a = {x...};
   ConstArray<N - x_length> c = {};

   INDEX a_idx = 0;
   INDEX c_idx = 0;
   for(INDEX i=0; i<N; ++i) {
      if( i < a[a_idx] ) {
         c[c_idx] = i;
         ++c_idx;
      } else {
         ++a_idx;
      }
   }
   return c;
}

// class for recursive partial loop generation.
// Generate ITER_SEL.size() loops with limits ITER_LIMITS. Those are picked by ITER_SEL.
// Also compute the indices 
// assume that ITER_SEL is in ascending order and enumerates a subset of dimensions of ITER_LIMITS
template<class FUNC, class ITER_LIMIT, INDEX CUR_DIM, INDEX ITER_SEL, INDEX ... ITER_SEL_REST>
inline void loopRec(FUNC func, ITER_LIMIT& il, INDEX sel_idx, INDEX full_idx)
{
   static_assert(CUR_DIM <= ITER_SEL, "CUR_DIM <= ITER_SEL");

   if(CUR_DIM == ITER_SEL) {
      for(INDEX i=0; i<il[ITER_SEL]; ++i) {
         loopRec<FUNC, ITER_LIMIT,CUR_DIM+1, ITER_SEL_REST...>(func, il, i + sel_idx * il[ITER_SEL], i+ full_idx * il[ITER_SEL]);
      }
   } else {
      loopRec<FUNC, ITER_LIMIT,CUR_DIM+1, ITER_SEL, ITER_SEL_REST ...>(func, il, sel_idx, full_idx * il[ITER_SEL]);
   }
}

template<class FUNC, class ITER_LIMIT, INDEX CUR_DIM, INDEX ... ITER_SEL_REST>
inline void loopRec(FUNC func, ITER_LIMIT&il, INDEX sel_idx, INDEX full_idx)
{
   func(sel_idx, full_idx);
}


template<INDEX DIM, class ITER_LIMIT = std::array<INDEX,DIM>, INDEX ... COMMON_IDX>
class GeneralLoop
{
   public:
   GeneralLoop(ITER_LIMIT iter_limit) : iter_limit_(iter_limit) { }
   GeneralLoop() { } // e.g. for non-type template for constant limits -> loop unrolling etc.
   ~GeneralLoop() {
      static_assert(0 <= DIM, "DIM must be >= 0");
   }

   template<typename Lambda1, typename Lambda2, typename Lambda3> 
   void loop(const Lambda1 f1, const Lambda2 f2, const Lambda3 f3)
      {
         auto func = [&](const INDEX sel_idx, INDEX full_idx) {
            f1(sel_idx);

            auto f_inner = [&]() { 
               constexpr auto common_idx_complement = complementIntegerPack<COMMON_IDX...>(); 
               auto idx = std::make_INDEXeger_sequence<INDEX, 3>();
               //auto idx = std::make_INDEXeger_sequence<INDEX, DIM - sizeof...(COMMON_IDX)>();
               //loopRec<decltype(f3), ITER_LIMIT, 0, common_idx_complement[idx]...>(f3,il,0,0);
            
            };
            //loopRec<delctype(>();
            f3(full_idx);
         };
         loopRec<decltype(func), ITER_LIMIT, 0, COMMON_IDX ...>(func,il,0,0);

      }

private:
   ITER_LIMIT iter_limit_;
};

*/


// general marginalization loop
// Let A be a 0/1-matric such that columns sum to 1.
// Require A*mu_1 = mu_2, where mu_{1,2} are multiplex-factors.
/*
class ProjectionLoop {

private:
   const INDEX leftDim_;
   const INDEX rightDim_;

};
*/

} // end namespace LP_MP

#endif // LP_MP_MESSAGE_LOOP_HXX
