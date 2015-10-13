#ifndef LP_MP_MULTIPLEX_MESSAGE_SIMD
#define LP_MP_MULTIPLEX_MESSAGE_SIMD

#include "LP_MP.h"
#include "factors_messages.hxx"
#include "message_loops.hxx"

namespace LP_MP {


// exactly as MultiplexMessage, but with SIMD enabled
// {LEFT|RIGHT}_SIDE_ACTIVE activates with SFINAE message updates for left and right factor resp.
//
// do zrobienia: provide scalar and SIMD reduction operations etc, so that the message loops can get them if needed. 
// This obviates the need for the fixed lambdas for the message loops. Instead provide a static method dispatcher object
template<class LEFT_LOOP_TYPE, class RIGHT_LOOP_TYPE, bool LEFT_SIDE_ACTIVE, bool RIGHT_SIDE_ACTIVE>
class MultiplexMessageiSIMD
{
public:
   typedef LEFT_LOOP_TYPE LeftLoopType;
   typedef RIGHT_LOOP_TYPE RightLoopType;

   MultiplexMessageSIMD(LEFT_LOOP_TYPE loopLeft, RIGHT_LOOP_TYPE loopRight, const INDEX kL, const INDEX kR)
      : loopLeft_(loopLeft),
      loopRight_(loopRight),
      kL_(kL),
      kR_(kR)
   { 
      assert(kL_ > 0 && kR_ > 0);
      assert(kL_ == 1 && kR_ == 1); // do zrobienia: for now
   }
   ~MultiplexMessageSIMD() {
      //std::cout << ddd_tmp << "\n";
   }

   // standard functions which take all possible arguments, to be replaced with minimal ones
   template<typename RIGHT_FACTOR, typename G1, typename G2, bool LSA = LEFT_SIDE_ACTIVE>
   typename std::enable_if<LSA,void>::type
   ReceiveMessageFromRight(RIGHT_FACTOR* const r, const G1& rightPot, G2& msg) 
   {
      MinimizeRight(r,rightPot,msg); 
   }
   template<typename LEFT_FACTOR, typename G1, typename G2, bool RSA = RIGHT_SIDE_ACTIVE>
   typename std::enable_if<RSA,void>::type 
   ReceiveMessageFromLeft(LEFT_FACTOR* l, const G1& leftPot, G2& msg) 
   { 
      MaximizeLeft(l,leftPot,msg); 
   }
   template<typename LEFT_FACTOR, typename RIGHT_FACTOR, typename G1, typename G2, typename G3, bool LSA = LEFT_SIDE_ACTIVE>
   typename std::enable_if<LSA,void>::type
   SendMessageToRight(LEFT_FACTOR* const l, RIGHT_FACTOR* const r, const G1& leftPot, const G2& rightPot, G3& msg, const REAL omega)
   { 
      MaximizeLeft(l,leftPot,msg,omega); 
   }
   template<typename LEFT_FACTOR, typename RIGHT_FACTOR, typename G1, typename G2, typename G3, bool RSA = RIGHT_SIDE_ACTIVE>
   typename std::enable_if<RSA,void>::type
   SendMessageToLeft(LEFT_FACTOR* const l, RIGHT_FACTOR* const r, const G1& leftPot, const G2& rightPot, G3& msg, const REAL omega) 
   { 
      MinimizeRight(r,rightPot,msg,omega); 
   }

   // do zrobienia: also give a possibility to do this with SIMD, i.e.
   // template<typename G>
   // void RepamLeft(G& repamPot, const REAL_SIMD msg, cons INDEX_SIMD dim) {} etc.
   
   // reparametrize left potential for i-th entry of msg
   template<typename G>
   void RepamLeft(G& repamPot, const REAL msg, const INDEX dim)
   { 
      loopLeft_.loop(dim, [&](const INDEX i) { repamPot[i] = repamPot[i] - msg; });
   }
   template<typename G>
   void RepamRight(G& repamPot, const REAL msg, const INDEX msg_index) __attribute__ ((always_inline))
   {
      loopRight_.loop(msg_index, [&](const INDEX i) { repamPot[i] +=  msg; }); return; // possibly is already vectorized

      // explicit SIMD loop
      const INDEX start_index = msg_index*loopRight_.GetDim(0);
      const INDEX end_index   = start_index + loopRight_.GetDim(0);
      const INDEX first_aligned_idx = AlignmentUpper(start_index);
      const INDEX  last_aligned_idx = AlignmentLower(end_index); 

      REAL_SIMD msg_v(msg);

      // first go over the first unaligned scalars
      // build a mask with the overlapping scalars on and add the message conditionally on that mask
      // do the same at the end
      {
      const INDEX begin_scalar_idx = AlignmentLower(start_index);
      const INDEX offset = start_index - begin_scalar_idx;
      REAL_MASK_SIMD mask(INDEX_SIMD::IndexesFromZero() >= offset);

      const INDEX begin_vector_idx = toVectorIdx(begin_scalar_idx);
      //REAL_SIMD begin_v = repamPot.vector(begin_vector_idx);
      //begin_v(mask) += msg_v;
      repamPot.vector(begin_vector_idx) += msg_v;//  = begin_v;
      }

      msg_v = REAL_SIMD(msg);

      for(INDEX idx_v = toVectorIdx(first_aligned_idx); idx_v < toVectorIdx(last_aligned_idx); ++idx_v) {
         repamPot.vector(idx_v) += msg_v;
      }
      
      /*
      {
      const INDEX end_scalar_idx = AlignmentLower(end_index);
      const INDEX offset = end_index - end_scalar_idx;
      const REAL_MASK_SIMD mask = Vc::int_v::IndexesFromZero() < offset;

      const INDEX end_vector_idx = toVectorIdx(end_scalar_idx);
      REAL_SIMD end_v = repamPot.vector(end_vector_idx);
      end_v(mask) += msg_v;
      repamPot.vector(end_vector_idx) = end_v;
      }
      */
   }
   
private:

   template<typename LEFT_FACTOR, typename G1, typename G2>
   void MaximizeLeft(LEFT_FACTOR* const l, const G1& leftPot, G2& msg, const REAL omega = 1.0)
   {
      REAL delta;
      loopLeft_.loop( 
         [&](const INDEX outer_idx){ 
         delta = std::numeric_limits<REAL>::max(); 
         }, 
         [&](const INDEX full_idx, const INDEX outer_idx){ 
         delta = std::min(delta, leftPot[ full_idx ]);
         },
         [&](const INDEX outer_idx){ 
         msg[ outer_idx ] = msg[ outer_idx ] + omega*delta;
         } );

      //Optimize(
      //      [](const REAL y) -> REAL { return y; },
      //      leftPot, msg, l, kL_, loopLeft_, omega);
   }

   // get lower and upper alignment boundaries
   const INDEX AlignmentLower(const INDEX i) {
      assert(i>= 0);
      return i - i % REAL_SIMD::Size;
   }
   const INDEX AlignmentUpper(const INDEX i) {
      assert(i>= 0);
      return i - i % REAL_SIMD::Size + REAL_SIMD::Size;
   }
   const INDEX toVectorIdx(const INDEX idx_s) {
      assert(idx_s % REAL_SIMD::Size == 0);
      return idx_s / REAL_SIMD::Size;
   }
   // do zrobienia: write instructions to directly get lower and upper vector indices

   template<typename RIGHT_FACTOR, typename G1, typename G2>
   void MinimizeRight(RIGHT_FACTOR* r, const G1& rightPot, G2& msg, const REAL omega = 1.0) 
   {
      /*
      std::array<INDEX,2> i;
      constexpr INDEX COMMON_IDX = 1;
      for(i[COMMON_IDX]=0; i[COMMON_IDX]<loopRight_.GetDim(COMMON_IDX); ++i[COMMON_IDX]) {
         REAL delta = std::numeric_limits<double>::max();

         for(i[1-COMMON_IDX]=0; i[1-COMMON_IDX]<loopRight_.GetDim(1-COMMON_IDX); ++i[1-COMMON_IDX]) {
            delta = std::min(delta, rightPot[i[0] + i[1]*loopRight_.GetDim(0)]);
         }
         
         msg[ i[COMMON_IDX] ] = msg[ i[COMMON_IDX] ] - delta;
         //ddd_tmp += delta;
      }
      return;
      */


      // explicit SIMD loops

      // do it directly now and only allow pairwise right factor
      for(INDEX i1=0; i1<loopRight_.GetDim(1); ++i1) {
         // compute aligned boundaries which contain relevant data
         const INDEX start_index = i1*loopRight_.GetDim(0);
         const INDEX end_index   = start_index + loopRight_.GetDim(0);
         const INDEX first_aligned_idx = AlignmentUpper(start_index);
         const INDEX  last_aligned_idx = AlignmentLower(end_index); 
         // first go over the first unaligned scalars
         REAL_SIMD delta_v(std::numeric_limits<REAL>::max());
         // do zrobienia: the vector goes over the boundaries and wraps into the aligned memory to be processed later. This is ok, if start_index - end_index > REAL_SIMD::Size
         delta_v = Vc::min(rightPot.vector(toVectorIdx(AlignmentLower(start_index)), start_index - AlignmentLower(start_index)), delta_v); 
         for(INDEX idx_v = toVectorIdx(first_aligned_idx); idx_v < toVectorIdx(last_aligned_idx); ++idx_v) {
            delta_v = Vc::min(delta_v, rightPot.vector(idx_v));
         }
         // finally, go over the remaining unaligned scalars
         delta_v = Vc::min(rightPot.vector(toVectorIdx(last_aligned_idx-1), end_index - AlignmentLower(end_index)), delta_v); 
         // explicit method:
         //for(INDEX idx_s = last_aligned_idx; idx_s < end_index; ++idx_s) {
         //   delta_s = std::min(delta_s, rightPot[idx_s]);
         //}
         REAL delta_s = delta_v.min();
         //ddd_tmp += delta_s;
         msg[ i1 ] = msg[ i1 ] - omega*delta_s;
      }



      //Optimize(
      //      [](const REAL y) -> REAL { return -y; }, 
      //      rightPot, msg, r, kR_, loopRight_, omega);
   }

   LEFT_LOOP_TYPE loopLeft_;
   RIGHT_LOOP_TYPE loopRight_;
   // do zrobienia: replace this with templates or get from factor
   const INDEX kL_, kR_; // multipliers, see explanation of the model above

   double ddd_tmp;
};

/*
template<class LEFT_FACTOR_TYPE, class RIGHT_FACTOR_TYPE, bool LEFT_SIDE_ACTIVE, bool RIGHT_SIDE_ACTIVE>
//template<typename MultiplexMessageSIMD<LEFT_FACTOR_TYPE,RIGHT_FACTOR_TYPE>::OptimizeOp Op, class FACTOR_TYPE, typename G1, typename G2, class LOOP>
template<typename OP, class FACTOR_TYPE, typename G1, typename G2, class LOOP>
void 
MultiplexMessageSIMD<LEFT_FACTOR_TYPE,RIGHT_FACTOR_TYPE,LEFT_SIDE_ACTIVE,RIGHT_SIDE_ACTIVE>::Optimize(OP op, const G1& pot, G2& msg, FACTOR_TYPE* fac, const INDEX k, LOOP& l, const REAL omega)
//MultiplexMessageSIMD::Optimize(const G1& pot, G2& msg, FACTOR_TYPE* fac, const INDEX k, LOOP& l, const REAL omega)
{
   //const REAL breakPoINDEX = fac->GetBreakpoINDEXCost(pot);
   assert(k == 1); // do zrobienia: templatize for k such that it may be a constexpr, thus fast, possibly get k directly from factor

   REAL delta;
   l.loop( 
         [&](const INDEX outer_idx){ 
         delta = std::numeric_limits<REAL>::max(); 
         }, 
         [&](const INDEX full_idx, const INDEX outer_idx){ 
         delta = std::min(delta, pot[ full_idx ]);
         //delta = std::min(delta, (pot[ full_idx ] - breakPoINDEX)/REAL(k)); // this is very slow, makes algorithm run two times longer
         },
         [&](const INDEX outer_idx){ 
         assert( outer_idx < msg.size() );
         msg[ outer_idx ] = msg[ outer_idx ] + omega*op(delta);
         } );
}
*/

} // end namespace LP_MP

#endif // LP_MP_MULTIPLEX_MESSAGE_SIMD


