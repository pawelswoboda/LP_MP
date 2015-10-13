#ifndef LP_MP_MULTIPLEX_MESSAGE_SIMD
#define LP_MP_MULTIPLEX_MESSAGE_SIMD

#include "LP_MP.h"
#include "factors_messages.hxx"
#include "message_loops.hxx"
#include <chrono>

namespace LP_MP {

#define LP_MP_SIMD

// exactly as MultiplexMessage, but with SIMD enabled
// {LEFT|RIGHT}_SIDE_ACTIVE activates with SFINAE message updates for left and right factor resp.
//
// do zrobienia: provide scalar and SIMD reduction operations etc, so that the message loops can get them if needed. 
// This obviates the need for the fixed lambdas for the message loops. Instead provide a static method dispatcher object
template<class LEFT_LOOP_TYPE, class RIGHT_LOOP_TYPE, bool LEFT_SIDE_ACTIVE, bool RIGHT_SIDE_ACTIVE, bool LEFT_SIMD = false, bool RIGHT_SIMD = false>
class MultiplexMessageSIMD
{
public:
   typedef LEFT_LOOP_TYPE LeftLoopType;
   typedef RIGHT_LOOP_TYPE RightLoopType;

   MultiplexMessageSIMD(LEFT_LOOP_TYPE loopLeft, RIGHT_LOOP_TYPE loopRight, const INDEX kL, const INDEX kR)
      : loopLeft_(loopLeft),
      loopRight_(loopRight),
      kL_(kL),
      kR_(kR),
      delta_(loopLeft_.GetDim(0)),
      t_hot({}), t_cold({})
   { 
      assert(kL_ > 0 && kR_ > 0);
      assert(kL_ == 1 && kR_ == 1); // do zrobienia: for now
   }
   ~MultiplexMessageSIMD() {
      //std::cout << "Cold memory access times = " << std::chrono::duration_cast<std::chrono::nanoseconds>(t_cold).count() << "\n";
      //std::cout << " Hot memory access times = " << std::chrono::duration_cast<std::chrono::nanoseconds>(t_hot).count() << "\n";
      //std::cout << "===================================================\n";
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

   // reparametrize left potential for i-th entry of msg
   template<typename G>
   void RepamLeft(G& repamPot, const REAL msg, const INDEX dim)
   { 
      loopLeft_.loop(dim, [&](const INDEX i) { repamPot[i] -= msg; });
   }
   template<typename G>
   void RepamRight(G& repamPot, const REAL msg, const INDEX msg_index) //__attribute__ ((always_inline))
   {
      if(RIGHT_SIMD) {
         for(INDEX i0=0; i0<loopRight_.GetDim(0); ++i0) {
            repamPot[i0 + msg_index*loopRight_.GetDim(0)] += msg;
         }
      } else if(LEFT_SIMD) {
         for(INDEX i1=0; i1<loopRight_.GetDim(1); ++i1) {
            repamPot[msg_index + i1*loopRight_.GetDim(0)] += msg;
         }
      }
   }
   // do batch reparametrization: more efficient, as contiguous memory access possible.
   template<typename RIGHT_FACTOR_TYPE, typename MSG_VAL_ARRAY>
   void RepamRight(RIGHT_FACTOR_TYPE& rightPot, const MSG_VAL_ARRAY& msg_val)
   {
      if(RIGHT_SIMD) {
#ifndef LP_MP_SIMD
         for(INDEX i1=0; i1<loopRight_.GetDim(1); ++i1) {
            for(INDEX i0=0; i0<loopRight_.GetDim(0); ++i0) {
               rightPot[i0 + i1*loopRight_.GetDim(0)] += msg_val[i1];
            }
         }
#endif
#ifdef LP_MP_SIMD
         // SIMD
         INDEX c=0;
         for(INDEX i1=0; i1<loopRight_.GetDim(1); ++i1) {
            assert(msg_val.vectorsCount() == rightPot.vectorsCount(0));
            assert(msg_val.size() == loopRight_.GetDim(0));
            for(INDEX i0=0; i0<msg_val.vectorsCount(); ++i0) {
               //rightPot.vector(i0 + i1*rightPot.vectorsCount(0)) += msg_val[i1];
               rightPot.vector(c) += msg_val[i1];
               ++c;
            }
         }
#endif
      } else if(LEFT_SIMD) {
#ifndef LP_MP_SIMD
         for(INDEX i1=0; i1<loopRight_.GetDim(1); ++i1) {
            for(INDEX i0=0; i0<loopRight_.GetDim(0); ++i0) {
               rightPot[i0 + i1*loopRight_.GetDim(0)] += msg_val[i0];
            }
         }
#endif
#ifdef LP_MP_SIMD
         // SIMD
         INDEX c=0;
         for(INDEX i1=0; i1<loopRight_.GetDim(1); ++i1) {
            assert(msg_val.vectorsCount() == rightPot.vectorsCount(0));
            assert(msg_val.size() == loopRight_.GetDim(0));
            for(INDEX i0=0; i0<msg_val.vectorsCount(); ++i0) {
               //rightPot.vector(i0 + i1*rightPot.vectorsCount(0)) += msg_val.vector(i0);
               rightPot.vector(c) += msg_val.vector(i0);
               ++c;
            }
         }
#endif
      }

      /*
      // explicit SIMD loop, unaligned version
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
      //for(INDEX i=0; i<loopLeft_.GetDim(0); ++i) {
      //   msg[i] += omega*leftPot[i];
      //}

      // message update is inefficient unless performed in batch!
      for(INDEX i=0; i<loopLeft_.GetDim(0); ++i) {
         delta_[i] = omega*leftPot[i];
      }
      msg += delta_;
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
      // exaclty one must be true
      static_assert(!(LEFT_SIMD && RIGHT_SIMD), "");
      static_assert(LEFT_SIMD || RIGHT_SIMD, "");
      if(RIGHT_SIMD) {
#ifndef LP_MP_SIMD
         {
         // measure time in the first run: "cold access"
         auto t_begin_cold = std::chrono::steady_clock::now();
         std::array<INDEX,2> i;
         constexpr INDEX COMMON_IDX = 1;
         for(i[1]=0; i[1]<loopRight_.GetDim(1); ++i[1]) {
            REAL delta = std::numeric_limits<double>::max();

            for(i[0]=0; i[0]<loopRight_.GetDim(0); ++i[0]) {
               delta = std::min(delta, rightPot[i[0] + i[1]*loopRight_.GetDim(0)]);
            }

            msg[ i[1] ] = msg[ i[1] ] - omega*delta;
         }
         auto t_end_cold = std::chrono::steady_clock::now();
         t_cold += t_end_cold - t_begin_cold;
         }

         {
         // measure "hot" access
         auto t_begin_hot = std::chrono::steady_clock::now();
         std::array<INDEX,2> i;
         constexpr INDEX COMMON_IDX = 1;
         for(i[1]=0; i[1]<loopRight_.GetDim(1); ++i[1]) {
            REAL delta = std::numeric_limits<double>::max();

            for(i[0]=0; i[0]<loopRight_.GetDim(0); ++i[0]) {
               delta = std::min(delta, rightPot[i[0] + i[1]*loopRight_.GetDim(0)]);
            }

            msg[ i[1] ] = msg[ i[1] ] - omega*delta;
         }
         auto t_end_hot = std::chrono::steady_clock::now();
         t_hot += t_end_hot - t_begin_hot;
         }
         return;
#endif
#ifdef LP_MP_SIMD
         {
         // measure time in the first run: "cold access"
         //auto t_begin_cold = std::chrono::steady_clock::now();
         //Vc::Memory<REAL_SIMD> delta(loopRight_.GetDim(1));
         INDEX c=0;
         for(INDEX i1=0; i1<loopRight_.GetDim(1); ++i1) {
            REAL_SIMD delta_tmp(std::numeric_limits<REAL>::max());
            for(INDEX i0=0; i0<rightPot.vectorsCount(0); ++i0) {
               //delta_tmp = Vc::min(delta_tmp, rightPot.vector(i0+i1*rightPot.vectorsCount(0)));
               delta_tmp = Vc::min(delta_tmp, rightPot.vector(c));
               ++c;
            }
            delta_[i1] = omega*delta_tmp.min();
         }
         msg -= delta_;
         //auto t_end_cold = std::chrono::steady_clock::now();
         //t_cold += t_end_cold - t_begin_cold;
         return;
         }

         {
         // measure time in the first run: "cold access"
         auto t_begin_hot = std::chrono::steady_clock::now();
         //Vc::Memory<REAL_SIMD> delta(loopRight_.GetDim(1));
         INDEX c=0;
         for(INDEX i1=0; i1<loopRight_.GetDim(1); ++i1) {
            REAL_SIMD delta_tmp(std::numeric_limits<REAL>::max());
            for(INDEX i0=0; i0<rightPot.vectorsCount(0); ++i0) {
               //delta_tmp = Vc::min(delta_tmp, rightPot.vector(i0+i1*rightPot.vectorsCount(0)));
               delta_tmp = Vc::min(delta_tmp, rightPot.vector(c));
               ++c;
            }
            delta_[i1] = omega*delta_tmp.min();
         }
         msg -= delta_;
         auto t_end_hot = std::chrono::steady_clock::now();
         t_hot += t_end_hot - t_begin_hot;
         }
#endif


         /*
         // do zrobienia: make it fit for aligned data
         // explicit SIMD loops
         for(INDEX i1=0; i1<loopRight_.GetDim(1); ++i1) {
            // compute aligned boundaries which contain relevant data
            const INDEX start_index = i1*loopRight_.GetDim(0);
            const INDEX end_index   = start_index + loopRight_.GetDim(0);
            const INDEX first_aligned_idx = AlignmentUpper(start_index);
            const INDEX  last_aligned_idx = AlignmentLower(end_index); 
            // first go over the first unaligned scalars
            // do zrobienia: the vector goes over the boundaries and wraps into the aligned memory to be processed later. This is ok, if start_index - end_index > REAL_SIMD::Size. Check for this
            REAL_SIMD delta_v = rightPot.vector(toVectorIdx(AlignmentLower(start_index)), start_index - AlignmentLower(start_index));
            for(INDEX idx_v = toVectorIdx(first_aligned_idx); idx_v < toVectorIdx(last_aligned_idx); ++idx_v) {
               delta_v = Vc::min(delta_v, rightPot.vector(idx_v));
            }
            // finally, go over the remaining unaligned scalars
            delta_v = Vc::min(rightPot.vector(toVectorIdx(last_aligned_idx-1), end_index - AlignmentLower(end_index)), delta_v); 
            const REAL delta_s = delta_v.min();
            msg[ i1 ] = msg[ i1 ] - omega*delta_s;
         }
         */

      } else if(LEFT_SIMD) {
#ifndef LP_MP_SIMD
         {
         // measure time in the first run: "cold access"
         auto t_begin_cold = std::chrono::steady_clock::now();
         std::array<INDEX,2> i;
         constexpr INDEX COMMON_IDX = 0;
         for(i[0]=0; i[0]<loopRight_.GetDim(0); ++i[0]) {
            REAL delta = std::numeric_limits<double>::max();

            for(i[1]=0; i[1]<loopRight_.GetDim(1); ++i[1]) {
               delta = std::min(delta, rightPot[i[0] + i[1]*loopRight_.GetDim(0)]);
            }

            msg[ i[0] ] = msg[ i[0] ] - omega*delta;
         }
         auto t_end_cold = std::chrono::steady_clock::now();
         t_cold += t_end_cold - t_begin_cold;
         }

         {
         // measure time in the first run: "cold access"
         auto t_begin_hot = std::chrono::steady_clock::now();
         std::array<INDEX,2> i;
         constexpr INDEX COMMON_IDX = 0;
         for(i[0]=0; i[0]<loopRight_.GetDim(0); ++i[0]) {
            REAL delta = std::numeric_limits<double>::max();

            for(i[1]=0; i[1]<loopRight_.GetDim(1); ++i[1]) {
               delta = std::min(delta, rightPot[i[0] + i[1]*loopRight_.GetDim(0)]);
            }

            msg[ i[0] ] = msg[ i[0] ] - omega*delta;
         }
         auto t_end_hot = std::chrono::steady_clock::now();
         t_hot += t_end_hot - t_begin_hot;
         }
         return;
#endif
#ifdef LP_MP_SIMD
         //Vc::Memory<REAL_SIMD> delta(loopRight_.GetDim(0));
         {
         // measure time in the first run: "cold access"
         //auto t_begin_cold = std::chrono::steady_clock::now();
         // loop for i1=0
         INDEX c=0;
         for(INDEX i0=0; i0<delta_.vectorsCount(); ++i0) {
            delta_.vector(i0) = rightPot.vector(i0);
            ++c;
         }
         for(INDEX i1=1; i1<loopRight_.GetDim(1); ++i1) {
            for(INDEX i0=0; i0<rightPot.vectorsCount(0); ++i0) {
               //delta_.vector(i0) = Vc::min(delta_.vector(i0), rightPot.vector(i0+i1*rightPot.vectorsCount(0)));
               delta_.vector(i0) = Vc::min(delta_.vector(i0), rightPot.vector(c));
               ++c;
            }
         }
         msg -= delta_;
         //auto t_end_cold = std::chrono::steady_clock::now();
         //t_cold += t_end_cold - t_begin_cold;
         return;
         }

         {
         // measure time in the first run: "cold access"
         auto t_begin_hot = std::chrono::steady_clock::now();
         // loop for i1=0
         INDEX c=0;
         for(INDEX i0=0; i0<delta_.vectorsCount(); ++i0) {
            delta_.vector(i0) = rightPot.vector(i0);
            ++c;
         }
         for(INDEX i1=1; i1<loopRight_.GetDim(1); ++i1) {
            for(INDEX i0=0; i0<rightPot.vectorsCount(0); ++i0) {
               //delta_.vector(i0) = Vc::min(delta_.vector(i0), rightPot.vector(i0+i1*rightPot.vectorsCount(0)));
               delta_.vector(i0) = Vc::min(delta_.vector(i0), rightPot.vector(c));
               ++c;
            }
         }
         msg -= delta_;
         auto t_end_hot = std::chrono::steady_clock::now();
         t_hot += t_end_hot - t_begin_hot;
         }
#endif


         /*
         // use aligned version, where we store the message temporarily, possibly get some static storage for this
         // get memory of length dim0 into which we record the maxima
         // is already better than naive version
         std::valarray<double> delta(std::numeric_limits<double>::max(),loopRight_.GetDim(0));
         for(INDEX i1=0; i1<loopRight_.GetDim(1); ++i1) {
            for(INDEX i0=0; i0<loopRight_.GetDim(0); ++i0) {
               const REAL delta_tmp = delta[i0];
               const REAL pot_tmp = rightPot[i0 + i1*loopRight_.GetDim(0)];
               delta[i0] = std::min(delta_tmp, pot_tmp);
               //delta[i0] = std::min(delta[i0], rightPot[i0 + i1*loopRight_.GetDim(0)]);
            }
         }
         msg -= delta; // batch message update might be more efficient than individual updates due to strides
         return;
         */


         
         /*

         for(i0=0; i0 + REAL_SIMD::Size <= loopRight_.GetDim(0); i0+=REAL_SIMD::Size) {
            REAL_SIMD delta_v(std::numeric_limits<double>::max());
            for(int i1=0; i1<loopRight_.GetDim(1); ++i1) {
               INDEX vector_no = toVectorIdx(i0 + i1*loopRight_.GetDim(0));
               delta_v = Vc::min(delta_v, rightPot.vector(vector_no));
            }
            for(int i0_offset=0; i0_offset<REAL_SIMD::Size; i0_offset++) {
               msg[ i0 + i0_offset ] = msg[ i0 + i0_offset ] - delta_v[i0_offset];
            }
         }
         i0 += REAL_SIMD::Size;
         // treat remaining rows extra
         // do zrobienia: do this also with SIMD and some overlap
         for(; i0<loopRight_.GetDim(0); ++i0) {
            REAL delta = std::numeric_limits<double>::max();
            for(int i1=0; i1<loopRight_.GetDim(1); ++i1) {
               delta = std::min(delta, rightPot[ i0 + i1*loopRight_.GetDim(0) ]);
            }
            msg[ i0 ] = msg[ i0 ] - delta; 
         }
         */
      }
   }

   LEFT_LOOP_TYPE loopLeft_;
   RIGHT_LOOP_TYPE loopRight_;
   // do zrobienia: replace this with templates or get from factor
   const INDEX kL_, kR_; // multipliers, see explanation of the model above

   Vc::Memory<REAL_SIMD> delta_; // for storing intermediate msg changes

   std::chrono::duration<size_t, std::nano> t_hot, t_cold;
};

} // end namespace LP_MP

#endif // LP_MP_MULTIPLEX_MESSAGE_SIMD

