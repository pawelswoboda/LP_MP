#ifndef LP_MP_SIMPLEX_MARGINALIZATION_MESSAGE_HXX
#define LP_MP_SIMPLEX_MARGINALIZATION_MESSAGE_HXX

#include "LP_MP.h"
#include "factors_messages.hxx"
#include "vector.hxx"
#include "memory_allocator.hxx"
#include <cmath>
#include "config.hxx"
#ifdef WITH_SAT
#include "sat_solver.hxx"
#endif

namespace LP_MP {

// specialized messages between UnarySimplexFactor and PairwiseSimplexFactor
template<Chirality CHIRALITY, bool SUPPORT_INFINITY = true>
class UnaryPairwiseMessage {
public:
   UnaryPairwiseMessage(const INDEX i1, const INDEX i2) {} // obsolete
   UnaryPairwiseMessage() {} 

   static constexpr INDEX pairwise_index_ = CHIRALITY == Chirality::left ? 0 : 1;

   template<typename RIGHT_FACTOR, typename G2>
   void ReceiveRestrictedMessageFromRight(const RIGHT_FACTOR& r, G2& msg) 
   {
      // we assume that only r.right_primal was assigned, r.left_primal not
      //assert(r.primal_[0] == i1_);
      if(r.primal()[0] < r.dim(0) && r.primal()[1] >= r.dim(1)) {
         vector<REAL> msgs(r.dim(pairwise_index_));
         for(INDEX x=0; x<r.dim(pairwise_index_); ++x) {
            msgs[x] = CHIRALITY == Chirality::left ? r(x,r.primal()[1]) : r(r.primal()[0],x);
         }
         msg -= msgs;
      }
   }

   template<typename G>
   void RepamLeft(G& r, const REAL msg, const INDEX msg_dim)
   {
      assert(!std::isnan(msg));
      if(SUPPORT_INFINITY) {
         r[msg_dim] += normalize( msg );
      } else {
         r[msg_dim] += msg;
      }
      assert(!std::isnan(r[msg_dim]));
   }
   template<typename A1, typename A2>
   void RepamRight(A1& r, const A2& msgs)
   {
      for(INDEX x=0; x<r.dim(pairwise_index_); ++x) {
         assert(!std::isnan(msgs[x]));
         const REAL val = SUPPORT_INFINITY ? normalize(msgs[x]) : msgs[x];
         if(CHIRALITY == Chirality::left) {
            r.msg1(x) += val;
            assert(!std::isnan(r.msg1(x)));
         } else {
            r.msg2(x) += val;
            assert(!std::isnan(r.msg2(x)));
         }
      }
   }
   template<typename G>
   void RepamRight(G& r, const REAL msg, const INDEX dim)
   {
      assert(!std::isnan(msg));
      const REAL val = SUPPORT_INFINITY ? normalize(msg) : msg;
      if(CHIRALITY == Chirality::left) {
         r.msg1(dim) += val;
         assert(!std::isnan(r.msg1(dim)));
      } else {
         r.msg2(dim) += val;
         assert(!std::isnan(r.msg2(dim)));
      }
   }

   template<typename LEFT_FACTOR, typename RIGHT_FACTOR>
   bool ComputeRightFromLeftPrimal(const LEFT_FACTOR& l, RIGHT_FACTOR& r)
   {
      //assert(l.primal() < l.size());
      if(l.primal() < l.size()) {
         const bool changed = (l.primal() != r.primal()[pairwise_index_]);
         r.primal()[pairwise_index_] = l.primal();
         return changed;
      } else {
         return false;
      }
   }

   template<typename LEFT_FACTOR, typename RIGHT_FACTOR>
   bool ComputeLeftFromRightPrimal(LEFT_FACTOR& l, const RIGHT_FACTOR& r)
   {
      //assert(r.primal()[pairwise_index_] < l.size());
      if(r.primal()[pairwise_index_] < l.size()) {
         const bool changed = (l.primal() != r.primal()[pairwise_index_]);
         l.primal() = r.primal()[pairwise_index_];
         return changed; 
      } else {
         return false;
      }
   }

#ifdef WITH_SAT
    template<typename SAT_SOLVER, typename LEFT_FACTOR, typename RIGHT_FACTOR>
    void construct_sat_clauses(SAT_SOLVER& s, const LEFT_FACTOR& l, const RIGHT_FACTOR& r, const sat_var left_begin, const sat_var right_begin) const
    {
       sat_literal_vector left_literals(l.size());
       load_sat_literals(left_begin, left_literals);

       sat_literal_vector right_literals_left(r.dim1());
       sat_literal_vector right_literals_right(r.dim2());
       load_sat_literals(right_begin, right_literals_left, right_literals_right);

       if(CHIRALITY == Chirality::left) {
          s.make_equal(left_literals.begin(), left_literals.end(), right_literals_left.begin(), right_literals_left.end());
       } else {
          assert(CHIRALITY == Chirality::right);
          s.make_equal(left_literals.begin(), left_literals.end(), right_literals_right.begin(), right_literals_right.end());
       }
    }
#endif 

    template<typename LEFT_FACTOR, typename G2>
    void send_message_to_right(const LEFT_FACTOR& l, G2& msg, const REAL omega = 1.0)
    {
      for(INDEX x=0; x<l.size(); ++x) {
        msg[x] -= omega*l[x];
      }
    }
    template<typename RIGHT_FACTOR, typename G2>
    void send_message_to_left(const RIGHT_FACTOR& r, G2& msg, const REAL omega = 1.0)
    {
#ifndef NDEBUG
       const REAL before_lb = r.LowerBound();
#endif
       vector<REAL> msgs(r.dim(pairwise_index_),std::numeric_limits<REAL>::infinity());
       if(CHIRALITY == Chirality::left) {
          r.min_marginal_1(msgs);
       } else {
          r.min_marginal_2(msgs); 
       }
       msg -= omega*msgs;
#ifndef NDEBUG
       const REAL after_lb = r.LowerBound();
       //assert(before_lb <= after_lb); 
#endif
    }

   template<typename LEFT_FACTOR, typename RIGHT_FACTOR>
   bool CheckPrimalConsistency(const LEFT_FACTOR& l, const RIGHT_FACTOR& r) const
   {
      return l.primal() == r.primal()[pairwise_index_];
   } 
};

template<INDEX I1, INDEX I2>
class PairwiseTripletMessage {
public:
   PairwiseTripletMessage() {} 
   PairwiseTripletMessage(const INDEX, const INDEX, const INDEX) {} 
   ~PairwiseTripletMessage() {
      static_assert(I1 < I2 && I2 < 3,""); 
   } 

   // for primal computation as in TRW-S, we need to compute restricted messages as well
   template<typename RIGHT_FACTOR, typename G2>
      void ReceiveRestrictedMessageFromRight(const RIGHT_FACTOR& r, G2& msg, typename PrimalSolutionStorage::Element rightPrimal) 
      {
         throw std::runtime_error("rounding on pairwise factors is not currently supported");
         assert(false);
      }

   template<typename A1, typename A2>
      void RepamLeft(A1& l, const A2& msgs)
      {
         // do zrobienia: possibly use counter
         for(INDEX x1=0; x1<l.dim1(); ++x1) {
            for(INDEX x2=0; x2<l.dim2(); ++x2) {
               l.cost(x1,x2) += normalize( msgs(x1,x2) );
               assert(!std::isnan(l(x1,x2)));
            }
         }
      }
   template<typename A1, typename A2>
      void RepamRight(A1& r, const A2& msgs)
      {
         // do zrobienia: possibly use counter
         if(I1 == 0 && I2 == 1) {
            for(INDEX x1=0; x1<r.dim1(); ++x1) {
               for(INDEX x2=0; x2<r.dim2(); ++x2) {
                  assert(!std::isnan(msgs(x1,x2)));
                  r.msg12(x1,x2) += normalize( msgs(x1,x2) );
                  assert(!std::isnan(r.msg12(x1,x2)));
               }
            }
         } else
            if(I1 == 0 && I2 == 2) {
               for(INDEX x1=0; x1<r.dim1(); ++x1) {
                  for(INDEX x2=0; x2<r.dim3(); ++x2) {
                     assert(!std::isnan(msgs(x1,x2)));
                     r.msg13(x1,x2) += normalize( msgs(x1,x2) );
                     assert(!std::isnan(r.msg13(x1,x2)));
                  }
               }
            } else
               if(I1 == 1 && I2 == 2) {
                  for(INDEX x1=0; x1<r.dim2(); ++x1) {
                     for(INDEX x2=0; x2<r.dim3(); ++x2) {
                        assert(!std::isnan(msgs(x1,x2)));
                        r.msg23(x1,x2) += normalize( msgs(x1,x2) );
                        assert(!std::isnan(r.msg23(x1,x2)));
                     }
                  }
               } else {
                  assert(false);
               }
      }

   template<typename LEFT_FACTOR, typename RIGHT_FACTOR>
      void ComputeRightFromLeftPrimal(const LEFT_FACTOR& l, RIGHT_FACTOR& r)
      {
         if(I1 == 0 && I2 == 1) {

            if(l.primal()[0] < l.dim1()) {
               r.primal()[0] = l.primal()[0];
            }
            if(l.primal()[1] < l.dim2()) {
               r.primal()[1] = l.primal()[1];
            }

         } else
            if(I1 == 0 && I2 == 2) {

               if(l.primal()[0] < l.dim1()) {
                  r.primal()[0] = l.primal()[0];
               }
               if(l.primal()[1] < l.dim2()) {
                  r.primal()[2] = l.primal()[1];
               }

            } else
               if(I1 == 1 && I2 == 2) {

                  if(l.primal()[0] < l.dim1()) {
                     r.primal()[1] = l.primal()[0];
                  }
                  if(l.primal()[1] < l.dim2()) {
                     r.primal()[2] = l.primal()[1];
                  }


               } else {
                  assert(false);
               }
      }

#ifdef WITH_SAT
   template<typename SAT_SOLVER, typename LEFT_FACTOR, typename RIGHT_FACTOR>
      void construct_sat_clauses(SAT_SOLVER& s, const LEFT_FACTOR& l, const RIGHT_FACTOR& r, const sat_var left_begin, const sat_var right_begin) const
      {
         sat_literal_vector left_literals_left(l.dim1());
         sat_literal_vector left_literals_right(l.dim2());
         sat_literal_matrix left_literals_pairwise(l.dim1(), l.dim2());
         load_sat_literals(left_begin, left_literals_left, left_literals_right, left_literals_pairwise);

         sat_literal_matrix right_literals_12(r.dim1(), r.dim2());
         sat_literal_matrix right_literals_13(r.dim1(), r.dim3());
         sat_literal_matrix right_literals_23(r.dim2(), r.dim3());
         load_sat_literals(right_begin, right_literals_12, right_literals_13, right_literals_23);

         if(I1 == 0 && I2 == 1) {
            s.make_equal(left_literals_pairwise.begin(), left_literals_pairwise.end(), right_literals_12.begin(), right_literals_12.end());
         } else if(I1 == 0 && I2 == 2) {
            s.make_equal(left_literals_pairwise.begin(), left_literals_pairwise.end(), right_literals_13.begin(), right_literals_13.end());
         } else if(I1 == 1 && I2 == 2) {
            s.make_equal(left_literals_pairwise.begin(), left_literals_pairwise.end(), right_literals_23.begin(), right_literals_23.end());
         } else {
            assert(false); // not possible
         }
      }
#endif // WITH_SAT

   template<typename LEFT_FACTOR, typename G2>
      void send_message_to_right(const LEFT_FACTOR& l, G2& msg, const REAL omega = 1.0)
      {
         msg -= omega*l;
      }
   template<typename RIGHT_FACTOR, typename G2>
      void send_message_to_left(const RIGHT_FACTOR& r, G2& msg, const REAL omega = 1.0)
      {
         matrix<REAL> msgs(r.dim(I1), r.dim(I2));
         if(I1 == 0 && I2 == 1) {
            r.min_marginal12(msgs);
         } else if(I1 == 0 && I2 == 2) {
               r.min_marginal13(msgs);
         } else if(I1 == 1 && I2 == 2) {
            r.min_marginal23(msgs);
         } else {
            assert(false);
         }
         msg -= omega*msgs;
      }
};

} // end namespace LP_MP

#endif // LP_MP_SIMPLEX_MARGINALIZATION_MESSAGE_HXX
