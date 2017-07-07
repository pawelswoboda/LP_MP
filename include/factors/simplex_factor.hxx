#ifndef LP_MP_SIMPLEX_FACTOR_HXX
#define LP_MP_SIMPLEX_FACTOR_HXX

#include "LP_MP.h"
#include "memory_allocator.hxx"
#include "vector.hxx"
//#include "cereal/types/array.hpp"
#ifdef WITH_SAT
#include "sat_interface.hxx"
#endif

// Investigate contingency tables (Knuth) for more general tabular structures.

namespace LP_MP {

// the polytope {x >= 0 : x_1 + ... + x_n = 1 }
// dual problem:
//      max_{z} z 
//      s.t.    z <= repamPot[i]

class UnarySimplexFactor : public vector<REAL> {
public:
   UnarySimplexFactor(const std::vector<REAL>& cost) : vector<REAL>(cost.begin(), cost.end()) {}
   UnarySimplexFactor(const INDEX n) : vector<REAL>(n, 0.0) {}

   REAL LowerBound() const { 
      const REAL lb = *std::min_element(this->begin(), this->end()); 
      assert(std::isfinite(lb));
      return lb;
   }

   REAL EvaluatePrimal() const 
   { 
      if(primal_ >= size()) {
         return std::numeric_limits<REAL>::infinity();
      }
      return (*this)[primal_]; 
   }
   void MaximizePotentialAndComputePrimal() 
   {
      if(primal_ >= size()) {
         primal_ = std::min_element(this->begin(), this->end()) - this->begin();
         assert(primal_ < size());
      }
   }

   // load/store function for the primal value
   template<class ARCHIVE> void serialize_primal(ARCHIVE& ar) { ar(primal_); }
   template<class ARCHIVE> void serialize_dual(ARCHIVE& ar) { ar( *static_cast<vector<REAL>*>(this) ); }

   void init_primal() { primal_ = size(); }
   INDEX primal() const { return primal_; }
   INDEX& primal() { return primal_; }
   void primal(const INDEX p) { primal_ = p; }

   INDEX subgradient(double* w) const
   {
      assert(primal_ < size());
      std::fill(w, w+this->size(), 0.0);
      w[primal_] = 1.0;
      return this->size();
   }
   REAL dot_product(double* w) const
   {
      return w[primal_];
   }

#ifdef WITH_SAT
   template<typename SAT_SOLVER>
   void construct_sat_clauses(SAT_SOLVER& s) const
   {
      auto vars = create_sat_variables(s, size());
      add_simplex_constraint_sat(s, vars.begin(), vars.end());
   }

   template<typename VEC>
   void reduce_sat(VEC& assumptions, const REAL th, sat_var begin) const
   {
      const REAL lb = LowerBound();
      for(INDEX i=0; i<this->size(); ++i) {
         if((*this)[i] > lb + th) { 
            assumptions.push_back(-to_literal(begin+i));
         }
      } 
   }

   template<typename SAT_SOLVER>
   void convert_primal(SAT_SOLVER& s, sat_var first)
   {
      for(INDEX i=first; i<first+this->size(); ++i) {
         if(lglderef(s,to_literal(i)) == 1) {
            primal_ = i-first;
         }
      }
   }
#endif

private:
   INDEX primal_;
};

// do zrobienia: if pairwise was supplied to us (e.g. external factor, then reflect this in constructor and only allocate space for messages.
// When tightening, we can simply replace pairwise pointer to external factor with an explicit copy. Reallocate left_msg_ and right_msg_ to make memory contiguous? Not sure, depends whether we use block_allocator, which will not acually release the memory
// when factor is copied, then pairwise_ must only be copied if it is actually modified. This depends on whether we execute SMRP or MPLP style message passing. Templatize for this possibility
class PairwiseSimplexFactor : public matrix_expression<REAL, PairwiseSimplexFactor> {
public:
   PairwiseSimplexFactor(const INDEX dim1, const INDEX dim2) : dim1_(dim1), dim2_(dim2)
   {
      const INDEX size = dim1_*dim2_ + dim1_ + dim2_;
      pairwise_ = global_real_block_allocator.allocate(size);
      //pairwise_ = new REAL[size]; // possibly use block allocator!
      assert(pairwise_ != nullptr);
      std::fill(pairwise_, pairwise_ + size, 0.0);
      left_msg_ = pairwise_ + dim1_*dim2_;
      right_msg_ = left_msg_ + dim1_;
   }

   template<typename MATRIX>
   PairwiseSimplexFactor(const INDEX dim1, const INDEX dim2, const MATRIX& m) 
   : PairwiseSimplexFactor(dim1,dim2)
   {
      for(INDEX x1=0; x1<this->dim1(); ++x1) {
         for(INDEX x2=0; x2<this->dim2(); ++x2) {
            this->cost(x1,x2) = m(x1,x2);
         }
      }
   }

   ~PairwiseSimplexFactor() {
      global_real_block_allocator.deallocate(pairwise_,1);
   }
   PairwiseSimplexFactor(const PairwiseSimplexFactor& o) : dim1_(o.dim1_), dim2_(o.dim2_) {
      const INDEX size = dim1_*dim2_ + dim1_ + dim2_;
      pairwise_ = global_real_block_allocator.allocate(size);
      //pairwise_ = new REAL[dim1_*dim2_ + dim1_ + dim2_]; // possibly use block allocator!
      assert(pairwise_ != nullptr);
      left_msg_ = pairwise_ + dim1_*dim2_;
      right_msg_ = left_msg_ + dim1_;
      for(INDEX i=0; i<dim1_*dim2_; ++i) { pairwise_[i] = o.pairwise_[i]; }
      for(INDEX i=0; i<dim1_; ++i) { left_msg_[i] = o.left_msg_[i]; }
      for(INDEX i=0; i<dim2_; ++i) { right_msg_[i] = o.right_msg_[i]; }
   }
   void operator=(const PairwiseSimplexFactor& o) {
      assert(dim1_ == o.dim1_ && dim2_ == o.dim2_);
      for(INDEX i=0; i<dim1_*dim2_; ++i) { pairwise_[i] = o.pairwise_[i]; }
      for(INDEX i=0; i<dim1_; ++i) { left_msg_[i] = o.left_msg_[i]; }
      for(INDEX i=0; i<dim2_; ++i) { right_msg_[i] = o.right_msg_[i]; }
   }

   REAL operator[](const INDEX x) const {
      const INDEX x1 = x/dim2_;
      const INDEX x2 = x%dim2_;
      assert(x1 < dim1_ && x2 < dim2_);
      return pairwise_[x] + left_msg_[x1] + right_msg_[x2];
   }
   // below is not nice: two different values, only differ by const!
   REAL operator()(const INDEX x1, const INDEX x2) const {
      assert(x1 < dim1_ && x2 < dim2_);
      return pairwise_[x1*dim2_ + x2] + left_msg_[x1] + right_msg_[x2];
   }
   REAL& cost(const INDEX x1, const INDEX x2) {
      assert(x1 < dim1_ && x2 < dim2_);
      return pairwise_[x1*dim2_ + x2];
   }
   REAL LowerBound() const {
      REAL lb = std::numeric_limits<REAL>::infinity();
      for(INDEX x1=0; x1<dim1_; ++x1) {
         for(INDEX x2=0; x2<dim2_; ++x2) {
            lb = std::min(lb, (*this)(x1,x2));
         }
      }
      assert(std::isfinite(lb));
      return lb;
   }

   
   INDEX dim1() const { return dim1_; }
   INDEX dim2() const { return dim2_; }
     
   REAL& pairwise(const INDEX x1, const INDEX x2) { return pairwise_[x1*dim2_ + x2]; }
   REAL& msg1(const INDEX x1) { return left_msg_[x1]; }
   REAL& msg2(const INDEX x2) { return right_msg_[x2]; }

   INDEX size() const { return dim1_*dim2_; }

   void init_primal() 
   {
      primal_[0] = dim1();
      primal_[1] = dim2();
   }
   REAL EvaluatePrimal() const { 
      if(primal_[0] >= dim1() || primal_[1] >= dim2()) {
         return std::numeric_limits<REAL>::infinity();
      }
      assert(primal_[0] < dim1());
      assert(primal_[1] < dim2());
      const REAL val = (*this)(primal_[0], primal_[1]); 
      //assert(val < std::numeric_limits<REAL>::infinity());
      return val;
   }
   void MaximizePotentialAndComputePrimal() 
   {
      REAL min_val = std::numeric_limits<REAL>::infinity();
      for(INDEX x1=0; x1<dim1(); ++x1) {
         for(INDEX x2=0; x2<dim2(); ++x2) {
            if(min_val >= (*this)(x1,x2)) {
               min_val = (*this)(x1,x2);
               primal_[0] = x1;
               primal_[1] = x2;
            }
         }
      }
   }

   template<class ARCHIVE> void serialize_primal(ARCHIVE& ar) { ar( primal_[0], primal_[1] ); }
   //template<class ARCHIVE> void serialize_dual(ARCHIVE& ar) { ar( cereal::binary_data( pairwise_, sizeof(REAL)*(size()+dim1()+dim2()) ) ); }
   template<class ARCHIVE> void serialize_dual(ARCHIVE& ar) { ar( binary_data<REAL>( pairwise_, sizeof(REAL)*(size()+dim1()+dim2()) ) ); }

   template<typename VECTOR>
   void min_marginal_1(VECTOR& m) const
   {
      assert(m.size() == dim1());
      std::fill(m.begin(), m.end(), std::numeric_limits<REAL>::infinity());
      for(INDEX x1=0; x1<dim1(); ++x1) {
          for(INDEX x2=0; x2<dim2(); ++x2) {
             m[x1] = std::min(m[x1],(*this)(x1,x2));
          }
       } 
   }

   template<typename VECTOR>
   void min_marginal_2(VECTOR& m) const
   {
      assert(m.size() == dim2());
      std::fill(m.begin(), m.end(), std::numeric_limits<REAL>::infinity());
      for(INDEX x1=0; x1<dim1(); ++x1) {
          for(INDEX x2=0; x2<dim2(); ++x2) {
             m[x2] = std::min(m[x2],(*this)(x1,x2));
          }
       } 
   }

   INDEX subgradient(double* w) const
   {
      assert(primal_[0] < dim1() && primal_[1] < dim2());
      std::fill(w, w+this->size(), 0.0);
      w[primal_[0]*dim2_ + primal_[1]] = 1.0;
      return this->size();
   }
   REAL dot_product(double* w) const
   {
      return w[primal_[0]*dim2_ + primal_[1]];
   }

#ifdef WITH_SAT
   template<typename SAT_SOLVER>
   void construct_sat_clauses(SAT_SOLVER& s) const
   {
      auto vars = create_sat_variables(s, dim1() + dim2() + dim1()*dim2());
      sat_var pairwise_var_begin = vars[dim1() + dim2()];

      std::vector<sat_var> tmp_vars;
      tmp_vars.reserve(std::max(dim1(), dim2()));

      for(INDEX x1=0; x1<dim1(); ++x1) {
         for(INDEX x2=0; x2<dim2(); ++x2) {
            tmp_vars.push_back(pairwise_var_begin + x1*dim2() + x2);
         }
         sat_var c = add_at_most_one_constraint_sat(s, tmp_vars.begin(), tmp_vars.end());
         make_sat_var_equal(s, to_literal(c), to_literal(vars[x1]));
         tmp_vars.clear();
      }
      for(INDEX x2=0; x2<dim2(); ++x2) {
         for(INDEX x1=0; x1<dim1(); ++x1) {
            tmp_vars.push_back(pairwise_var_begin + x1*dim2() + x2);
         }
         sat_var c = add_at_most_one_constraint_sat(s, tmp_vars.begin(), tmp_vars.end());
         make_sat_var_equal(s, to_literal(c), to_literal(vars[dim1() + x2]));
         tmp_vars.clear();
      }

      // is superfluous: summation constraints must be active due to unary simplex factors
      //add_simplex_constraint_sat(sat_var.begin(), sat_var.begin()+dim1());
      //add_simplex_constraint_sat(sat_var.begin()+dim1(), sat_var.begin()+dim1()+dim2());
   }

   template<typename VEC>
   void reduce_sat(VEC& assumptions, const REAL th, sat_var begin) const
   {
      begin += dim1() + dim2();
      const REAL lb = LowerBound();
      for(INDEX x1=0; x1<this->dim1(); ++x1) {
         for(INDEX x2=0; x2<this->dim2(); ++x2) {
            if((*this)(x1,x2) > lb + th) { 
               assumptions.push_back(-to_literal(begin + x1*dim2() + x2));
            }
         }
      } 
   }

   template<typename SAT_SOLVER>
   void convert_primal(SAT_SOLVER& s, sat_var first)
   {
      for(INDEX x1=0; x1<this->dim1(); ++x1) {
         if(lglderef(s,to_literal(first + x1)) == 1) {
            primal_[0] = x1;
         } 
      }
      for(INDEX x2=0; x2<this->dim2(); ++x2) {
         if(lglderef(s,to_literal(first + dim1() + x2)) == 1) {
            primal_[1] = x2;
         } 
      }
   }
#endif

   const std::array<INDEX,2>& primal() const { return primal_; }
   std::array<INDEX,2>& primal() { return primal_; }

   //INDEX primal_[0], primal_[1]; // not so nice: make getters and setters!
private:
   // those three pointers should lie contiguously in memory.
   REAL* pairwise_;
   REAL* left_msg_;
   REAL* right_msg_;
   std::array<INDEX,2> primal_;
   const INDEX dim1_, dim2_;

};

// factor assumes that triplet potentials is empty and holds only messages to pairwise factors, i.e. is latently factorizable
class SimpleTighteningTernarySimplexFactor : public tensor3_expression<REAL, SimpleTighteningTernarySimplexFactor> {
public:
   SimpleTighteningTernarySimplexFactor(const INDEX dim1, const INDEX dim2, const INDEX dim3) : dim1_(dim1), dim2_(dim2), dim3_(dim3) 
   {
      const INDEX size = dim1_*dim2_ + dim1_*dim3_ + dim2_*dim3_;
      
      msg12_ = global_real_block_allocator.allocate(size);
      assert(msg12_ != nullptr);
      msg13_ = msg12_ + dim1_*dim2_;
      msg23_ = msg13_ + dim1_*dim3_;
      std::fill(msg12_, msg12_ + size, 0.0);
   }
   ~SimpleTighteningTernarySimplexFactor() {
      global_real_block_allocator.deallocate(msg12_,1);
   }
   SimpleTighteningTernarySimplexFactor(const SimpleTighteningTernarySimplexFactor& o) : dim1_(o.dim1_), dim2_(o.dim2_), dim3_(o.dim3_) {
      const INDEX size = dim1_*dim2_ + dim1_*dim3_ + dim2_*dim3_;
      
      msg12_ = global_real_block_allocator.allocate(size);
      assert(msg12_ != nullptr);
      msg13_ = msg12_ + dim1_*dim2_;
      msg23_ = msg13_ + dim1_*dim3_;
      for(INDEX i=0; i<dim1_*dim2_; ++i) { msg12_[i] = o.msg12_[i]; }
      for(INDEX i=0; i<dim1_*dim3_; ++i) { msg13_[i] = o.msg13_[i]; }
      for(INDEX i=0; i<dim2_*dim3_; ++i) { msg23_[i] = o.msg23_[i]; }
   }
   void operator=(const SimpleTighteningTernarySimplexFactor& o) {
      assert(dim1_ == o.dim1_ && dim2_ == o.dim2_ && dim3_ == o.dim3_);
      for(INDEX i=0; i<dim1_*dim2_; ++i) { msg12_[i] = o.msg12_[i]; }
      for(INDEX i=0; i<dim1_*dim3_; ++i) { msg13_[i] = o.msg13_[i]; }
      for(INDEX i=0; i<dim2_*dim3_; ++i) { msg23_[i] = o.msg23_[i]; }
   }

   REAL LowerBound() const {
      REAL lb = std::numeric_limits<REAL>::infinity();
      for(INDEX x1=0; x1<dim1_; ++x1) {
         for(INDEX x2=0; x2<dim2_; ++x2) {
            for(INDEX x3=0; x3<dim3_; ++x3) {
               lb = std::min(lb, (*this)(x1,x2,x3));
            }
         }
      }
      assert(std::isfinite(lb));
      return lb;
   }

   REAL EvaluatePrimal() const
   {
      if(primal_[0] >= dim1_ || primal_[1] >= dim2_ || primal_[2] >= dim3_) {
         return std::numeric_limits<REAL>::infinity();
      }
      //assert((*this)(primal_[0], primal_[1], primal_[2]) < std::numeric_limits<REAL>::infinity());
      return (*this)(primal_[0], primal_[1], primal_[2]);
   }

   REAL operator()(const INDEX x1, const INDEX x2, const INDEX x3) const {
      return msg12(x1,x2) + msg13(x1,x3) + msg23(x2,x3);
   }

   /*
   REAL operator[](const INDEX x) const {
      const INDEX x1 = x / (dim2_*dim3_);
      const INDEX x2 = ( x % (dim2_*dim3_) ) / dim3_;
      const INDEX x3 = x % dim3_;
      return msg12(x1,x2) + msg13(x1,x3) + msg23(x2,x3);
   }
   */

   REAL min_marginal12(const INDEX x1, const INDEX x2) const {
      REAL marg = (*this)(x1,x2,0);
      for(INDEX x3=1; x3<dim3_; ++x3) {
         marg = std::min(marg, (*this)(x1,x2,x3));
      }
      return marg;
   }
   template<typename MSG>
   void min_marginal12(MSG& msg) const {
      assert(msg.dim1() == dim1_);
      assert(msg.dim2() == dim2_);
      //for(INDEX x1=0; x1<dim1_; ++x1) {
      //   for(INDEX x2=0; x2<dim2_; ++x2) {
      //      msg(x1,x2) = std::numeric_limits<REAL>::infinity();
      //   }
      //}
      for(INDEX x1=0; x1<dim1_; ++x1) {
         for(INDEX x2=0; x2<dim2_; ++x2) {
            msg(x1,x2) = min_marginal12(x1,x2);
         }
      }
   }

   template<typename MSG>
   void min_marginal13(MSG& msg) const {
      assert(msg.dim1() == dim1_);
      assert(msg.dim2() == dim3_);
      for(INDEX x1=0; x1<dim1_; ++x1) {
         for(INDEX x3=0; x3<dim3_; ++x3) {
            msg(x1,x3) = std::numeric_limits<REAL>::infinity();
         }
      }
      for(INDEX x1=0; x1<dim1_; ++x1) {
         for(INDEX x2=0; x2<dim2_; ++x2) {
            for(INDEX x3=0; x3<dim3_; ++x3) {
               msg(x1,x3) = std::min(msg(x1,x3), (*this)(x1,x2,x3));
            }
         }
      }
   }
   template<typename MSG>
   void min_marginal23(MSG& msg) const {
      assert(msg.dim1() == dim2_);
      assert(msg.dim2() == dim3_);
      for(INDEX x2=0; x2<dim2_; ++x2) {
         for(INDEX x3=0; x3<dim3_; ++x3) {
            msg(x2,x3) = std::numeric_limits<REAL>::infinity();
         }
      }
      for(INDEX x1=0; x1<dim1_; ++x1) {
         for(INDEX x2=0; x2<dim2_; ++x2) {
            for(INDEX x3=0; x3<dim3_; ++x3) {
               msg(x2,x3) = std::min(msg(x2,x3), (*this)(x1,x2,x3));
            }
         }
      }
   }

   INDEX size() const { return dim1()*dim2()*dim3(); }

   INDEX dim1() const { return dim1_; }
   INDEX dim2() const { return dim2_; }
   INDEX dim3() const { return dim3_; }

   REAL msg12(const INDEX x1, const INDEX x2) const { return msg12_[x1*dim2_ + x2]; }
   REAL msg13(const INDEX x1, const INDEX x3) const { return msg13_[x1*dim3_ + x3]; }
   REAL msg23(const INDEX x2, const INDEX x3) const { return msg23_[x2*dim3_ + x3]; }

   REAL& msg12(const INDEX x1, const INDEX x2) { return msg12_[x1*dim2_ + x2]; }
   REAL& msg13(const INDEX x1, const INDEX x3) { return msg13_[x1*dim3_ + x3]; }
   REAL& msg23(const INDEX x2, const INDEX x3) { return msg23_[x2*dim3_ + x3]; }

   void init_primal() {
      primal_[0] = dim1_;
      primal_[1] = dim2_;
      primal_[2] = dim3_;
   }
   template<class ARCHIVE> void serialize_primal(ARCHIVE& ar) { ar(primal_); }
   template<class ARCHIVE> void serialize_dual(ARCHIVE& ar) 
   { 
      //ar( cereal::binary_data( msg12_, sizeof(REAL)*(dim1()*dim2()) ) );
      //ar( cereal::binary_data( msg13_, sizeof(REAL)*(dim1()*dim3()) ) );
      //ar( cereal::binary_data( msg23_, sizeof(REAL)*(dim2()*dim3()) ) );
      ar( binary_data<REAL>( msg12_, sizeof(REAL)*(dim1()*dim2()) ) );
      ar( binary_data<REAL>( msg13_, sizeof(REAL)*(dim1()*dim3()) ) );
      ar( binary_data<REAL>( msg23_, sizeof(REAL)*(dim2()*dim3()) ) );
   }


#ifdef WITH_SAT
   template<typename SAT_SOLVER>
   void construct_sat_clauses(SAT_SOLVER& s) const
   {
      //auto vars = create_sat_variables(s, dim1()*dim2() + dim1()*dim3() + dim2()*dim3() + dim1()*dim2()*dim3());
      auto vars_12 = create_sat_variables(s, dim1()*dim2());
      auto vars_13 = create_sat_variables(s, dim1()*dim3());
      auto vars_23 = create_sat_variables(s, dim2()*dim3());
      auto vars_123 = create_sat_variables(s, dim1()*dim2()*dim3());
      //auto triplet_var_begin = vars[dim1()*dim2() + dim1()*dim3() + dim2()*dim3()];

      std::vector<sat_var> tmp_vars;
      tmp_vars.reserve(std::max({dim1()*dim2(), dim1()*dim3(), dim2()*dim3()}));

      for(INDEX x1=0; x1<dim1(); ++x1) {
         for(INDEX x2=0; x2<dim2(); ++x2) {
            for(INDEX x3=0; x3<dim3(); ++x3) {
               tmp_vars.push_back(vars_123[x1*dim2()*dim3() + x2*dim3() + x3]);
            }
            auto c = add_at_most_one_constraint_sat(s, tmp_vars.begin(), tmp_vars.end());
            make_sat_var_equal(s, to_literal(c), to_literal(vars_12[x1*dim2() + x2]));
            tmp_vars.clear();
         }
      }

      for(INDEX x1=0; x1<dim1(); ++x1) {
         for(INDEX x3=0; x3<dim3(); ++x3) {
            for(INDEX x2=0; x2<dim2(); ++x2) {
               tmp_vars.push_back(vars_123[x1*dim2()*dim3() + x2*dim3() + x3]);
            }
            auto c = add_at_most_one_constraint_sat(s, tmp_vars.begin(), tmp_vars.end());
            make_sat_var_equal(s, to_literal(c), to_literal(vars_13[x1*dim3() + x3]));
            tmp_vars.clear();
         }
      }

      for(INDEX x2=0; x2<dim2(); ++x2) {
         for(INDEX x3=0; x3<dim3(); ++x3) {
            for(INDEX x1=0; x1<dim1(); ++x1) {
               tmp_vars.push_back(vars_123[x1*dim2()*dim3() + x2*dim3() + x3]);
            }
            auto c = add_at_most_one_constraint_sat(s, tmp_vars.begin(), tmp_vars.end());
            make_sat_var_equal(s, to_literal(c), to_literal(vars_23[x2*dim3() + x3]));
            tmp_vars.clear();
         }
      }

      // summation constraints over triplet variables are superfluous: they are enforced via messages
   }

   template<typename VEC>
   void reduce_sat(VEC& assumptions, const REAL th, sat_var begin) const
   {
      begin += dim1()*dim2() + dim1()*dim3() + dim2()*dim3();
      const REAL lb = LowerBound();
      for(INDEX x1=0; x1<this->dim1(); ++x1) {
         for(INDEX x2=0; x2<this->dim2(); ++x2) {
            for(INDEX x3=0; x3<this->dim3(); ++x3) {
               if((*this)(x1,x2,x3) > lb + th) { 
                  assumptions.push_back(-to_literal(begin + x1*dim2()*dim3() + x2*dim3() + x3));
               }
            }
         }
      } 
   }

   template<typename SAT_SOLVER>
   void convert_primal(SAT_SOLVER& s, sat_var first)
   {
      for(INDEX x1=0; x1<this->dim1(); ++x1) {
         for(INDEX x2=0; x2<this->dim2(); ++x2) {
            if(lglderef(s,to_literal(first + x1*dim2() + x2)) == 1) {
               primal_[0] = x1;
               primal_[1] = x2;
            }
         } 
      }

      const INDEX pairwise_var_begin_13 = first + dim1()*dim2(); 
      for(INDEX x1=0; x1<dim1(); ++x1) {
         for(INDEX x3=0; x3<dim3(); ++x3) {
            if(lglderef(s,to_literal(pairwise_var_begin_13 + x1*dim3() + x3)) == 1) {
               assert(primal_[0] == x1);
               primal_[2] = x3;
            }
         }
      }
   }
#endif

   const std::array<INDEX,3>& primal() const { return primal_; }
   std::array<INDEX,3>& primal() { return primal_; }

protected:
   std::array<INDEX,3> primal_;
   const INDEX dim1_, dim2_, dim3_; // do zrobienia: possibly use 32 bit, for primal as well
   REAL *msg12_, *msg13_, *msg23_;
};


/*
class TighteningTernarySimplexFactor {
public:
   TighteningTernarySimplexFactor(const INDEX dim1, const INDEX dim2, const INDEX dim3) : dim1_(dim1), dim2_(dim2), dim3_(dim3) 
   {
      const INDEX size = dim1_*dim2_ + dim1_*dim3_ + dim2_*dim3_;
      
      msg12_ = global_real_block_allocator.allocate(size + std::ceil((1.0*sizeof(INDEX))/(1.0*sizeof(REAL))*size));
      assert(msg12_ != nullptr);
      msg13_ = msg12_ + dim1_*dim2_;
      msg23_ = msg13_ + dim1_*dim3_;
      std::fill(msg12_, msg12_ + size, 0.0);

      msg12_sorted_ = reinterpret_cast<INDEX*>(msg12_ + size);
      msg13_sorted_ = msg12_sorted_ + dim1_*dim2_;
      msg23_sorted_ = msg13_sorted_ + dim1_*dim3_;
      std::iota(msg12_sorted_, msg12_sorted_ + dim1_*dim2_, 0);
      std::iota(msg13_sorted_, msg13_sorted_ + dim1_*dim3_, 0);
      std::iota(msg23_sorted_, msg23_sorted_ + dim1_*dim2_, 0);
   }

   // slightly extended algorithms as used by Caetano et al for faster message passing
   void sort_indices() const { // do zrobienia: use insertion sort as it is faster on nearly sorted data
      std::sort(msg12_sorted_, msg12_sorted_ + dim1_*dim2_, [&](const INDEX a, const INDEX b) { return msg12_[msg12_sorted_[a]] < msg12_[msg12_sorted_[b]]; });
      std::sort(msg13_sorted_, msg13_sorted_ + dim1_*dim3_, [&](const INDEX a, const INDEX b) { return msg13_[msg13_sorted_[a]] < msg13_[msg13_sorted_[b]]; });
      std::sort(msg23_sorted_, msg23_sorted_ + dim1_*dim2_, [&](const INDEX a, const INDEX b) { return msg23_[msg23_sorted_[a]] < msg23_[msg23_sorted_[b]]; });
   }

   REAL FindBest3(const INDEX x1, const INDEX x2) const {
      REAL min = std::numeric_limits<REAL>::infinity();
      INDEX end13 = msg13_sorted_inverse(x1,msg23_sorted(x2,0));
      INDEX end23 = msgmsg23_sorted_inverse(msg13_sorted(x1,0));
      for(INDEX l3=0; l3<dim3_; ++l3) {
         if(l3 >= std::min(end13, end23)) break;
         {
            const INDEX x3 = msg13_sorted(x1,l3);
            min_val = std::min(min_val, msg_13(x1,x3) + msg_23(x2,x3));
            end12 = std::min(end12,msg23_sorted_inverse(x2, x3) );
         }
         {
            const INDEX x3 = msg23_sorted(x2,l3);
            min_val = std::min(min_val, msg_13(x1,x3) + msg_23(x2,x3));
            end23 = std::min(end23,msg13_sorted_inverse(x2, x3) );
         }
      }
      return min_val;
   }

   REAL LowerBound() const {
      REAL val = std::numeric_limits<REAL>::infinity();
      sort_indices();
      for(INDEX l1=0; l1<dim1_*dim2_; ++l1) { // add early stopping by checking whether val < msg12_sorted(l1)
         const INDEX x1 = msg12_sorted_1(l1);
         const INDEX x2 = msg12_sorted_2(l1);
         // find best index x3 such that
         val = std::min(val, msg12(x1,x2), FindBest3(x1,x2));
      }

      return val;
   }

   REAL& msg12(const INDEX x1, const INDEX x2) { return msg12_[x1*dim2_ + x2]; }
   REAL& msg13(const INDEX x1, const INDEX x3) { return msg13_[x1*dim3_ + x3]; }
   REAL& msg23(const INDEX x2, const INDEX x3) { return msg12_[x2*dim3_ + x3]; }
   INDEX msg12_sorted_1(const INDEX x) const { return msg12_sorted_[x]/dim1_; }
   INDEX msg12_sorted_2(const INDEX x) const { return msg12_sorted_[x]%dim1_; }
private:
   const INDEX dim1_, dim2_, dim3_;
   REAL *msg12_, *msg13_, *msg23_;
   mutable INDEX *msg12_sorted_, *msg13_sorted_, *msg23_sorted_;

public:

};
*/

class pairwise_potts_factor : public vector<REAL> {
   template<typename MATRIX>
   static bool is_potts(const MATRIX& c)
   {
      if(c.dim1() != c.dim2()) { return false; }
      if(c.dim1() <= 1) { return false; }
      for(INDEX x1=0; x1<c.dim1(); ++x1) {
         for(INDEX x2=0; x2<c.dim2(); ++x2) {
            if(x1 == x2) {
              if(c(x1,x2) != 0.0) return false;
            } else {
               if(c(x1,x2) != c(0,1)) return false;
            }
         }
      }
      return true;
   }
public:
   template<typename VEC>
   pairwise_potts_factor(const INDEX dim1, const INDEX dim2, const VEC& cost)
   : pairwise_potts_factor(dim1, cost(0,1))
   {
      assert(dim1 == dim2);
      assert(is_potts(cost));
      assert(cost(0,0) == 0.0);
   }

   pairwise_potts_factor(const INDEX dim1, const INDEX dim2)
      : pairwise_potts_factor(dim1, REAL(0.0))
   {
      assert(dim1 == dim2);
   }

   pairwise_potts_factor(const INDEX dim1, const INDEX dim2, const REAL diff_cost)
      : pairwise_potts_factor(dim1, diff_cost)
   {
      assert(dim1 == dim2); 
   } 

   pairwise_potts_factor(const INDEX dim, const REAL diff_cost)
      : 
         vector<REAL>(2*dim,0.0), // holds reparametrisation
         diff_cost_(diff_cost) 
   {}

   REAL operator()(const INDEX x1, const INDEX x2) const
   {
      const REAL msg_val = msg1(x1) + msg2(x2);
      const REAL potts_val = x1 != x2 ? diff_cost() : 0.0;
      return msg_val + potts_val;
   }

   // min cost for same label and different label
   std::array<REAL,2> min_values() const
   {
      const auto smallest2 = two_smallest_elements<REAL>(msg2_begin(), msg2_end());

      REAL min_same_label = std::numeric_limits<REAL>::infinity();
      REAL min_diff_label = std::numeric_limits<REAL>::infinity();
      for(INDEX i=0; i<dim(); ++i) {
         const REAL same_label = (*this)[i] + (*this)[i+dim()];
         min_same_label = std::min(min_same_label, same_label);
         const REAL diff_label = (*this)[i] + diff_cost() + ((*this)[i+dim()] == smallest2[0] ? smallest2[1] : smallest2[0]);
         min_diff_label = std::min(min_diff_label, diff_label);
      }
      return {min_same_label, min_diff_label}; 
   }

   REAL LowerBound() const {
      const auto v = min_values();
      return std::min(v[0], v[1]);
   }

   REAL EvaluatePrimal() const
   {
      if(primal_[0] < dim() && primal_[1] < dim()) {
         const REAL pairwise_cost = primal_[0] == primal_[1] ? 0.0 : diff_cost_;
         const REAL msg_cost = (*this)[primal_[0]] + (*this)[primal_[1] + dim()];
         return pairwise_cost + msg_cost; 
      } else {
         return std::numeric_limits<REAL>::infinity();
      }
   }

   template<typename VECTOR>
   void min_marginal_1(VECTOR& m) const
   {
      assert(m.size() == dim());

      const auto smallest2 = two_smallest_elements<REAL>(msg2_begin(), msg2_end());

      for(INDEX i=0; i<dim(); ++i) {
         const REAL same_label = (*this)[i+dim()];
         const REAL diff_label = diff_cost() + ((*this)[i+dim()] == smallest2[0] ? smallest2[1] : smallest2[0]);
         m[i] = (*this)[i] + std::min(same_label, diff_label); 
      } 
   }

   template<typename VECTOR>
   void min_marginal_2(VECTOR& m) const
   {
      assert(m.size() == dim());

      const auto smallest2 = two_smallest_elements<REAL>(msg1_begin(), msg1_end());

      for(INDEX i=0; i<dim(); ++i) {
         const REAL same_label = (*this)[i];
         const REAL diff_label = diff_cost() + ((*this)[i] == smallest2[0] ? smallest2[1] : smallest2[0]);
         m[i] = (*this)[i+dim()] + std::min(same_label, diff_label); 
      } 
   }

   REAL min_marginal_cut() const
   {
      const auto v = min_values();
      return v[1] - v[0];
      //return v[0] - v[1];
   }

   INDEX dim() const { return this->size()/2; }
   INDEX dim1() const { return dim(); }
   INDEX dim2() const { return dim(); }

   REAL* msg1_begin() const { return this->begin(); }
   REAL* msg1_end() const { return this->begin() + dim(); }
   REAL msg1(const INDEX x1) const { assert(x1 < dim()); return (*this)[x1]; }
   REAL& msg1(const INDEX x1) { assert(x1 < dim()); return (*this)[x1]; }

   REAL* msg2_begin() const { return this->begin() + dim(); }
   REAL* msg2_end() const { return this->end(); }
   REAL msg2(const INDEX x2) const { return (*this)[dim() + x2]; }
   REAL& msg2(const INDEX x2) { return (*this)[dim() + x2]; }

   REAL diff_cost() const { return diff_cost_; }
   REAL& diff_cost() { return diff_cost_; }

   void init_primal() { primal_[0] = std::numeric_limits<INDEX>::max(); primal_[1] = std::numeric_limits<INDEX>::max(); }
   auto& primal() { return primal_; }
   const auto& primal() const { return primal_; }

   template<typename ARCHIVE> void serialize_dual(ARCHIVE& ar) { ar( *static_cast<vector<REAL>*>(this), diff_cost_ ); }
   template<typename ARCHIVE> void serialize_primal(ARCHIVE& ar) { ar( primal_ ); }

protected:
   REAL diff_cost_;
   std::array<INDEX,2> primal_;

};
} // end namespace LP_MP

#endif // LP_MP_SIMPLEX_FACTOR_HXX

