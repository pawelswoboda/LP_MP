#ifndef LP_MP_SIMPLEX_FACTOR_HXX
#define LP_MP_SIMPLEX_FACTOR_HXX

#include "LP_MP.h"
#include "memory_allocator.hxx"
#include "vector.hxx"

// Investigate contingency tables (Knuth) for more general tabular structures.

namespace LP_MP {

// the polytope {x >= 0 : x_1 + ... + x_n = 1 }
// dual problem:
//      max_{z} z 
//      s.t.    z <= repamPot[i]

// do zrobienia: replace std::vector by dynamic vector that is not resizable and gets memory from factor allocator
template<typename T = REAL> // artifact: do zrobienia: remove
class SimplexFactor : public std::vector<REAL>
{
public:
   using VECTOR = std::vector<REAL>;
   SimplexFactor(const std::vector<REAL>& cost)
      : VECTOR(cost.begin(), cost.end()) // better use iterators and allocator at same time
   {}
   SimplexFactor(const INDEX n)
      : VECTOR(n,0.0)
   {}

   void MaximizePotentialAndComputePrimal(typename PrimalSolutionStorage::Element primal) const
   {
      INDEX min_element;
      REAL min_value = std::numeric_limits<REAL>::infinity();
      for(INDEX i=0; i<this->size(); ++i) {
         if(primal[i] == true) {
            return;
         } else if(primal[i] == false) {
            // do nothing
         } else {
            assert(primal[i] == unknownState);
            primal[i] = false;
            if(min_value >= (*this)[i]) {
               min_value = (*this)[i];
               min_element = i;
            }
         }
      }
      assert(min_element < this->size());
      primal[min_element] = true;
   };

   REAL LowerBound() const {
      return *std::min_element(this->begin(), this->end());
   }

   // if there is exactly one unknownIndex (rest false), make it true
   // if there is already one true index, make all other false
   void PropagatePrimal(PrimalSolutionStorage::Element primal)
   {
      INDEX noTrue = 0;
      INDEX noUnknown = 0;
      INDEX unknownIndex;
      for(INDEX i=0; i<this->size(); ++i) {
         if(primal[i] == true) { 
            ++noTrue;
         }
         if(primal[i] == unknownState) {
            ++noUnknown;
            unknownIndex = i;
         }
      }
      if(noTrue == 0 && noUnknown == 1) {
         primal[unknownIndex] = true;
      } else if(noTrue == 1 && noUnknown > 0) {
         for(INDEX i=0; i<this->size(); ++i) {
            if(primal[i] == unknownState) {
               primal[i] = false;
            }
         }
      }
   }

   REAL EvaluatePrimal(const PrimalSolutionStorage::Element primal) const
   {
      REAL cost;
      INDEX primalSum = 0;
      for(INDEX i=0; i<this->size(); ++i) {
         if(primal[i] == true) { 
            cost = (*this)[i];
         }
         primalSum += primal[i];
      }
      if(primalSum == 1) {
         return cost;
      } else {
         //std::cout << "primal not inferred correctly: " << primalSum << "\n";
         return std::numeric_limits<REAL>::infinity();
      }
   }

  //INDEX GetNumberOfAuxVariables() const { return 0; }

  /*
  void ReduceLp(LpInterfaceAdapter* lp) const {
    REAL lb = LowerBound();
    REAL epsi = lp->GetEpsilon();
    for(INDEX i=0;i<this->size();i++){
      if((*this)[i] >= lb + epsi){
        lp->SetVariableBound(lp->GetVariable(i),0.0,0.0,false);
      }
    }
  }
  
  void CreateConstraints(LpInterfaceAdapter* lp) const { 
    LinExpr lhs = lp->CreateLinExpr();
    for(INDEX i=0;i<lp->GetFactorSize();i++){
      lhs += lp->GetVariable(i);
    }
    LinExpr rhs = lp->CreateLinExpr();
    rhs += 1;
    lp->addLinearEquality(lhs,rhs);
  }
  */
  
};

   
class UnarySimplexFactor : public vector {
public:
   UnarySimplexFactor(const std::vector<REAL>& cost)
      : vector(cost.begin(), cost.end())
   {}
   UnarySimplexFactor(const INDEX n)
      : vector(n, 0.0)
   {}

   void MaximizePotentialAndComputePrimal(typename PrimalSolutionStorage::Element primal) const
   {
      INDEX min_element;
      REAL min_value = std::numeric_limits<REAL>::infinity();
      for(INDEX i=0; i<this->size(); ++i) {
         if(primal[i] == true) {
            return;
         } else if(primal[i] == false) {
            // do nothing
         } else {
            assert(primal[i] == unknownState);
            primal[i] = false;
            if(min_value >= (*this)[i]) {
               min_value = (*this)[i];
               min_element = i;
            }
         }
      }
      assert(min_element < this->size());
      primal[min_element] = true;
   };

   REAL LowerBound() const {
      return *std::min_element(this->begin(), this->end());
   }

   // if there is exactly one unknownIndex (rest false), make it true
   // if there is already one true index, make all other false
   void PropagatePrimal(PrimalSolutionStorage::Element primal)
   {
      INDEX noTrue = 0;
      INDEX noUnknown = 0;
      INDEX unknownIndex;
      for(INDEX i=0; i<this->size(); ++i) {
         if(primal[i] == true) { 
            ++noTrue;
         }
         if(primal[i] == unknownState) {
            ++noUnknown;
            unknownIndex = i;
         }
      }
      if(noTrue == 0 && noUnknown == 1) {
         primal[unknownIndex] = true;
      } else if(noTrue == 1 && noUnknown > 0) {
         for(INDEX i=0; i<this->size(); ++i) {
            if(primal[i] == unknownState) {
               primal[i] = false;
            }
         }
      }
   }

   REAL EvaluatePrimal(const PrimalSolutionStorage::Element primal) const
   {
      REAL cost;
      INDEX primalSum = 0;
      for(INDEX i=0; i<this->size(); ++i) {
         if(primal[i] == true) { 
            cost = (*this)[i];
         }
         primalSum += primal[i];
      }
      if(primalSum == 1) {
         return cost;
      } else {
         //std::cout << "primal not inferred correctly: " << primalSum << "\n";
         return std::numeric_limits<REAL>::infinity();
      }
   }
};

// do zrobienia: if pairwise was supplied to us (e.g. external factor, then reflect this in constructor and only allocate space for messages.
// When tightening, we can simply replace pairwise pointer to external factor with an explicit copy. Reallocate left_msg_ and right_msg_ to make memory contiguous? Not sure, depends whether we use block_allocator, which will not acually release the memory
// when factor is copied, then pairwise_ must only be copied if it is actually modified. This depends on whether we execute SMRP or MPLP style message passing. Templatize for this possibility
class PairwiseSimplexFactor {
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
   // below is not nice: two values, only differ by const!
   REAL operator()(const INDEX x1, const INDEX x2) const {
      assert(x1 < dim1_ && x2 < dim2_);
      return pairwise_[x1*dim2_ + x2] + left_msg_[x1] + right_msg_[x2];
   }
   REAL& operator()(const INDEX x1, const INDEX x2) {
      assert(x1 < dim1_ && x2 < dim2_);
      return pairwise_[x1*dim2_ + x2];
   }
   REAL LowerBound() const {
      REAL lb = std::numeric_limits<REAL>::infinity();
      for(INDEX x1=0; x1<dim1_; ++x1) {
         for(INDEX x2=0; x2<dim2_; ++x2) {
            //std::cout << (*this)(x1,x2) << ", ";
            lb = std::min(lb, (*this)(x1,x2));
         }
      }
      //std::cout << "\n";
      return lb;
   }

   REAL EvaluatePrimal(const PrimalSolutionStorage::Element primal) const
   {
      REAL cost;
      INDEX primalSum = 0;
      for(INDEX i=0; i<dim1_*dim2_; ++i) {
         if(primal[i] == true) { 
            cost = (*this)[i];
         }
         primalSum += primal[i];
      }
      if(primalSum == 1) {
         return cost;
      } else {
         return std::numeric_limits<REAL>::infinity();
      }
   }

   // do zrobienia: propagate primal is needed

   
   INDEX dim1() const { return dim1_; }
   INDEX dim2() const { return dim2_; }
     
   REAL& pairwise(const INDEX x1, const INDEX x2) {
      return pairwise_[x1*dim2_ + x2];
   }
   REAL& left(const INDEX x1) {
      return left_msg_[x1];
   }
   REAL& right(const INDEX x2) {
      return right_msg_[x2];
   }

   INDEX size() const {
      return dim1_*dim2_;
   }

private:
   // those three pointers should lie contiguously in memory
   REAL* pairwise_;
   REAL* left_msg_;
   REAL* right_msg_;
   const INDEX dim1_, dim2_;
};

// factor assumes that triplet potentials is empty and holds only messages to pairwise factors, i.e. is latently factorizable
class SimpleTighteningTernarySimplexFactor {
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
      REAL val = std::numeric_limits<REAL>::infinity();
      for(INDEX x1=0; x1<dim1_; ++x1) {
         for(INDEX x2=0; x2<dim2_; ++x2) {
            for(INDEX x3=0; x3<dim3_; ++x3) {
               val = std::min(val, (*this)(x1,x2,x3));
            }
         }
      }
      return val;
   }

   REAL EvaluatePrimal(const PrimalSolutionStorage::Element primal) const
   {
      REAL cost;
      INDEX primalSum = 0;

      for(INDEX i=0; i<dim1_*dim2_; ++i) {
         if(primal[i] == true) { 
            cost = (*this)[i];
         }
         primalSum += primal[i];
      }
      if(primalSum == 1) {
         return cost;
      } else {
         return std::numeric_limits<REAL>::infinity();
      }
   }

   REAL operator()(const INDEX x1, const INDEX x2, const INDEX x3) const {
      return msg12(x1,x2) + msg13(x1,x3) + msg23(x2,x3);
   }

   REAL operator[](const INDEX x) const {
      const INDEX x1 = x / (dim2_*dim3_);
      const INDEX x2 = ( x % (dim2_*dim3_) ) / dim3_;
      const INDEX x3 = x % dim3_;
      return msg12(x1,x2) + msg13(x1,x3) + msg23(x2,x3);
   }

   REAL min_marginal12(const INDEX x1, const INDEX x2) const {
      REAL marg = (*this)(x1,x2,0);
      for(INDEX x3=1; x3<dim3_; ++x3) {
         marg = std::min(marg, (*this)(x1,x2,x3));
      }
      return marg;
   }
   template<typename MSG>
   void min_marginal12(MSG& msg) const {
      for(INDEX x1=0; x1<dim1_; ++x1) {
         for(INDEX x2=0; x2<dim2_; ++x2) {
            msg(x1,x2) = std::numeric_limits<REAL>::infinity();
         }
      }
      for(INDEX x1=0; x1<dim1_; ++x1) {
         for(INDEX x2=0; x2<dim2_; ++x2) {
            msg(x1,x2) = min_marginal12(x1,x2);
         }
      }
   }

   template<typename MSG>
   void min_marginal13(MSG& msg) const {
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
   REAL& msg23(const INDEX x2, const INDEX x3) { return msg12_[x2*dim3_ + x3]; }
protected:
   const INDEX dim1_, dim2_, dim3_;
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

} // end namespace LP_MP

#endif // LP_MP_SIMPLEX_FACTOR_HXX

