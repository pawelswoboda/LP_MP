#ifndef LP_MP_DISCRETE_TOMOGRAPHY_SEQUENTIAL_HXX
#define LP_MP_DISCRETE_TOMOGRAPHY_SEQUENTIAL_HXX

#include "config.hxx"

namespace LP_MP {

// class holding state of variable and sum of chain up to that variable
class dt_sum_state_factor {
public:
   constexpr static INDEX no_primal_decision = std::numeric_limits<INDEX>::max();

   dt_sum_state_factor(const INDEX no_labels, const INDEX sum_size) 
      : pot_(no_labels, sum_size, 0.0)
   {
   assert(pot_.dim1() > 0);
   }

   template<typename ARRAY>
   void summation_cost(const ARRAY& cost) {
      assert(cost.size() <= sum_size());
      for(INDEX x=0; x<pot_.dim1(); ++x) {
         for(INDEX sum=0; sum<pot_.dim2(); ++sum) {
            if(x+sum < cost.size()) {
               pot_(x,sum) = cost[x+sum];
            } else {
               pot_(x,sum) = std::numeric_limits<REAL>::infinity();
            } 
         }
      }
   }

   REAL LowerBound() const {
      return *std::min_element(pot_.begin(), pot_.end());
      //REAL min_value = std::numeric_limits<REAL>::infinity();
      //for(auto it=pot_.begin(); it!=pot_.end(); ++it) {
      //   min_value = std::min(min_value, *it);
      //}
      //return min_value;
   }

   REAL EvaluatePrimal() const {
      return std::numeric_limits<REAL>::infinity();
   }

   void MaximizePotentialAndComputePrimal() {
      REAL min_value = std::numeric_limits<REAL>::infinity();
      assert(sum_ == no_primal_decision && state_ == no_primal_decision);
      for(INDEX state=0; state<no_labels(); ++state) {
         for(INDEX sum=0; sum<sum_size(); ++sum) {
            if((*this)(state, sum) <= min_value) {
               state_ = state;
               sum_ = sum;
               min_value = (*this)(state,sum);
            }
         }
      } 
   }

   REAL operator()(const INDEX state, const INDEX sum) const { return pot_(state,sum); }
   REAL& operator()(const INDEX state, const INDEX sum) { return pot_(state,sum); }
   const matrix<REAL>& pot() const { return pot_; }

   INDEX no_labels() const { return pot_.dim1(); }
   INDEX sum_size() const { return pot_.dim2(); }
   INDEX size() const { return pot_.size(); }

   void init_primal() { state_ = no_primal_decision; sum_ = no_primal_decision; }
   template<class ARCHIVE> void serialize_primal(ARCHIVE& ar) { ar( state_, sum_ ); }
   template<class ARCHIVE> void serialize_dual(ARCHIVE& ar) { ar(pot_); }

   INDEX state_;
   INDEX sum_;
private:
   matrix<REAL> pot_; // first dimension is state, second one is sum
};

class dt_sum_state_pairwise_factor {
public:
   constexpr static INDEX no_primal_decision = std::numeric_limits<INDEX>::max();

   dt_sum_state_pairwise_factor(const INDEX no_labels, const INDEX prev_sum_size, const INDEX next_sum_size) 
      : prev_(no_labels, prev_sum_size, 0.0),
      next_(no_labels, next_sum_size, 0.0),
      reg_(no_labels, no_labels, 0.0)
   {
      assert(prev_sum_size <= next_sum_size);
   }

   template<typename LAMBDA>
   void for_each_label_sum(const LAMBDA&& f) const {
      assert(prev_sum_size() <= next_sum_size());
      assert(prev_sum_size() + no_labels()-1 >= next_sum_size());
      for(INDEX x1=0; x1<no_labels(); ++x1) {
         for(INDEX x2=0; x2<no_labels(); ++x2) {
            const INDEX max_sum = std::min(prev_sum_size(), next_sum_size()-x1);
            for(INDEX sum=0; sum<max_sum; ++sum) {
               f(x1,x2,sum);
            }
         }
      }
   }

   REAL eval(const INDEX x1, const INDEX x2, const INDEX sum) const {
      assert(x1 < no_labels() && x2 < no_labels() && sum < prev_sum_size() && sum+x1 < next_sum_size());
      return prev_(x1,sum) + reg_(x1,x2) + next_(x2, sum+x1);
   }
   REAL LowerBound() const {
      REAL min_value = std::numeric_limits<REAL>::infinity();
      for_each_label_sum([&](INDEX x1, INDEX x2, INDEX sum) { min_value = std::min(min_value, eval(x1,x2,sum)); });
      return min_value;
   }

   REAL EvaluatePrimal() const {
      return std::numeric_limits<REAL>::infinity();
   }

   void marginalize_pairwise(matrix<REAL>& msg) const {
      std::fill(msg.begin(), msg.end(), std::numeric_limits<REAL>::infinity());
      for_each_label_sum([&](INDEX x1, INDEX x2, INDEX sum) { msg(x1,x2) = std::min(msg(x1,x2), eval(x1,x2,sum)); });
   }

   REAL& prev(const INDEX state, const INDEX sum) { return prev_(state,sum); }
   REAL& next(const INDEX state, const INDEX sum) { return next_(state,sum); }
   REAL& reg(const INDEX x1, const INDEX x2) { return reg_(x1,x2); }
   REAL prev(const INDEX state, const INDEX sum) const { return prev_(state,sum); }
   REAL next(const INDEX state, const INDEX sum) const { return next_(state,sum); }
   REAL reg(const INDEX x1, const INDEX x2) const { return reg_(x1,x2); }

   INDEX no_labels() const { return reg_.dim2(); }
   INDEX prev_sum_size() const { return prev_.dim2(); }
   INDEX next_sum_size() const { return next_.dim2(); }
   INDEX size() const { return no_labels()*no_labels()*prev_sum_size(); }

   void init_primal() { state_[0] = no_primal_decision; state_[1] = no_primal_decision; sum_[0] = no_primal_decision; sum_[1] = no_primal_decision; }
   template<class ARCHIVE> void serialize_primal(ARCHIVE& ar) { ar( state_, sum_ ); }
   template<class ARCHIVE> void serialize_dual(ARCHIVE& ar) { ar( prev_, next_, reg_ ); }

   std::array<INDEX,2> state_;
   std::array<INDEX,2> sum_; // both sums are not strictly needed, they help however in labeling
private:
   matrix<REAL> prev_; // first dimension is state, second one is sum
   matrix<REAL> next_; // first dimension is state, second one is sum
   matrix<REAL> reg_; // regularizer
};

// message between unaries and sum, not necessary for the relaxation. Also does not seem to improve convergence
/*
class dt_unary_sum_message {
public:
   dt_unary_sum_message() {}

   template<typename RIGHT_FACTOR, typename G2>
   void ReceiveMessageFromRight(const RIGHT_FACTOR& f_right, G2& msg){
      MakeRightFactorUniform(f_right, msg, 1.0);
   }

   //template<typename RIGHT_FACTOR, typename G2>
   //void SendMessageToLeft(const RIGHT_FACTOR& f_right, G2& msg, const REAL omega){
   //   MakeRightFactorUniform(f_right, msg, omega);
   //}

   template<typename LEFT_FACTOR, typename G3>
   void SendMessageToRight(const LEFT_FACTOR& f_left, G3& msg, const REAL omega){
      MakeLeftFactorUniform(f_left, msg, omega);
   }

   //template<typename LEFT_FACTOR, typename G3>
   //void ReceiveMessageFromLeft(const LEFT_FACTOR& f_left, G3& msg){
   //   MakeLeftFactorUniform(f_left, msg, 1.0);
   //}

   template<typename LEFT_FACTOR, typename MSG>
   void MakeLeftFactorUniform(const LEFT_FACTOR& f_left, MSG& msg, const REAL omega){
      auto& pot = f_left;
      msg -= omega*pot;
   }

   template<typename RIGHT_FACTOR, typename MSG>
   void MakeRightFactorUniform(const RIGHT_FACTOR& f_right, MSG& msg, const REAL omega){
      vector marg(f_right.no_labels(), std::numeric_limits<REAL>::infinity());
      for(INDEX x=0; x<f_right.no_labels(); ++x) {
         for(INDEX sum=0; sum<f_right.sum_size(); ++sum) {
            marg[x] = std::min(marg[x], f_right(x,sum));
         }
      } 
      msg -= omega*marg;
   }

   template<typename LEFT_FACTOR, typename MSG>
   void RepamLeft(LEFT_FACTOR& f_left, const MSG& msg)
   {
      for(INDEX x=0; x<f_left.size(); ++x) {
         f_left[x] += normalize( msg[x] );
      }
   }

   template<typename RIGHT_FACTOR, typename MSG>
   void RepamRight(RIGHT_FACTOR& f_right, const MSG& msg)
   {
      for(INDEX x=0; x<f_right.no_labels(); ++x) {
         for(INDEX sum=0; sum<f_right.sum_size(); ++sum) {
            f_right(x,sum) += normalize( msg[x] );
         }
      }
   }
};
*/

// left factor is dt_sum_state_factor, right one is dt_sum_state_pairwise_factor
// DIRECTION signifies: left = sum_state_factor is previous
//                      right = sum_state_factor is next
// possibly make functions static
template<Chirality DIRECTION>
class dt_sum_pairwise_message {
public:
   dt_sum_pairwise_message(bool transpose) : transpose_(transpose) {
      //assert(transpose_ == false);
   }

   ~dt_sum_pairwise_message() {
      static_assert(DIRECTION == Chirality::left || DIRECTION == Chirality::right, "");
   }

   template<typename RIGHT_FACTOR, typename G2>
   void ReceiveMessageFromRight(const RIGHT_FACTOR& f_right, G2& msg){
      MakeRightFactorUniform(f_right, msg, 1.0);
   }

   //template<typename RIGHT_FACTOR, typename G2>
   //void SendMessageToLeft(const RIGHT_FACTOR& f_right, G2& msg, const REAL omega){
   //   MakeRightFactorUniform(f_right, msg, omega);
   //}

   template<typename LEFT_FACTOR, typename G3>
   void SendMessageToRight(const LEFT_FACTOR& f_left, G3& msg, const REAL omega){
      //std::cout << "sum to sum pairwise weight = " << omega << "\n";
      MakeLeftFactorUniform(f_left, msg, omega);
   }

   //template<typename LEFT_FACTOR, typename G3>
   //void ReceiveMessageFromLeft(const LEFT_FACTOR& f_left, G3& msg){
   //   MakeLeftFactorUniform(f_left, msg, 1.0);
   //}

   template<typename LEFT_FACTOR, typename MSG>
   void MakeLeftFactorUniform(const LEFT_FACTOR& f_left, MSG& msg, const REAL omega){
      auto& pot = f_left.pot();
      msg -= omega*pot;
   }

   template<typename RIGHT_FACTOR, typename MSG>
   void MakeRightFactorUniform(const RIGHT_FACTOR& f_right, MSG& msg, const REAL omega){
      matrix<REAL> msgs(f_right.no_labels(), DIRECTION == Chirality::left ? f_right.prev_sum_size() : f_right.next_sum_size(), std::numeric_limits<REAL>::infinity());
      for(INDEX x1=0; x1<f_right.no_labels(); ++x1) {
         for(INDEX x2=0; x2<f_right.no_labels(); ++x2) {
            if(DIRECTION == Chirality::left) {
               for(INDEX sum=0; sum<f_right.prev_sum_size(); ++sum) {
                  if(sum+x1 < f_right.next_sum_size()) {
                     msgs(x1,sum) = std::min(msgs(x1,sum), f_right.eval(x1,x2,sum));
                  }
               }
            } else {
               for(INDEX sum=0; sum<f_right.prev_sum_size(); ++sum) {
                  if(sum+x1 < f_right.next_sum_size()) {
                     msgs(x2,sum+x1) = std::min(msgs(x2,sum+x1), f_right.eval(x1,x2,sum));
                  }
               } 
            }
         }
      } 
      msg -= omega*msgs;
   }

   template<typename RIGHT_FACTOR, typename MSG>
   void ReceiveRestrictedMessageFromRight(const RIGHT_FACTOR& r, MSG& msg)
   {
      matrix<REAL> msgs(r.no_labels(), DIRECTION == Chirality::left ? r.prev_sum_size() : r.next_sum_size(), std::numeric_limits<REAL>::infinity());

      if(DIRECTION == Chirality::left) {
         assert(r.state_[1] < r.no_labels() && r.sum_[1] < r.next_sum_size());
         for(INDEX state=0; state<std::min(r.no_labels(), r.sum_[1]+1); ++state) {
            msgs(state, r.sum_[1] - state) = r.eval(state, r.state_[1], r.sum_[1] - state);
         }
      } else {
         assert(false);
      } 
      msg -= msgs;
   }

   template<typename LEFT_FACTOR, typename MSG>
   void RepamLeft(LEFT_FACTOR& f_left, const MSG& msg)
   {
      for(INDEX x=0; x<f_left.no_labels(); ++x) {
         for(INDEX sum=0; sum<f_left.sum_size(); ++sum) {
            f_left(x,sum) += normalize( msg(x,sum) );
         }
      }
   }

   template<typename RIGHT_FACTOR, typename MSG>
   void RepamRight(RIGHT_FACTOR& f_right, const MSG& msg)
   {
      for(INDEX x=0; x<f_right.no_labels(); ++x) {
         if(DIRECTION == Chirality::left) {
            for(INDEX sum=0; sum<f_right.prev_sum_size(); ++sum) {
               f_right.prev(x,sum) += normalize( msg(x,sum) );
            }
         } else {
            for(INDEX sum=0; sum<f_right.next_sum_size(); ++sum) {
               f_right.next(x,sum) += normalize( msg(x,sum) );
            }
         }
      }
   }

   template<typename LEFT_FACTOR, typename RIGHT_FACTOR>
   void ComputeRightFromLeftPrimal(const LEFT_FACTOR& l, RIGHT_FACTOR& r)
   {
      if(DIRECTION == Chirality::left) {
         r.state_[0] = l.state_;
         r.sum_[0] = l.sum_;
      } else {
         r.state_[1] = l.state_;
         r.sum_[1] = l.sum_;
      } 
   }

   /*
   template<typename LEFT_FACTOR, typename RIGHT_FACTOR>
   void ComputeLeftFromRightPrimal(LEFT_FACTOR& l, const RIGHT_FACTOR& r) {
      if(DIRECTION == Chirality::left) {
         assert(r.state_[0] < r.no_labels() && r.sum_[0] < r.prev_sum_size());
         l.state_ = r.state_[0];
         l.sum_ = r.sum_[0];
      } else {
         assert(r.state_[1] < r.no_labels() && r.sum_[1] < r.prev_sum_size());
         l.state_ = r.state_[1]; 
         l.sum_ = r.sum_[1];
      } 
   }
   */
private:
   bool transpose_; // not needed
};

// connect sequential dt_sum_factor (left) with recursive DiscreteTomograhpyCountingFactor (right)
template<Chirality DIRECTION>
class dt_sum_counting_message {
public:
   ~dt_sum_counting_message() {
      static_assert(DIRECTION == Chirality::left || DIRECTION == Chirality::right, "");
   }
   template<typename RIGHT_FACTOR, typename G2>
   void ReceiveMessageFromRight(const RIGHT_FACTOR& f_right, G2& msg){
      MakeRightFactorUniform(f_right, msg, 1.0);
   }

   //template<typename RIGHT_FACTOR, typename G2>
   //void SendMessageToLeft(const RIGHT_FACTOR& f_right, G2& msg, const REAL omega){
   //   MakeRightFactorUniform(f_right, msg, omega);
   //}

   template<typename LEFT_FACTOR, typename G3>
   void SendMessageToRight(const LEFT_FACTOR& f_left, G3& msg, const REAL omega){
      MakeLeftFactorUniform(f_left, msg, omega);
   }

   //template<typename LEFT_FACTOR, typename G3>
   //void ReceiveMessageFromLeft(const LEFT_FACTOR& f_left, G3& msg){
   //   MakeLeftFactorUniform(f_left, msg, 1.0);
   //}

   template<typename LEFT_FACTOR, typename MSG>
   void MakeLeftFactorUniform(const LEFT_FACTOR& f_left, MSG& msg, const REAL omega){
      auto& pot = f_left.pot();
      msg -= omega*pot;
   }

   template<typename RIGHT_FACTOR, typename MSG>
   void MakeRightFactorUniform(const RIGHT_FACTOR& f_right, MSG& msg, const REAL omega){
      std::array<INDEX,3> dim;

      if(DIRECTION == Chirality::left) {
         assert(1 == f_right.no_left_labels());
         dim = {1, f_right.no_center_left_labels(), f_right.left_sum_size()};
      } else {
         assert(1 == f_right.no_right_labels());
         dim = {f_right.no_center_right_labels(), 1, f_right.right_sum_size()}; 
      }
      tensor3<REAL> msgs(dim[0], dim[1], dim[2]);

      if(DIRECTION == Chirality::left) {
         f_right.MessageCalculation_Left(msgs);
         matrix_view_of_tensor<0,REAL> msgs_reduced(msgs,0);
         msg -= omega*msgs_reduced;
      } else {
         f_right.MessageCalculation_Right(msgs);
         matrix_view_of_tensor<1,REAL> msgs_reduced(msgs,0);
         msg -= omega*msgs_reduced;
      }
   }

   template<typename LEFT_FACTOR, typename MSG>
   void RepamLeft(LEFT_FACTOR& f_left, const MSG& msg)
   {
      //static_if<std::is_base_of<matrix,msg>([&](auto f) {});
      assert(msg.dim1() == f_left.no_labels());
      assert(msg.dim2() == f_left.sum_size());
      for(INDEX x=0; x<f_left.no_labels(); ++x) {
         for(INDEX sum=0; sum<f_left.sum_size(); ++sum) {
            f_left(x,sum) += normalize( msg(x,sum) );
         }
      }
   }

   template<typename RIGHT_FACTOR, typename MSG>
   void RepamRight(RIGHT_FACTOR& f_right, const MSG& msg)
   {
      if(DIRECTION == Chirality::left) {
         assert(f_right.no_left_labels() == 1);
         assert(msg.dim1() == f_right.no_center_left_labels());
         assert(msg.dim2() == f_right.left_sum_size());
         for(INDEX x=0; x<f_right.no_center_left_labels(); ++x) {
            for(INDEX sum=0; sum<f_right.left_sum_size(); ++sum) {
               f_right.left(0,x,sum) += normalize( msg(x,sum) );
            }
         }
      } else {
         assert(f_right.no_right_labels() == 1);
         assert(msg.dim1() == f_right.no_center_right_labels());
         assert(msg.dim2() == f_right.right_sum_size());
         for(INDEX x=0; x<f_right.no_center_right_labels(); ++x) {
            for(INDEX sum=0; sum<f_right.right_sum_size(); ++sum) {
               f_right.right(x,0,sum) += normalize( msg(x,sum) );
            }
         }
      }
   }
};

// connect mrf pairwise factor (left) with dt_sequential_pairwise_factor (right)
class dt_sum_pairwise_pairwise_message {
public:
   dt_sum_pairwise_pairwise_message (const bool transpose) : transpose_(transpose) {}

   template<typename RIGHT_FACTOR, typename G2>
   void ReceiveMessageFromRight(const RIGHT_FACTOR& f_right, G2& msg){
      MakeRightFactorUniform(f_right, msg, 1.0);
   }

   //template<typename RIGHT_FACTOR, typename G2>
   //void SendMessageToLeft(const RIGHT_FACTOR& f_right, G2& msg, const REAL omega){
   //   MakeRightFactorUniform(f_right, msg, omega);
   //}

   template<typename LEFT_FACTOR, typename G3>
   void SendMessageToRight(const LEFT_FACTOR& f_left, G3& msg, const REAL omega){
      //std::cout << "pairwise to sum pairwise weight = " << omega << "\n";
      MakeLeftFactorUniform(f_left, msg, omega);
   }

   //template<typename LEFT_FACTOR, typename G3>
   //void ReceiveMessageFromLeft(const LEFT_FACTOR& f_left, G3& msg){
   //   MakeLeftFactorUniform(f_left, msg, 1.0);
   //}

   template<typename LEFT_FACTOR, typename MSG>
   void MakeLeftFactorUniform(const LEFT_FACTOR& f_left, MSG& msg, const REAL omega){
      msg -= omega*f_left;
   }

   template<typename RIGHT_FACTOR, typename MSG>
   void MakeRightFactorUniform(const RIGHT_FACTOR& f_right, MSG& msg, const REAL omega){
      matrix<REAL> msgs(f_right.no_labels(), f_right.no_labels());
      f_right.marginalize_pairwise(msgs);
      if(transpose_) {
         msgs.transpose();
      }
      msg -= omega*msgs;
   }

   template<typename LEFT_FACTOR, typename MSG>
   void RepamLeft(LEFT_FACTOR& f_left, const MSG& msg)
   {
      for(INDEX x1=0; x1<f_left.dim1(); ++x1) {
         for(INDEX x2=0; x2<f_left.dim2(); ++x2) {
            f_left(x1,x2) += normalize( msg(x1,x2) );
         }
      }
   }

   template<typename RIGHT_FACTOR, typename MSG>
   void RepamRight(RIGHT_FACTOR& f_right, const MSG& msg)
   {
      if(!transpose_) {
         for(INDEX x1=0; x1<f_right.no_labels(); ++x1) {
            for(INDEX x2=0; x2<f_right.no_labels(); ++x2) {
               f_right.reg(x1,x2) += normalize( msg(x1,x2) );
            }
         }
      } else {
         for(INDEX x1=0; x1<f_right.no_labels(); ++x1) {
            for(INDEX x2=0; x2<f_right.no_labels(); ++x2) {
               f_right.reg(x1,x2) += normalize( msg(x2,x1) );
            }
         }
      }
   }

private:
   const bool transpose_; 
};

} // end namespace LP_MP

#endif // LP_MP_DISCRETE_TOMOGRAPHY_SEQUENTIAL_HXX
