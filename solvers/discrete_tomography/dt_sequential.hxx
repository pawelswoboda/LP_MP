#ifndef LP_MP_DISCRETE_TOMOGRAPHY_SEQUENTIAL_HXX
#define LP_MP_DISCRETE_TOMOGRAPHY_SEQUENTIAL_HXX

#include "config.hxx"

namespace LP_MP {

// class holding state of variable and sum of chain up to that variable
class dt_sum_state_factor {
public:
   dt_sum_state_factor(const INDEX no_labels, const INDEX sum_size) 
      : pot_(no_labels, sum_size, 0.0)
   {}

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
      REAL min_value = std::numeric_limits<REAL>::infinity();
      for(auto it=pot_.begin(); it!=pot_.end(); ++it) {
         min_value = std::min(min_value, *it);
      }
      return min_value;
   }

   REAL EvaluatePrimal(PrimalSolutionStorage::Element primal) const {
      return std::numeric_limits<REAL>::infinity();
   }

   REAL& operator()(const INDEX state, const INDEX sum) {
      return pot_(state,sum);
   }
   const matrix& pot() const { return pot_; }

   INDEX no_labels() const { return pot_.dim1(); }
   INDEX sum_size() const { return pot_.dim2(); }
   INDEX size() const { return pot_.size(); }
private:
   matrix pot_; // first dimension is state, second one is sum
};

class dt_sum_state_pairwise_factor {
public:
   dt_sum_state_pairwise_factor(const INDEX no_labels, const INDEX prev_sum_size, const INDEX next_sum_size) 
      : prev_(no_labels, prev_sum_size, 0.0),
      next_(no_labels, next_sum_size, 0.0),
      reg_(no_labels, no_labels, 0.0)
   {
      assert(prev_sum_size <= next_sum_size);
   }

   REAL eval(const INDEX x1, const INDEX x2, const INDEX sum) const {
      assert(x1 < no_labels() && x2 < no_labels() && sum < prev_.dim2());
      return prev_(x1,sum) + reg_(x1,x2) + next_(x2, sum+x1);
   }
   REAL LowerBound() const {
      REAL min_value = std::numeric_limits<REAL>::infinity();
      for(INDEX x1=0; x1<no_labels(); ++x1) {
         for(INDEX x2=0; x2<no_labels(); ++x2) {
            for(INDEX sum=0; sum<prev_.dim2(); ++sum) {
               min_value = std::min(min_value, eval(x1,x2,sum));
            }
         }
      }
      return min_value;
   }

   REAL EvaluatePrimal(PrimalSolutionStorage::Element primal) const {
      return std::numeric_limits<REAL>::infinity();
   }

   void marginalize_pairwise(matrix& msg) const {
      std::fill(msg.begin(), msg.end(), std::numeric_limits<REAL>::infinity());

      for(INDEX x1=0; x1<no_labels(); ++x1) {
         for(INDEX x2=0; x2<no_labels(); ++x2) {
            for(INDEX sum=0; sum<prev_.dim2(); ++sum) {
               msg(x1,x2) = std::min(msg(x1,x2), eval(x1,x2,sum));
            }
         }
      }
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
private:
   matrix prev_; // first dimension is state, second one is sum
   matrix next_; // first dimension is state, second one is sum
   matrix reg_; // regularizer
};

// left factor is dt_sum_state_factor, right one is dt_sum_state_pairwise_factor
// DIRECTION signifies: left = sum_state_factor is previous
//                      right = sum_state_factor is next
// possibly make functions static
template<Chirality DIRECTION>
class dt_sum_pairwise_message {
public:
   dt_sum_pairwise_message(bool transpose) : transpose_(transpose) {
      assert(transpose_ == false);
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
      matrix msgs(f_right.no_labels(), DIRECTION == Chirality::left ? f_right.prev_sum_size() : f_right.next_sum_size(), std::numeric_limits<REAL>::infinity());
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
private:
   bool transpose_; 
};

// connect sequential dt_sum_factor (left) with recursive DiscreteTomograhpyCountingFactor (right)
template<Chirality DIRECTION>
class dt_sum_counting_message {
public:
   dt_sum_counting_message(const bool transpose) : transpose_(transpose) {}

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
      MakeLeftFactorUniform(f_left, msg, 0.5*omega);
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
      tensor3 msgs(f_right.no_labels(), f_right.no_labels(), DIRECTION == Chirality::left ? f_right.left_sum_size() : f_right.right_sum_size());
      if(DIRECTION == Chirality::left) {
         f_right.MessageCalculation_Left(msgs);
      } else {
         f_right.MessageCalculation_Right(msgs);
      }
      msg -= omega*msgs;
   }

   template<typename LEFT_FACTOR, typename MSG>
   void RepamLeft(LEFT_FACTOR& f_left, const MSG& msg)
   {
      //static_if<std::is_base_of<matrix,msg>([&](auto f) {});
      matrix msg_marg(f_left.no_labels(),f_left.sum_size(), std::numeric_limits<REAL>::infinity());
      for(INDEX x1=0; x1<f_left.no_labels(); ++x1) {
         for(INDEX x2=0; x2<f_left.no_labels(); ++x2) {
            for(INDEX sum=0; sum<f_left.sum_size(); ++sum) {
               if(DIRECTION == Chirality::left) {
                  msg_marg(x1,sum) += std::min( msg_marg(x1,sum), normalize( msg(x1,x2,sum) ));
               } else {
                  msg_marg(x1,sum) += std::min( msg_marg(x1,sum), normalize( msg(x1,x2,sum) )); 
               }
            }
         }
      }
      for(INDEX x=0; x<f_left.no_labels(); ++x) {
         for(INDEX sum=0; sum<f_left.sum_size(); ++sum) {
            f_left(x,sum) += msg_marg(x,sum);
         }
      }
   }

   template<typename RIGHT_FACTOR, typename MSG>
   void RepamRight(RIGHT_FACTOR& f_right, const MSG& msg)
   {
      /*
      for(INDEX x1=0; x1<f_right.no_labels(); ++x1) {
         for(INDEX x2=0; x2<f_right.no_labels(); ++x2) {
            if(DIRECTION == Chirality::left) {
               for(INDEX sum=0; sum<f_right.left_sum_size(); ++sum) {
                  f_right.left()(x1, x2, sum) += normalize( msg(x1, sum) );
               }
            } else {
               for(INDEX sum=0; sum<f_right.right_sum_size(); ++sum) {
                  f_right.right()(x1, x2, sum) += normalize( msg(x2, sum) );
               }
            }
         }
      }
      */
   }
private:
   bool transpose_; 
};

// connect mrf pairwise factor (left) with dt_sequential_pairwise_factor (right)
class dt_sum_pairwise_pairwise_message {
public:
   dt_sum_pairwise_pairwise_message (const bool transpose) : transpose_(transpose) {}

   //template<typename RIGHT_FACTOR, typename G2>
   //void ReceiveMessageFromRight(const RIGHT_FACTOR& f_right, G2& msg){
   //   MakeRightFactorUniform(f_right, msg, 1.0);
   //}

   template<typename RIGHT_FACTOR, typename G2>
   void SendMessageToLeft(const RIGHT_FACTOR& f_right, G2& msg, const REAL omega){
      MakeRightFactorUniform(f_right, msg, 0.5*omega);
   }

   //template<typename LEFT_FACTOR, typename G3>
   //void SendMessageToRight(const LEFT_FACTOR& f_left, G3& msg, const REAL omega){
   //   MakeLeftFactorUniform(f_left, msg, omega);
   //}

   template<typename LEFT_FACTOR, typename G3>
   void ReceiveMessageFromLeft(const LEFT_FACTOR& f_left, G3& msg){
      MakeLeftFactorUniform(f_left, msg, 1.0);
   }

   template<typename LEFT_FACTOR, typename MSG>
   void MakeLeftFactorUniform(const LEFT_FACTOR& f_left, MSG& msg, const REAL omega){
      msg -= omega*f_left;
   }

   template<typename RIGHT_FACTOR, typename MSG>
   void MakeRightFactorUniform(const RIGHT_FACTOR& f_right, MSG& msg, const REAL omega){
      matrix msgs(f_right.no_labels(), f_right.no_labels());
      f_right.marginalize_pairwise(msgs);
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
      for(INDEX x1=0; x1<f_right.no_labels(); ++x1) {
         for(INDEX x2=0; x2<f_right.no_labels(); ++x2) {
            f_right.reg(x1,x2) += normalize( msg(x1,x2) );
         }
      }
   }

private:
   bool transpose_; 
};

} // end namespace LP_MP

#endif // LP_MP_DISCRETE_TOMOGRAPHY_SEQUENTIAL_HXX
