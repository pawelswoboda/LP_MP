#ifndef LP_MP_DISCRETE_TOMOGRAPHY_SEQUENTIAL_HXX
#define LP_MP_DISCRETE_TOMOGRAPHY_SEQUENTIAL_HXX

#include "config.hxx"

namespace LP_MP {

// class holding state of variable and sum of chain up to that variable
/*
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
*/

class dt_sum_state_pairwise_factor {
public:
   constexpr static INDEX no_primal_decision = std::numeric_limits<INDEX>::max();

   dt_sum_state_pairwise_factor(const INDEX no_labels, const INDEX prev_sum_size, const INDEX next_sum_size) 
      : prev_(no_labels, prev_sum_size, 0.0),
      next_(no_labels, next_sum_size, 0.0),
      reg_(no_labels, no_labels, 0.0)
   {
      assert(this->prev_sum_size() <= this->next_sum_size());
      assert(this->next_sum_size() <= this->prev_sum_size() + this->no_labels());
      assert(this->next_sum_size() >= this->no_labels());
   }

   template<typename ITERATOR>
   void summation_cost(ITERATOR begin, ITERATOR end) {
      assert(std::distance(begin,end) <= next_sum_size()+no_labels());
      for(INDEX x2=0; x2<no_labels(); ++x2) {
         for(INDEX sum=0; sum<next_sum_size(); ++sum) {
            if(sum+x2 < std::distance(begin,end)) {
               next_(x2,sum) = *(begin + sum+x2);
            } else {
               next_(x2,sum) = std::numeric_limits<REAL>::infinity(); 
            }
         }
      }
   }

   template<typename LAMBDA>
   void for_each_label_sum(const LAMBDA&& f) const {
      assert(prev_sum_size() <= next_sum_size());
      assert(prev_sum_size() + no_labels()-1 >= next_sum_size());
      for(INDEX x1=0; x1<no_labels(); ++x1) {
         const INDEX max_sum = std::min(prev_sum_size(), next_sum_size()-x1);
         for(INDEX x2=0; x2<no_labels(); ++x2) {
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
      if( !(state_[0] < no_labels() && state_[1] < no_labels()) ) {
         return std::numeric_limits<REAL>::infinity();
      }
      if( !(sum_[0] + state_[0] == sum_[1]) ) {
         return std::numeric_limits<REAL>::infinity();
      }
      if( !(sum_[0] < prev_sum_size() && sum_[1] < next_sum_size()) ) {
         return std::numeric_limits<REAL>::infinity();
      }
      return eval(state_[0], state_[1], sum_[0]);
   }

   void marginalize_pairwise(matrix<REAL>& msg) const {
      std::fill(msg.begin(), msg.end(), std::numeric_limits<REAL>::infinity());
      for_each_label_sum([&](INDEX x1, INDEX x2, INDEX sum) { msg(x1,x2) = std::min(msg(x1,x2), eval(x1,x2,sum)); });

      //std::cout << "(" << no_labels() << "," << prev_sum_size() << "," << next_sum_size() <<")\n";
      //for(INDEX x1=0; x1<no_labels(); ++x1) {
      //   for(INDEX x2=0; x2<no_labels(); ++x2) {
      //      std::cout << msg(x1,x2) << ",";
      //   }
      //std::cout << "\n";
      //}
      //std::cout << "\n";
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



   template<typename SAT_SOLVER>
   void construct_sat_clauses(SAT_SOLVER& s) const
   {
      auto prev_vars = create_sat_variables(s, no_labels()*prev_sum_size());
      auto next_vars = create_sat_variables(s, no_labels()*next_sum_size());
      auto pairwise_vars = create_sat_variables(s, no_labels()*no_labels());
      auto gluing_vars = create_sat_variables(s, no_labels()*no_labels()*prev_sum_size());

      add_simplex_constraint_sat(s, prev_vars.begin(), prev_vars.end());
      add_simplex_constraint_sat(s, next_vars.begin(), next_vars.end());
      // pairwise vars sum to one through message to pairwise simplex factor
      //add_simplex_constraint_sat(s, pairwise_vars.begin(), pairwise_vars.end());
      // not needed, forced by simplex constraints on marginalized variables
      //add_simplex_constraint_sat(s, gluing_vars.begin(), gluing_vars.end()); 

      std::vector<sat_var> tmp_vars;
      tmp_vars.reserve(no_labels()*prev_sum_size());

      for(INDEX x1=0; x1<no_labels(); ++x1) {
         for(INDEX sum=0; sum<prev_sum_size(); ++sum) {
            for(INDEX x2=0; x2<no_labels(); ++x2) {
               tmp_vars.push_back(gluing_vars[x1*no_labels()*prev_sum_size() + x2*prev_sum_size() + sum]);
            }
            const auto prev_var = prev_vars[x1*prev_sum_size() + sum];
            const auto sum_var = add_at_most_one_constraint_sat(s, tmp_vars.begin(), tmp_vars.end());
            make_sat_var_equal(s, to_literal(prev_var), to_literal(sum_var));
            tmp_vars.clear();
         }
      }

      for(INDEX x2=0; x2<no_labels(); ++x2) {
         for(INDEX next_sum=0; next_sum<next_sum_size(); ++next_sum) {
            for(INDEX x1=0; x1<no_labels(); ++x1) {
               if(x1 <= next_sum  && next_sum-x1 < prev_sum_size()) {
                  tmp_vars.push_back(gluing_vars[x1*no_labels()*prev_sum_size() + x2*prev_sum_size() + (next_sum-x1)]);
               }
            }
            const auto next_var = next_vars[x2*next_sum_size() + next_sum];
            const auto sum_var = add_at_most_one_constraint_sat(s, tmp_vars.begin(), tmp_vars.end());
            if(!tmp_vars.empty()) {
               make_sat_var_equal(s, to_literal(next_var), to_literal(sum_var));
            }
            tmp_vars.clear();
         }
      }

      for(INDEX x1=0; x1<no_labels(); ++x1) {
         for(INDEX x2=0; x2<no_labels(); ++x2) {
            for(INDEX sum=0; sum<prev_sum_size(); ++sum) {
               tmp_vars.push_back(gluing_vars[x1*no_labels()*prev_sum_size() + x2*prev_sum_size() + sum]);
            }
            const auto pairwise_var = pairwise_vars[x1*no_labels() + x2];
            const auto sum_var = add_at_most_one_constraint_sat(s, tmp_vars.begin(), tmp_vars.end());
            make_sat_var_equal(s, to_literal(pairwise_var), to_literal(sum_var));
            tmp_vars.clear(); 
         }
      } 
   }

   template<typename VEC>
   void reduce_sat(VEC& assumptions, const REAL th, sat_var begin) const
   {
      const REAL lb = LowerBound();
      const sat_var gluing_vars_begin = begin + no_labels()*prev_sum_size() + no_labels()*next_sum_size() + no_labels()*no_labels();
      INDEX no_active_labels = 0;
      for_each_label_sum([this,lb,th,gluing_vars_begin,&assumptions,&no_active_labels](const INDEX x1, const INDEX x2, const INDEX sum) {
            if(this->eval(x1,x2,sum) > lb + th) {
              assumptions.push_back(-to_literal(gluing_vars_begin + x1*this->no_labels()*this->prev_sum_size() + x2*this->prev_sum_size() + sum));
            } else {
              no_active_labels++;
            } 
      });
      assert(no_active_labels > 0);
      if(no_active_labels == 1) { // possibly directly force label if only one active label survives
      }
   }

   template<typename SAT_SOLVER>
   void convert_primal(SAT_SOLVER& s, sat_var first)
   {
      for(INDEX x1=0; x1<no_labels(); ++x1) {
         for(INDEX sum=0; sum<prev_sum_size(); ++sum) {
            if(lglderef(s,to_literal(first + x1*prev_sum_size() + sum)) == 1) {
               state_[0] = x1;
               sum_[0] = sum;
            }
         }
      }

      for(INDEX x2=0; x2<no_labels(); ++x2) {
         for(INDEX sum=0; sum<next_sum_size(); ++sum) {
            if(lglderef(s,to_literal(first + no_labels()*prev_sum_size() + x2*next_sum_size() + sum)) == 1) {
               state_[1] = x2;
               sum_[1] = sum;
            }
         }
      }
      assert(state_[0] + sum_[0] == sum_[1]);
   }

   std::array<INDEX,2> state_;
   std::array<INDEX,2> sum_; // both sums are not strictly needed, they help however in labeling
private:
   matrix<REAL> prev_; // first dimension is state, second one is sum
   matrix<REAL> next_; // first dimension is state, second one is sum
   matrix<REAL> reg_; // regularizer
};


// message between two dt_sequential_pairwise_factors. Left factor is smaller one
class dt_pairwise_message {
public:
   template<typename RIGHT_FACTOR, typename G2>
   void ReceiveMessageFromRight(const RIGHT_FACTOR& f_right, G2& msg){
      MakeRightFactorUniform(f_right, msg, 01.0);
   }

   template<typename RIGHT_FACTOR, typename G2>
   void SendMessageToLeft(const RIGHT_FACTOR& f_right, G2& msg, const REAL omega){
      MakeRightFactorUniform(f_right, msg, 01.0*omega);
   }

   template<typename LEFT_FACTOR, typename G3>
   void SendMessageToRight(const LEFT_FACTOR& f_left, G3& msg, const REAL omega){
      //std::cout << "sum to sum pairwise weight = " << omega << "\n";
      MakeLeftFactorUniform(f_left, msg, 01.0*omega);
   }

   template<typename LEFT_FACTOR, typename G3>
   void ReceiveMessageFromLeft(const LEFT_FACTOR& f_left, G3& msg){
      MakeLeftFactorUniform(f_left, msg, 01.0);
   }

   template<typename LEFT_FACTOR, typename MSG>
   void MakeLeftFactorUniform(const LEFT_FACTOR& f_left, MSG& msg, const REAL omega)
   {
      matrix<REAL> msgs(f_left.no_labels(), f_left.next_sum_size(), std::numeric_limits<REAL>::infinity());
      f_left.for_each_label_sum([&msgs, f_left](const INDEX x1, const INDEX x2, const INDEX sum) {
            assert(sum+x1 < f_left.next_sum_size());
            msgs(x2,sum+x1) = std::min(msgs(x2,sum+x1), f_left.eval(x1,x2,sum));
      });
      msg -= omega*msgs;
   }

   template<typename RIGHT_FACTOR, typename MSG>
   void MakeRightFactorUniform(const RIGHT_FACTOR& f_right, MSG& msg, const REAL omega)
   {
      matrix<REAL> msgs(f_right.no_labels(), f_right.prev_sum_size(), std::numeric_limits<REAL>::infinity());
      f_right.for_each_label_sum([&msgs, f_right](const INDEX x1, const INDEX x2, const INDEX sum) {
            msgs(x1,sum) = std::min(msgs(x1,sum), f_right.eval(x1,x2,sum));
      });
      msg -= omega*msgs;
   }

   template<typename LEFT_FACTOR, typename MSG>
   void RepamLeft(LEFT_FACTOR& f_left, const MSG& msg)
   {
      assert(msg.dim1() == f_left.no_labels());
      assert(msg.dim2() == f_left.next_sum_size());
      for(INDEX x=0; x<f_left.no_labels(); ++x) {
         for(INDEX sum=0; sum<f_left.next_sum_size(); ++sum) {
            f_left.next(x,sum) += normalize( msg(x,sum) );
         }
      }
   }

   template<typename RIGHT_FACTOR, typename MSG>
   void RepamRight(RIGHT_FACTOR& f_right, const MSG& msg)
   {
      assert(msg.dim1() == f_right.no_labels());
      assert(msg.dim2() == f_right.prev_sum_size()); 
      for(INDEX x=0; x<f_right.no_labels(); ++x) {
         for(INDEX sum=0; sum<f_right.prev_sum_size(); ++sum) {
            f_right.prev(x,sum) += normalize( msg(x,sum) );
         }
      }
   }

   template<typename SAT_SOLVER, typename LEFT_FACTOR, typename RIGHT_FACTOR>
   void construct_sat_clauses(SAT_SOLVER& s, LEFT_FACTOR& l, RIGHT_FACTOR& r, sat_var left_begin, sat_var right_begin) const
   {
      assert(l.next_sum_size() == r.prev_sum_size());
      for(INDEX x=0; x<l.no_labels(); ++x) {
         for(INDEX sum=0; sum<l.next_sum_size(); ++sum) {
            const auto left_var = left_begin + l.no_labels()*l.prev_sum_size()  + x*l.next_sum_size() + sum;
            const auto right_var = right_begin + x*r.prev_sum_size() + sum;
            make_sat_var_equal(s, to_literal(left_var), to_literal(right_var));
         }
      }
   } 

};

// message from mrf unary to sequential pairwise factor
template<Chirality DIRECTION>
class dt_sum_unary_message {
public:
   template<typename RIGHT_FACTOR, typename G2>
   void ReceiveMessageFromRight(const RIGHT_FACTOR& f_right, G2& msg){
      MakeRightFactorUniform(f_right, msg, 01.0);
   }

   template<typename RIGHT_FACTOR, typename G2>
   void SendMessageToLeft(const RIGHT_FACTOR& f_right, G2& msg, const REAL omega){
      MakeRightFactorUniform(f_right, msg, 01.0*omega);
   }

   template<typename LEFT_FACTOR, typename G3>
   void SendMessageToRight(const LEFT_FACTOR& f_left, G3& msg, const REAL omega){
      MakeLeftFactorUniform(f_left, msg, 01.0*omega);
   }

   template<typename LEFT_FACTOR, typename G3>
   void ReceiveMessageFromLeft(const LEFT_FACTOR& f_left, G3& msg){
      MakeLeftFactorUniform(f_left, msg, 01.0);
   }

   template<typename LEFT_FACTOR, typename MSG>
   void MakeLeftFactorUniform(const LEFT_FACTOR& f_left, MSG& msg, const REAL omega)
   {
      msg -= omega*f_left;
   }

   template<typename RIGHT_FACTOR, typename MSG>
   void MakeRightFactorUniform(const RIGHT_FACTOR& f_right, MSG& msg, const REAL omega)
   {
      vector<REAL> msg_val(f_right.no_labels(), std::numeric_limits<REAL>::infinity());
      f_right.for_each_label_sum([&msg_val,&f_right](const INDEX x1, const INDEX x2, const INDEX sum) {
            if(DIRECTION == Chirality::left) {
              msg_val[x1] = std::min( msg_val[x1], f_right.eval(x1,x2,sum) ); 
            } else {
              msg_val[x2] = std::min( msg_val[x2], f_right.eval(x1,x2,sum) ); 
            }
      });
      msg -= omega*msg_val;
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
      for(INDEX x1=0; x1<f_right.no_labels(); ++x1) {
         for(INDEX x2=0; x2<f_right.no_labels(); ++x2) {
            if(DIRECTION == Chirality::left) {
               f_right.reg(x1,x2) += normalize( msg[x1] );
            } else {
               f_right.reg(x1,x2) += normalize( msg[x2] );
            }
         }
      }
   } 
   template<typename SAT_SOLVER, typename LEFT_FACTOR, typename RIGHT_FACTOR>
   void construct_sat_clauses(SAT_SOLVER& s, LEFT_FACTOR& l, RIGHT_FACTOR& r, sat_var left_begin, sat_var right_begin) const
   {
   }
};

// left factor is dt_sum_state_factor, right one is dt_sum_state_pairwise_factor
// DIRECTION signifies: left = sum_state_factor is previous
//                      right = sum_state_factor is next
// possibly make functions static
/*
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
      MakeRightFactorUniform(f_right, msg, 01.0);
   }

   template<typename RIGHT_FACTOR, typename G2>
   void SendMessageToLeft(const RIGHT_FACTOR& f_right, G2& msg, const REAL omega){
      MakeRightFactorUniform(f_right, msg, 01.0*omega);
   }

   template<typename LEFT_FACTOR, typename G3>
   void SendMessageToRight(const LEFT_FACTOR& f_left, G3& msg, const REAL omega){
      //std::cout << "sum to sum pairwise weight = " << omega << "\n";
      MakeLeftFactorUniform(f_left, msg, 01.0*omega);
   }

   template<typename LEFT_FACTOR, typename G3>
   void ReceiveMessageFromLeft(const LEFT_FACTOR& f_left, G3& msg){
      MakeLeftFactorUniform(f_left, msg, 01.0);
   }

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

   
   //template<typename LEFT_FACTOR, typename RIGHT_FACTOR>
   //void ComputeLeftFromRightPrimal(LEFT_FACTOR& l, const RIGHT_FACTOR& r) {
   //   if(DIRECTION == Chirality::left) {
   //      assert(r.state_[0] < r.no_labels() && r.sum_[0] < r.prev_sum_size());
   //      l.state_ = r.state_[0];
   //      l.sum_ = r.sum_[0];
   //   } else {
   //      assert(r.state_[1] < r.no_labels() && r.sum_[1] < r.prev_sum_size());
   //      l.state_ = r.state_[1]; 
   //      l.sum_ = r.sum_[1];
   //   } 
   //}
   
private:
   bool transpose_; // not needed
};
*/

// connect sequential dt_sum_state_pairwise_factor (left) with recursive DiscreteTomograhpyCountingFactor (right)
template<Chirality DIRECTION>
class dt_sum_counting_message {
public:
   ~dt_sum_counting_message() {
      static_assert(DIRECTION == Chirality::left || DIRECTION == Chirality::right, "");
   }

   template<typename RIGHT_FACTOR, typename G2>
   void ReceiveMessageFromRight(const RIGHT_FACTOR& f_right, G2& msg){
      MakeRightFactorUniform(f_right, msg, 01.0);
   }

   template<typename RIGHT_FACTOR, typename G2>
   void SendMessageToLeft(const RIGHT_FACTOR& f_right, G2& msg, const REAL omega){
      MakeRightFactorUniform(f_right, msg, 01.0*omega);
   }

   template<typename LEFT_FACTOR, typename G3>
   void SendMessageToRight(const LEFT_FACTOR& f_left, G3& msg, const REAL omega){
      MakeLeftFactorUniform(f_left, msg, 01.0*omega);
   }

   template<typename LEFT_FACTOR, typename G3>
   void ReceiveMessageFromLeft(const LEFT_FACTOR& f_left, G3& msg){
      MakeLeftFactorUniform(f_left, msg, 01.0);
   }

   template<typename LEFT_FACTOR, typename MSG>
   void MakeLeftFactorUniform(const LEFT_FACTOR& f_left, MSG& msg, const REAL omega){
      matrix<REAL> msg_val(f_left.no_labels(), f_left.next_sum_size(), std::numeric_limits<REAL>::infinity());
      f_left.for_each_label_sum([&msg_val,&f_left](const INDEX x1, const INDEX x2, const INDEX sum) {
            msg_val(x2, x1+sum) = std::min( msg_val(x2, x1+sum), f_left.eval(x1,x2,sum) ); 
      });
      msg -= omega*msg_val;
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
      assert(msg.dim2() == f_left.next_sum_size());
      for(INDEX x=0; x<f_left.no_labels(); ++x) {
         for(INDEX sum=0; sum<f_left.next_sum_size(); ++sum) {
            f_left.next(x,sum) += normalize( msg(x,sum) );
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
      MakeRightFactorUniform(f_right, msg, 01.0);
   }

   template<typename RIGHT_FACTOR, typename G2>
   void SendMessageToLeft(const RIGHT_FACTOR& f_right, G2& msg, const REAL omega){
      MakeRightFactorUniform(f_right, msg, 01.0*omega);
   }

   template<typename LEFT_FACTOR, typename G3>
   void SendMessageToRight(const LEFT_FACTOR& f_left, G3& msg, const REAL omega){
      //std::cout << "pairwise to sum pairwise weight = " << omega << "\n";
      MakeLeftFactorUniform(f_left, msg, 01.0*omega);
   }

   template<typename LEFT_FACTOR, typename G3>
   void ReceiveMessageFromLeft(const LEFT_FACTOR& f_left, G3& msg){
      MakeLeftFactorUniform(f_left, msg, 01.0);
   }

   template<typename LEFT_FACTOR, typename MSG>
   void MakeLeftFactorUniform(const LEFT_FACTOR& f_left, MSG& msg, const REAL omega){
      msg -= omega*f_left;
   }

   template<typename RIGHT_FACTOR, typename MSG>
   void MakeRightFactorUniform(const RIGHT_FACTOR& f_right, MSG& msg, const REAL omega){
      matrix<REAL> msgs(f_right.no_labels(), f_right.no_labels(), std::numeric_limits<REAL>::infinity());
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
            f_left.cost(x1,x2) += normalize( msg(x1,x2) );
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

   template<typename SAT_SOLVER, typename LEFT_FACTOR, typename RIGHT_FACTOR>
   void construct_sat_clauses(SAT_SOLVER& s, LEFT_FACTOR& l, RIGHT_FACTOR& r, sat_var left_begin, sat_var right_begin) const
   {
      assert(!transpose_);
      assert(r.no_labels() == l.dim1());
      assert(r.no_labels() == l.dim2());
      // to do: take care of transposed_!
      for(INDEX x1=0; x1<l.dim1(); ++x1) {
         for(INDEX x2=0; x2<l.dim2(); ++x2) {
            const sat_var left_var = left_begin + l.dim1() + l.dim2() + x1*l.dim1() + x2;
            const sat_var right_var = right_begin + r.no_labels()*r.prev_sum_size() + r.no_labels()*r.next_sum_size() + x1*r.no_labels() + x2;
            make_sat_var_equal(s, to_literal(left_var), to_literal(right_var));
         }
      }
   }

private:
   const bool transpose_; 
};

} // end namespace LP_MP

#endif // LP_MP_DISCRETE_TOMOGRAPHY_SEQUENTIAL_HXX
