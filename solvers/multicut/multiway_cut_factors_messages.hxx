#ifndef LP_MP_MULTIWAY_CUT_FACTORS_MESSAGES_HXX
#define LP_MP_MULTIWAY_CUT_FACTORS_MESSAGES_HXX

#include "factors/labeling_list_factor.hxx"
#include "multicut_factors_messages.hxx"
#include "factors/simplex_factor.hxx"
#include "vector.hxx"
#include "help_functions.hxx"

namespace LP_MP {

 
// for each node v in the multicut graph, exactly one edge from the terminals must be non-cut
// corresponds to unary factor in Potts MRF
class one_terminal_edge_active_factor : public vector<REAL> {
public:
   one_terminal_edge_active_factor(const INDEX no_edge) 
      : vector<REAL>(no_edge, 0.0)
   {
      assert(no_edge >= 2);
   }

   REAL LowerBound() const
   {
      REAL sum = 0.0;
      REAL max_elem = -std::numeric_limits<REAL>::infinity();
      for(auto it=this->begin(); it!=this->end(); ++it) {
         sum += *it;
         max_elem = std::max(max_elem, *it);
      }
      return sum - max_elem;
   }

   REAL EvaluatePrimal() const
   {
      const REAL sum = std::accumulate(this->begin(), this->end(), 0.0);
      return sum - (*this)[non_cut_edge_]; 
   }

   void MaximizePotentialAndComputePrimal()
   {
      non_cut_edge_ = std::max_element(this->begin(), this->end()) - this->begin();
   }

   void init_primal() { non_cut_edge_ = std::numeric_limits<INDEX>::max(); }
   INDEX& primal() { return non_cut_edge_; }
   INDEX primal() const { return non_cut_edge_; }

   template<typename ARCHIVE> void serialize_dual(ARCHIVE& ar) { ar( *static_cast<vector<REAL>*>(this) ); }
   template<typename ARCHIVE> void serialize_primal(ARCHIVE& ar) { ar( non_cut_edge_ ); }

private:
   INDEX non_cut_edge_;
};

// message between multicut edge and one_terminal_edge_active_factor defined above
// to do: obsolete, remove
class edge_terminal_constraint_message {
public:
   edge_terminal_constraint_message(const INDEX idx) : right_index_(idx) {}
   template<typename LEFT_FACTOR, typename MSG>
   void RepamLeft(LEFT_FACTOR& l, const MSG& msg)
   {
      l[0] += msg;
   }
   template<typename RIGHT_FACTOR, typename MSG>
   void RepamRight(RIGHT_FACTOR& r, const MSG& msg)
   {
      r[right_index_] += msg;
   }

   template<typename RIGHT_FACTOR, typename MSG>
   void ReceiveMessageFromRight(const RIGHT_FACTOR& r, MSG& msg)
   {
      const REAL sum = std::accumulate(r.begin(), r.end(), 0.0);
      const auto largest = two_largest_elements<REAL>(r.begin(), r.end());
      const REAL min_cost = sum - largest[0];
      const REAL second_min_cost = sum - largest[1];
      if(r[right_index_] == largest[0]) {
         msg -= largest[0] - largest[1]; //min_cost - second_min_cost;
      } else {
         msg -= r[right_index_] - largest[0];
      } 
   }

   template<typename LEFT_FACTOR, typename MSG>
   void SendMessageToRight(const LEFT_FACTOR& l, MSG& msg, const REAL omega)
   {
      const REAL msg_val = omega*l[0];
      msg -= msg_val;
   }

   template<typename LEFT_FACTOR, typename RIGHT_FACTOR>
   void ComputeLeftFromRightPrimal(LEFT_FACTOR& l, const RIGHT_FACTOR& r) 
   {
      if(r.primal() == right_index_) {
         l.primal()[0] = false;
      } else {
         l.primal()[0] = true; 
      } 
   }

private:
   const INDEX right_index_; // entry in one_temrinal_edge_active_factor. first entry refers to multicut edge
};

// the convex hull of a given edge uv and all terminal edges tv and tu for t in T (the terminal node set).
// Corresponds to the pairwise factor in graphical models
class multi_terminal_factor : public pairwise_potts_factor {
public:
  multi_terminal_factor(const INDEX T) : pairwise_potts_factor(T, 0.0) {} 

  REAL LowerBound() const {
     cut_to_one_hot_encoding(msg1_begin(), msg1_end());
     cut_to_one_hot_encoding(msg2_begin(), msg2_end());
     const REAL lb = pairwise_potts_factor::LowerBound();
     one_hot_to_cut_encoding(msg1_begin(), msg1_end());
     one_hot_to_cut_encoding(msg2_begin(), msg2_end());
     return lb;
  }

  REAL EvaluatePrimal() const 
  {
     if(this->primal_[0] < this->dim() && this->primal_[1] < this->dim()) {
        const REAL sum_u = std::accumulate(this->msg1_begin(), this->msg1_end(), 0.0);
        const REAL sum_v = std::accumulate(this->msg2_begin(), this->msg2_end(), 0.0);

        const REAL local_costs = sum_u - *(this->msg1_begin() + primal_[0]) + sum_v - *(this->msg2_begin() + primal_[1]);
        if(primal_[0] == primal_[1]) {
           return local_costs;
        } else {
           return local_costs + this->diff_cost();
        }
     } else {
        return std::numeric_limits<REAL>::infinity();
     }
  }

  template<typename VECTOR>
  void min_marginal_1(VECTOR& m) const
  {
     cut_to_one_hot_encoding(msg1_begin(), msg1_end());
     cut_to_one_hot_encoding(msg2_begin(), msg2_end());
     pairwise_potts_factor::min_marginal_1(m);
     one_hot_to_cut_encoding(msg1_begin(), msg1_end());
     one_hot_to_cut_encoding(msg2_begin(), msg2_end());
  }

  template<typename VECTOR>
  void min_marginal_2(VECTOR& m) const
  {
     cut_to_one_hot_encoding(msg1_begin(), msg1_end());
     cut_to_one_hot_encoding(msg2_begin(), msg2_end());
     pairwise_potts_factor::min_marginal_2(m);
     one_hot_to_cut_encoding(msg1_begin(), msg1_end());
     one_hot_to_cut_encoding(msg2_begin(), msg2_end());
  }

  REAL min_marginal_cut() const
  {
     const auto m = this->min_values();
     return m[1] - m[0];
  }
  
protected:

  // transform costs
  template<typename ITERATOR>
  static void cut_to_one_hot_encoding(ITERATOR begin, ITERATOR end)
  {
     const REAL sum = std::accumulate(begin, end, 0.0);
     std::transform(begin, end, begin, [sum](const REAL x) { return sum - x; });
     //for(; begin != end; ++begin) { *begin = sum - *begin; }
  }

  template<typename ITERATOR>
  static void one_hot_to_cut_encoding(ITERATOR begin, ITERATOR end)
  {
     const INDEX n = std::distance(begin, end);
     assert(n > 1);
     const REAL sum = std::accumulate(begin, end, 0.0);
     std::transform(begin, end, begin, [n,sum](const REAL x) { return -x + sum/(REAL(n-1)); }); 
  }
};


template<INDEX POSITION>
class at_most_one_active_multi_terminal_message {
public:
  ~at_most_one_active_multi_terminal_message()
  { 
    static_assert(POSITION == 0 || POSITION == 1,"");
  }

  template<typename LEFT_FACTOR, typename MSG>
    void RepamLeft(LEFT_FACTOR& l, const MSG& msg)
    {
      for(INDEX i=0; i<l.size(); ++i) {
        l[i] += msg[i];
      }
    }
  template<typename RIGHT_FACTOR, typename MSG>
    void RepamRight(RIGHT_FACTOR& r, const MSG& msg)
    {  
      if(POSITION == 0) {
        for(INDEX i=0; i<r.dim(); ++i) {
          *(r.msg1_begin() + i) += msg[i];
        }
      } else {
        for(INDEX i=0; i<r.dim(); ++i) {
          *(r.msg2_begin() + i) += msg[i];
        }
      }
    }

  template<typename RIGHT_FACTOR, typename MSG>
    void ReceiveMessageFromRight(const RIGHT_FACTOR& r, MSG& msg)
    {
      // perform fast belief propagation similary as for Potts MRF.
      vector<REAL> m(r.dim());
      if(POSITION == 0) {
        r.min_marginal_1(m);
      } else {
        r.min_marginal_2(m); 
      }
      msg -= m;
    }

  template<typename LEFT_FACTOR, typename MSG>
    void SendMessageToRight(const LEFT_FACTOR& l, MSG& msg, const REAL omega)
    {
      msg -= omega*l;
    }

  template<typename LEFT_FACTOR, typename RIGHT_FACTOR>
    void ComputeRightFromLeftPrimal(const LEFT_FACTOR& l, RIGHT_FACTOR& r)
    {
      r.primal()[POSITION] = l.primal();
    }
};

// do zrobienia> wypelnij tu
// left factor is multicut edge, right one is multi_terminal_factor
class multicut_edge_potts_message
{
public:
  template<typename LEFT_FACTOR, typename MSG>
    void RepamLeft(LEFT_FACTOR& l, const MSG& msg)
    {
      l[0] += msg;
    }

  template<typename RIGHT_FACTOR, typename MSG>
    void RepamRight(RIGHT_FACTOR& r, const MSG& msg)
    {  
      r.diff_cost() += msg;
    }

  
  template<typename RIGHT_FACTOR, typename MSG>
    void ReceiveMessageFromRight(const RIGHT_FACTOR& r, MSG& msg)
    {
       msg -= r.min_marginal_cut();
    }
    

  template<typename LEFT_FACTOR, typename MSG>
    void ReceiveMessageFromLeft(const LEFT_FACTOR& l, MSG& msg)
    {
       msg -= l[0];
    }

  
  template<typename LEFT_FACTOR, typename MSG>
    void SendMessageToRight(const LEFT_FACTOR& l, MSG& msg, const REAL omega)
    {
      msg -= omega*l[0];
    }
    

  template<typename RIGHT_FACTOR, typename MSG>
    void SendMessageToLeft(const RIGHT_FACTOR& r, MSG& msg, const REAL omega)
    {
       msg -= omega*r.min_marginal_cut();
    }


};

} // end namespace LP_MP

#endif // LP_MP_MULTIWAY_CUT_FACTORS_MESSAGES_HXX
