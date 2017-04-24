#ifndef LP_MP_MULTIWAY_CUT_FACTORS_MESSAGES_HXX
#define LP_MP_MULTIWAY_CUT_FACTORS_MESSAGES_HXX

#include "factors/labeling_list_factor.hxx"
#include "multicut_factors_messages.hxx"
#include "vector.hxx"
#include "help_functions.hxx"

namespace LP_MP {

// for each node v in the multicut graph, exactly one edge from the terminals must be non-cut
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

// message between edge and one_temrinal_edge_active_factor defined above
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

// the convex hull of a given edge uv and all terminal edges tv and tu for t in T (the terminal node set)
class multi_terminal_factor : public vector<REAL> {
public:
  multi_terminal_factor(const INDEX T) : vector<REAL>(2*T+1) {} // vector stores cost of all terminal edges attached to u and v and of the edge uv

  INDEX no_classes() const { return this->size()/2; }

  // min cost for same label and different label
  std::array<REAL,2> min_values() const
  {
    const REAL sum_u = std::accumulate(this->begin(), this->begin() + no_classes(), 0.0);
    const REAL sum_v = std::accumulate(this->begin() + no_classes(), this->begin() + 2*no_classes(), 0.0);

    // v holds the unary costs for one hot encoding (e.g. as in MRF)
    vector<REAL> v(2*no_classes());
    for(INDEX i=0; i<no_classes(); ++i) {
      v[i] = sum_u - (*this)[i];
    }
    for(INDEX i=0; i<no_classes(); ++i) {
      v[i+no_classes()] = sum_v - (*this)[i+no_classes()];
    }
    const auto smallest2 = two_smallest_elements<REAL>(v.begin() + no_classes(), v.begin() + 2*no_classes());

    REAL min_same_label = std::numeric_limits<REAL>::infinity();
    REAL min_diff_label = std::numeric_limits<REAL>::infinity();
    for(INDEX i=0; i<no_classes(); ++i) {
      const REAL same_label = v[i] + v[i+no_classes()];
      min_same_label = std::min(min_same_label, same_label);
      const REAL diff_label = v[i] + v[2*no_classes()] + v[i+no_classes()] == smallest2[0] ? smallest2[1] : smallest2[0];
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
    const REAL sum_u = std::accumulate(this->begin(), this->begin() + no_classes(), 0.0);
    const REAL sum_v = std::accumulate(this->begin() + no_classes(), this->begin() + 2*no_classes(), 0.0);

    const REAL local_costs = sum_u - (*this)[label_[0]] + sum_v - (*this)[no_classes() + label_[1]];
    if(label_[0] == label_[1]) {
      return local_costs;
    } else {
      return local_costs + (*this)[2*no_classes()];
    }
  }

  template<typename VECTOR>
  void min_marginal_1(VECTOR& m) const
  {
    assert(m.size() == no_classes());
    const REAL sum_u = std::accumulate(this->begin(), this->begin() + no_classes(), 0.0);
    const REAL sum_v = std::accumulate(this->begin() + no_classes(), this->begin() + 2*no_classes(), 0.0);

    // v holds the unary costs for one hot encoding (e.g. as in MRF)
    vector<REAL> v(2*no_classes());
    for(INDEX i=0; i<no_classes(); ++i) {
      v[i] = sum_u - (*this)[i];
    }
    for(INDEX i=0; i<no_classes(); ++i) {
      v[i+no_classes()] = sum_v - (*this)[i+no_classes()];
    }
    const auto smallest2 = two_smallest_elements<REAL>(v.begin() + no_classes(), v.begin() + 2*no_classes());

    for(INDEX i=0; i<no_classes(); ++i) {
      const REAL same_label = v[i+no_classes()];
      const REAL diff_label = v[2*no_classes()] + v[i+no_classes()] == smallest2[0] ? smallest2[1] : smallest2[0];
      m[i] = v[i] + std::min(same_label, diff_label); 
    } 

    // transform back from one hot encoding to all except one encoding
    for(INDEX i=0; i<no_classes(); ++i) {

    }
  }

  template<typename VECTOR>
  void min_marginal_2(VECTOR& m) const
  {
    assert(m.size() == no_classes());
    const REAL sum_u = std::accumulate(this->begin(), this->begin() + no_classes(), 0.0);
    const REAL sum_v = std::accumulate(this->begin() + no_classes(), this->begin() + 2*no_classes(), 0.0);

    // v holds the unary costs for one hot encoding (e.g. as in MRF)
    vector<REAL> v(2*no_classes());
    for(INDEX i=0; i<no_classes(); ++i) {
      v[i] = sum_u - (*this)[i];
    }
    for(INDEX i=0; i<no_classes(); ++i) {
      v[i+no_classes()] = sum_v - (*this)[i+no_classes()];
    }
    const auto smallest2 = two_smallest_elements<REAL>(v.begin(), v.begin() + no_classes());

    for(INDEX i=0; i<no_classes(); ++i) {
      const REAL same_label = v[i];
      const REAL diff_label = v[2*no_classes()] + v[i] == smallest2[0] ? smallest2[1] : smallest2[0];
      m[i] = v[i+no_classes()] + std::min(same_label, diff_label); 
    } 
  }
  
  void init_primal() { label_[0] = std::numeric_limits<INDEX>::max(); label_[1] = std::numeric_limits<INDEX>::max(); }
  auto& primal() { return label_; }
  const auto& primal() const { return label_; }

  template<typename ARCHIVE> void serialize_dual(ARCHIVE& ar) { ar( *static_cast<vector<REAL>*>(this) ); }
  template<typename ARCHIVE> void serialize_primal(ARCHIVE& ar) { ar( label_ ); }

private:
  std::array<INDEX,2> label_;
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
        for(INDEX i=0; i<r.no_classes(); ++i) {
          r[i] += msg[i];
        }
      } else {
        for(INDEX i=0; i<r.no_classes(); ++i) {
          r[i + r.no_classes()] += msg[i];
        }
      }
    }

  template<typename RIGHT_FACTOR, typename MSG>
    void ReceiveMessageFromRight(const RIGHT_FACTOR& r, MSG& msg)
    {
      // perform fast belief propagation similary as for Potts MRF.
      vector<REAL> m(r.no_classes());
      if(POSITION == 0) {
        r.min_marginal_1(m);
      } else {
        r.min_marginal_1(m); 
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

class multicut_edge_multi_terminal_message
{
public:
  template<typename LEFT_FACTOR, typename MSG>
    void RepamLeft(LEFT_FACTOR& l, const MSG& msg)
    {
      l += msg;
    }
  template<typename RIGHT_FACTOR, typename MSG>
    void RepamRight(RIGHT_FACTOR& r, const MSG& msg)
    {  
      r[r.no_classes()] += msg;
    }

  template<typename RIGHT_FACTOR, typename MSG>
    void ReceiveMessageFromRight(const RIGHT_FACTOR& r, MSG& msg)
    {
      // perform fast belief propagation similary as for Potts MRF.
    }

  template<typename LEFT_FACTOR, typename MSG>
    void SendMessageToRight(const LEFT_FACTOR& l, MSG& msg, const REAL omega)
    {
      msg -= omega*l;
    }

};


// first edge is multicut edge, then come edges to terminal node
// it is only required that whenever multicut edge is cut, then at least one terminal edge must be cut as well.
using asymmetric_multiway_cut_triplet_labelings = labelings<
  labeling<0,0,1>,
  labeling<0,1,0>,
  labeling<0,1,1>,
  labeling<1,0,1>,
  labeling<1,1,0>,
  labeling<1,1,1>
   >;

using asymmetric_multiway_cut_triplet_factor = labeling_factor< asymmetric_multiway_cut_triplet_labelings, true >;
using edge_asymmetric_multiway_cut_triplet_message_0 = labeling_message< multicut_edge_labelings, asymmetric_multiway_cut_triplet_labelings, 0 >;
using edge_asymmetric_multiway_cut_triplet_message_1 = labeling_message< multicut_edge_labelings, asymmetric_multiway_cut_triplet_labelings, 1 >;
using edge_asymmetric_multiway_cut_triplet_message_2 = labeling_message< multicut_edge_labelings, asymmetric_multiway_cut_triplet_labelings, 2 >;

} // end namespace LP_MP

#endif // LP_MP_MULTIWAY_CUT_FACTORS_MESSAGES_HXX
