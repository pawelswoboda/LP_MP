#ifndef LP_MP_DETECTION_FACTOR_HXX
#define LP_MP_DETECTION_FACTOR_HXX

#include "config.hxx"
#include "vector.hxx"
#include "cereal/types/array.hpp"
#include "sat_interface.hxx"

namespace LP_MP {

enum class exit_constraint_position {lower,upper}; // detection is non-overlapping below/above some other

// factor containing whether cell hypothesis is true, all its possible predecessors and successors.

// primal has the following meaning:
// incoming_edge_ = std::numeric_limits<INDEX>::max() means no decision is taken yet.
//                = std::numeric_limits<REAL>::max()-1 means that no primal is to be taken.
//                < no_incoming_edges() means incoming edge is assigned
// same for outgoing_edge_
class detection_factor {

  friend class transition_message;
  friend class at_most_one_cell_message;
  template<exit_constraint_position POSITION> friend class exit_constraint_message;

public:

  constexpr static INDEX no_primal_decision = std::numeric_limits<INDEX>::max();
  constexpr static INDEX no_edge_taken = std::numeric_limits<INDEX>::max()-1;

  detection_factor(const INDEX no_incoming_edges, const INDEX no_outgoing_edges, const REAL detection_cost, const REAL appearance_cost, const REAL disappearance_cost)
    : no_incoming_edges_(no_incoming_edges+1),
    no_outgoing_edges_(no_outgoing_edges+1)
  {
    pot_ = global_real_block_allocator_array[stack_allocator_index].allocate(size());
    assert(pot_ != nullptr);
    std::fill(pot_+1, pot_ + size(), 0.0);
    pot_[0] = detection_cost;
    incoming( no_incoming_edges ) = appearance_cost;
    outgoing( no_outgoing_edges ) = disappearance_cost;
  }
  ~detection_factor() {
    global_real_block_allocator_array[stack_allocator_index].deallocate(pot_,1);
  }
  detection_factor(const detection_factor& o) 
    : no_incoming_edges_(o.no_incoming_edges_),
    no_outgoing_edges_(o.no_outgoing_edges_)
  {
    pot_ = global_real_block_allocator_array[stack_allocator_index].allocate(size());
    assert(pot_ != nullptr);
    for(INDEX i=0; i<size(); ++i) { pot_[i] = o.pot_[i]; }
  }
  void operator=(const detection_factor& o) {
    assert(no_incoming_edges_ == o.no_incoming_edges_);
    assert(no_outgoing_edges_ == o.no_outgoing_edges_);
    for(INDEX i=0; i<size(); ++i) { pot_[i] = o.pot_[i]; }
  }

  REAL operator[](const INDEX i) const {
    assert(i<size());
    return pot_[i];
  }
  REAL& operator[](const INDEX i) {
    assert(i<size());
    return pot_[i];
  }
  REAL incoming(const INDEX i) const {
    assert(i<no_incoming_edges_);
    return pot_[1+i];
  }
  REAL& incoming(const INDEX i) {
    assert(i<no_incoming_edges_);
    return pot_[1+i];
  }
  REAL outgoing(const INDEX i) const {
    assert(i<no_outgoing_edges_);
    return pot_[1+no_incoming_edges_+i];
  }
  REAL& outgoing(const INDEX i) {
    assert(i<no_outgoing_edges_);
    return pot_[1+no_incoming_edges_+i];
  }

  REAL cost_of_detection() const {
    //std::cout << pot_[0] << " : ";
    //for(auto it=incoming_begin(); it!= incoming_end(); ++it) std::cout << *it << ",";
    //std:cout << " : ";
    //for(auto it=outgoing_begin(); it!= outgoing_end(); ++it) std::cout << *it << ",";
    //std::cout << " = ";

    assert(no_incoming_edges_ > 0);
    assert(no_outgoing_edges_ > 0);
    REAL lb = pot_[0];
    lb += *std::min_element(incoming_begin(),incoming_end()); 
    lb += *std::min_element(outgoing_begin(),outgoing_end()); 
    //std::cout << lb << "\n";
    return lb;
  }

  void MaximizePotentialAndComputePrimal() {
    if(incoming_edge_ < no_incoming_edges() && outgoing_edge_ < no_outgoing_edges()) {
      return;
    }
    if(incoming_edge_ == no_edge_taken && no_incoming_edges() > 0) {
      outgoing_edge_ = no_edge_taken;
      return;
    }
    if(outgoing_edge_ == no_edge_taken && no_outgoing_edges() > 0) {
      incoming_edge_ = no_edge_taken;
      return;
    }

    INDEX incoming_cand = no_edge_taken;
    REAL lb = pot_[0];
    if(no_incoming_edges() > 0) {
      incoming_cand = std::min_element(incoming_begin(), incoming_end()) - incoming_begin();
      lb += incoming(incoming_cand); 
    }
    INDEX outgoing_cand = no_edge_taken;
    if(no_outgoing_edges() > 0) {
      outgoing_cand = std::min_element(outgoing_begin(), outgoing_end()) - outgoing_begin();
      lb += outgoing(outgoing_cand); 
    }
    // one edge already labelled
    if(incoming_edge_ < no_incoming_edges() && no_outgoing_edges() > 0 && outgoing_edge_ == no_primal_decision) {
      assert(outgoing_edge_ == no_primal_decision);
      outgoing_edge_ = outgoing_cand;
      return;
    } else if(outgoing_edge_ < no_outgoing_edges() && no_incoming_edges() > 0 && incoming_edge_ == no_primal_decision) {
      assert(incoming_edge_ == no_primal_decision);
      incoming_edge_ = incoming_cand;
      return;
    }
    // no edge labelled yet
    if(lb < 0.0) {
      incoming_edge_ = incoming_cand;
      outgoing_edge_ = outgoing_cand;
      return;
    } else {
      incoming_edge_ = no_edge_taken;
      outgoing_edge_ = no_edge_taken;
    }
  }
  REAL LowerBound() const {
    return std::min(cost_of_detection(), 0.0);
  }

  //REAL& detection_cost() { return *pot_; }

  REAL* incoming_begin() const { return pot_+1; };
  REAL* incoming_end() const { return pot_+1+no_incoming_edges_; };
  INDEX no_incoming_edges() const { return no_incoming_edges_; }

  REAL* outgoing_begin() const { return pot_+1+no_incoming_edges_; };
  REAL* outgoing_end() const { return pot_+size(); };
  INDEX no_outgoing_edges() const { return no_outgoing_edges_; }

  REAL EvaluatePrimal() const {
    assert(no_incoming_edges_ > 0);
    assert(no_outgoing_edges_ > 0);
    if(incoming_edge_ < no_incoming_edges() && outgoing_edge_ < no_outgoing_edges()) {
      return pot_[0] + incoming(incoming_edge_) + outgoing(outgoing_edge_);
    } else {
      return 0.0;
    }
  }

  INDEX size() const { return 1 + no_incoming_edges_ + no_outgoing_edges_; }

  void init_primal() 
  { 
    incoming_edge_ = no_primal_decision;
    outgoing_edge_ = no_primal_decision; 
  }
  template<class ARCHIVE> void serialize_primal(ARCHIVE& ar) { ar( incoming_edge_, outgoing_edge_ ); }
  template<class ARCHIVE> void serialize_dual(ARCHIVE& ar) { ar( cereal::binary_data( pot_, sizeof(REAL)*size()) ); }

  INDEX outgoing_edge() const { return outgoing_edge_; }
  INDEX incoming_edge() const { return incoming_edge_; }

  bool detection_active() const {
    return (incoming_edge_ < no_incoming_edges() || outgoing_edge_ < no_outgoing_edges());
  }

  template<typename SAT_SOLVER>
  void construct_sat_clauses(SAT_SOLVER& s) const
  {
    // create variables
    auto detection_var = lglmaxvar(s);
    lglincvar(s); 
    //auto detection_var = s.nVars();
    //s.new_var(); //detection will rather be false
    //std::cout << "first var in detection factor = " << detection_var << "\n";
    auto incoming_var = create_sat_variables(s, no_incoming_edges());
    auto outgoing_var = create_sat_variables(s, no_outgoing_edges());
    auto outgoing_sum = add_at_most_one_constraint_sat(s, outgoing_var.begin(), outgoing_var.end());
    auto incoming_sum = add_at_most_one_constraint_sat(s, incoming_var.begin(), incoming_var.end());

    // detection var must be equal to incoming and outgoing var
    make_sat_var_equal(s, to_literal(detection_var),to_literal(incoming_sum));
    make_sat_var_equal(s, to_literal(detection_var),to_literal(outgoing_sum));
  }

  template<typename VEC>
  void reduce_sat(VEC& assumptions, const REAL th, sat_var begin) const
  {
    const REAL cost = LowerBound();
    if(cost <= -th) {
      //std::cout << "in reduction: detection factor must be on\n";
      assumptions.push_back(to_literal(begin));
    }
    if(cost <= th) {
      const REAL incoming_min = *std::min_element(incoming_begin(), incoming_end());
      for(INDEX i=0; i<no_incoming_edges(); ++i) {
        if(incoming(i) > incoming_min + th) { 
           assumptions.push_back(-to_literal(begin+1+i));
         }
      }

      const REAL outgoing_min = *std::min_element(outgoing_begin(), outgoing_end());
      for(INDEX i=0; i<no_outgoing_edges(); ++i) {
        if(outgoing(i) > outgoing_min + th) { 
           assumptions.push_back(-to_literal(begin+1+no_incoming_edges()+i));
         }
      } 
    } else {
      for(auto i=0; i<this->size(); ++i) {
        assumptions.push_back(-to_literal(begin+i));
      }
    }
  }

  template<typename SAT_SOLVER>
  void convert_primal(SAT_SOLVER& s, sat_var first)
  {
    //assert(s.get_model()[first] != CMSat::l_Undef);
    //std::cout << lglmaxvar(s) << "\n";
    if(lglderef(s, to_literal(first)) == 1) {
    //if(s.get_model()[first] == CMSat::l_True) {
      //std::cout << "in primal conversion: detection factor is on\n";
      // find index of incoming and outgoing active edge
      for(INDEX i=first+1; i<first+1+no_incoming_edges(); ++i) {
        if(lglderef(s,to_literal(i)) == 1) {
        //if(s.get_model()[i] == CMSat::l_True) {
          incoming_edge_ = i-first-1;
          //std::cout << "incoming edge = " << i << ", ";
        }
      }
      for(INDEX i=first+1+no_incoming_edges(); i<first+size(); ++i) {
        if(lglderef(s,to_literal(i)) == 1) {
        //if(s.get_model()[i] == CMSat::l_True) {
          outgoing_edge_ = i-first-1-no_incoming_edges();
          //std::cout << "outgoing edge = " << i << "\n";
        }
      }
    } else {
      incoming_edge_ = no_edge_taken;
      outgoing_edge_ = no_edge_taken;
    }
  }


private:
  const INDEX no_incoming_edges_, no_outgoing_edges_;
  REAL* pot_;
  
  INDEX incoming_edge_, outgoing_edge_;
};

// message connecting outgoing edge to incoming edge of detection factors between consecutive timeframes
class transition_message {
public:
  transition_message(const bool split, const INDEX outgoing_edge_index, const INDEX incoming_edge_index) 
    : 
    outgoing_edge_index_(outgoing_edge_index),
    incoming_edge_index_(incoming_edge_index),
    split_(split)
  {}

  // send messages from detection factor along outgoing edges
  template<typename RIGHT_FACTOR, typename MSG_ARRAY, typename ITERATOR>
  static void SendMessagesToLeft(const RIGHT_FACTOR& rightFactor, MSG_ARRAY msg_begin, MSG_ARRAY msg_end, ITERATOR omega_begin)
  {
    REAL omega = 0.0;
    auto omega_it = omega_begin;
    for(auto it = msg_begin; it!=msg_end; ++it, ++omega_it) {
      omega += *omega_it;
    }
    assert(omega <= 1.0 + eps);
    assert(omega > 0.0);
    //std::cout << "send messages to left with omega = " << omega << "\n";

    // check #messages+1 = no incoming edges
    {
      INDEX c=0;
      for(auto it = msg_begin; it!=msg_end; ++it)  ++c;
      assert(c+1 == rightFactor.no_incoming_edges());
    }


    const REAL detection_outgoing_cost = rightFactor[0] + *std::min_element(rightFactor.outgoing_begin(), rightFactor.outgoing_end());

    std::vector<bool> edge_taken(rightFactor.no_incoming_edges(), false);
    omega_it = omega_begin;
    for(auto msg_it = msg_begin; msg_it!=msg_end; ++msg_it, ++omega_it) {
      if(*omega_it > 0.0) {
        edge_taken[(*msg_it).GetMessageOp().incoming_edge_index_] = true;
      }
    }

    REAL smallest_taken = std::numeric_limits<REAL>::infinity();
    REAL second_smallest_taken = std::numeric_limits<REAL>::infinity();
    REAL smallest_not_taken = std::numeric_limits<REAL>::infinity();
    for(INDEX i=0; i<rightFactor.no_incoming_edges(); ++i) {
      const REAL val = rightFactor.incoming(i);
      if(edge_taken[i]) {
        const REAL min = std::min(smallest_taken, val);
        const REAL max = std::max(smallest_taken, val);
        smallest_taken = min;
        second_smallest_taken = std::min(max, second_smallest_taken);
      } else {
        smallest_not_taken = std::min(val, smallest_not_taken); 
      }
    }

    REAL set_to_cost;
    if(smallest_not_taken < smallest_taken) {
      set_to_cost = std::min(detection_outgoing_cost + smallest_not_taken, 0.0);
    } else {
      set_to_cost = std::min({detection_outgoing_cost + second_smallest_taken, detection_outgoing_cost + smallest_not_taken, 0.0});
    }

    // compute smallest and second smallest value over all incoming edges
    assert(rightFactor.no_incoming_edges() > 1);
    //const auto smallest_incoming = two_smallest_elements<REAL>(rightFactor.incoming_begin(), rightFactor.incoming_end()-1); // do not take into account disappearance cost 
    //const auto smallest_incoming = two_smallest_elements<REAL>(rightFactor.incoming_begin(), rightFactor.incoming_end()); // do not take into account disappearance cost 

    //const REAL appearance_cost = rightFactor.incoming(rightFactor.no_incoming_edges()-1);
    //const REAL set_to_cost = std::min(detection_outgoing_cost + std::min(appearance_cost, smallest_incoming[1]), 0.0);
    //const REAL set_to_cost = std::min(detection_outgoing_cost + smallest_incoming[1], 0.0);
    omega_it = omega_begin;
    for(auto msg_it=msg_begin; msg_it!=msg_end; ++msg_it, ++omega_it) {
      if(*omega_it > 0.0) {
        const REAL msg = detection_outgoing_cost + rightFactor.incoming( (*msg_it).GetMessageOp().incoming_edge_index_ ) - set_to_cost;
        (*msg_it)[0] -= omega*msg;
      }
    } 
  }

  // send messages from detection factor along incoming edges
  template<typename LEFT_FACTOR, typename MSG_ARRAY, typename ITERATOR>
  static void SendMessagesToRight(const LEFT_FACTOR& leftFactor, MSG_ARRAY msg_begin, MSG_ARRAY msg_end, ITERATOR omega_begin)
  {
    REAL omega = 0.0;
    auto omega_it = omega_begin;
    for(auto it = msg_begin; it!=msg_end; ++it, ++omega_it) {
      omega += *omega_it;
    }
    assert(omega > 0.0);
    assert(omega <= 1.0 + eps);

    // do zrobienia: check whether no of messages is greater than no of outgoing edges. Multiple messages can act on same edge when divisions occur
    {
      INDEX c=0;
      for(auto it = msg_begin; it!=msg_end; ++it)  ++c;
      assert(c+1 >= std::distance(leftFactor.outgoing_begin(), leftFactor.outgoing_end()));
    }


    const REAL detection_incoming_cost = leftFactor[0] + *std::min_element(leftFactor.incoming_begin(), leftFactor.incoming_end());

    std::vector<bool> edge_taken(leftFactor.no_outgoing_edges(), false);
    omega_it = omega_begin;
    for(auto msg_it = msg_begin; msg_it!=msg_end; ++msg_it, ++omega_it) {
      if(*omega_it > 0.0) {
        edge_taken[(*msg_it).GetMessageOp().outgoing_edge_index_] = true;
      }
    }

    REAL smallest_taken = std::numeric_limits<REAL>::infinity();
    REAL second_smallest_taken = std::numeric_limits<REAL>::infinity();
    REAL smallest_not_taken = std::numeric_limits<REAL>::infinity();
    for(INDEX i=0; i<leftFactor.no_outgoing_edges(); ++i) {
      const REAL val = leftFactor.outgoing(i);
      if(edge_taken[i]) {
        const REAL min = std::min(smallest_taken, val);
        const REAL max = std::max(smallest_taken, val);
        smallest_taken = min;
        second_smallest_taken = std::min(max, second_smallest_taken);
      } else {
        smallest_not_taken = std::min(val, smallest_not_taken); 
      }
    }

    REAL set_to_cost;
    if(smallest_not_taken < smallest_taken) {
      set_to_cost = std::min(detection_incoming_cost + smallest_not_taken, 0.0);
    } else {
      set_to_cost = std::min({detection_incoming_cost + second_smallest_taken, detection_incoming_cost + smallest_not_taken, 0.0});
    }




    // compute smallest and second smallest value over all outgoing edges
    //const auto smallest_outgoing = two_smallest_elements<REAL>(leftFactor.outgoing_begin(), leftFactor.outgoing_end()-1); // omit disappearance term
    //const auto smallest_outgoing = two_smallest_elements<REAL>(leftFactor.outgoing_begin(), leftFactor.outgoing_end()); // omit disappearance term

    //const REAL disappearance_cost = leftFactor.outgoing(leftFactor.no_outgoing_edges()-1);
    //const REAL set_to_cost = std::min(detection_incoming_cost + std::min(disappearance_cost, smallest_outgoing[1]), 0.0);
    //const REAL set_to_cost = std::min(detection_incoming_cost + smallest_outgoing[1], 0.0);
    omega_it = omega_begin;
    for(; msg_begin!=msg_end; ++msg_begin, ++omega_it) {
      if(*omega_it > 0.0) {
        const REAL w = (*msg_begin).GetMessageOp().split_ ? 0.5 : 1.0;
        const REAL msg = w*(detection_incoming_cost + leftFactor.outgoing( (*msg_begin).GetMessageOp().outgoing_edge_index_ ) - set_to_cost);
        (*msg_begin)[0] -= omega*msg;
      }
    } 
    return;
  }

  template<typename LEFT_FACTOR, typename RIGHT_FACTOR>
  bool ComputeRightFromLeftPrimal(const LEFT_FACTOR& l, RIGHT_FACTOR& r)
  {
    if(l.outgoing_edge_ == outgoing_edge_index_ && r.incoming_edge_ != incoming_edge_index_) {
      r.incoming_edge_ = incoming_edge_index_;
      return true;
    } else {
      return false;
    }
  }
  template<typename LEFT_FACTOR, typename RIGHT_FACTOR>
  bool ComputeLeftFromRightPrimal(LEFT_FACTOR& l, const RIGHT_FACTOR& r) 
  {
    if(r.incoming_edge_ == incoming_edge_index_ && r.outgoing_edge_ != outgoing_edge_index_) {
      l.outgoing_edge_ = outgoing_edge_index_;
      return true;
    } else {
      return false;
    }
  }

  template<typename RIGHT_FACTOR, typename G2>
  void send_message_to_left(RIGHT_FACTOR& r, G2& msg, const REAL omega)
  {
    // to do: measure time! is this or the construction below faster?
    const REAL detection_outgoing_cost = r[0] + *std::min_element(r.outgoing_begin(), r.outgoing_end());

    const REAL incoming_val = r.incoming(incoming_edge_index_);
    r.incoming(incoming_edge_index_) = std::numeric_limits<REAL>::infinity();
    const REAL min_incoming_val = *std::min_element(r.incoming_begin(), r.incoming_end());
    r.incoming(incoming_edge_index_) = incoming_val;

    msg[0] -= omega*( detection_outgoing_cost + incoming_val - std::min(detection_outgoing_cost + min_incoming_val, 0.0) ); 
  }
  
  template<typename RIGHT_FACTOR, typename G2>
  void ReceiveRestrictedMessageFromRight(RIGHT_FACTOR& r, G2& msg)
  { 
    if(r.incoming_edge_ < r.no_incoming_edges()) { // incoming edge has already been picked
      if(r.incoming_edge_ == incoming_edge_index_) {
        const REAL val_prev = msg.GetLeftFactor()->GetFactor()->outgoing(outgoing_edge_index_);
        msg[0] -= -std::numeric_limits<REAL>::infinity(); 
        const REAL val = msg.GetLeftFactor()->GetFactor()->outgoing(outgoing_edge_index_);
        //assert(msg.GetLeftFactor()->GetFactor()->outgoing(outgoing_edge_index_) < -10000);
      } else {
        const REAL val_prev = msg.GetLeftFactor()->GetFactor()->outgoing(outgoing_edge_index_);
        msg[0] -= std::numeric_limits<REAL>::infinity(); 
        const REAL val = msg.GetLeftFactor()->GetFactor()->outgoing(outgoing_edge_index_);
        //assert(msg.GetLeftFactor()->GetFactor()->outgoing(outgoing_edge_index_) > 10000);
      }
    } else { // no incoming edge has been picked yet.
      return;
      //send_message_to_left(r,msg); 
    }
  }

  template<typename LEFT_FACTOR, typename G2>
  void send_message_to_right(LEFT_FACTOR& l, G2& msg, const REAL omega)
  { 
    const REAL detection_incoming_cost = l[0] + *std::min_element(l.incoming_begin(), l.incoming_end());

    const REAL outgoing_val = l.outgoing(outgoing_edge_index_);
    l.outgoing(outgoing_edge_index_) = std::numeric_limits<REAL>::infinity();
    const REAL min_outgoing_val = *std::min_element(l.outgoing_begin(), l.outgoing_end());
    l.outgoing(outgoing_edge_index_) = outgoing_val;

    msg[0] -= omega*( detection_incoming_cost + outgoing_val - std::min(detection_incoming_cost + min_outgoing_val, 0.0) );
  }

  template<typename LEFT_FACTOR, typename G2>
  void ReceiveRestrictedMessageFromLeft(LEFT_FACTOR& l, G2& msg)
  { 
    if(l.outgoing_edge_ < l.no_outgoing_edges()) { // outgoing edge has already been picked
      if(l.outgoing_edge_ == outgoing_edge_index_) {
        msg[0] -= -std::numeric_limits<REAL>::infinity(); 
        assert(msg.GetRightFactor()->GetFactor()->incoming(incoming_edge_index_) < -10000);
      } else {
        msg[0] -= std::numeric_limits<REAL>::infinity(); 
        assert(msg.GetRightFactor()->GetFactor()->incoming(incoming_edge_index_) > 10000);
      }
    } else { // no incoming edge has been picked yet.
      return;
      //send_message_to_right(l,msg); 
    }
  }

  template<typename G>
  void RepamLeft(G& l, const REAL msg, const INDEX msg_dim)
  {
    assert(msg_dim == 0);
    l.outgoing(outgoing_edge_index_) += msg;
  }
  template<typename G>
  void RepamRight(G& r, const REAL msg, const INDEX msg_dim)
  {
    assert(msg_dim == 0);
    r.incoming(incoming_edge_index_) += msg;
  }

  template<typename LEFT_FACTOR, typename RIGHT_FACTOR>
  bool CheckPrimalConsistency(const LEFT_FACTOR& l, const RIGHT_FACTOR& r) const
  {
    if(l.outgoing_edge_ == outgoing_edge_index_) {
      return r.incoming_edge_ == incoming_edge_index_;
    }
    if(r.incoming_edge_ == incoming_edge_index_) {
      return r.outgoing_edge_ == outgoing_edge_index_;
    } 
    return true;
  }

  template<typename SAT_SOLVER, typename LEFT_FACTOR, typename RIGHT_FACTOR>
  void construct_sat_clauses(SAT_SOLVER& s, LEFT_FACTOR& l, RIGHT_FACTOR& r, sat_var left_begin, sat_var right_begin) const
  {
    auto left_var = left_begin+1+l.no_incoming_edges() + outgoing_edge_index_;
    auto right_var = right_begin+1+ incoming_edge_index_;
    make_sat_var_equal(s, to_literal(left_var), to_literal(right_var));
  }
private:
  const INDEX outgoing_edge_index_;
  const INDEX incoming_edge_index_;
  const bool split_; // is split really needed? Whether it is a split transition can be found out in SendMessagesTo... by checking whether there is two outgoing_edge_indices
};

// multiple cell detection hypotheses can be mutually exclusive. 
// simplex x1 + ... + xn <= 1
// to account for overlapping detections: only one can be active
class at_most_one_cell_factor : public vector<REAL> {

  friend class at_most_one_cell_message;

public:

  constexpr static INDEX no_primal_decision = std::numeric_limits<INDEX>::max();
  constexpr static INDEX no_primal_active = std::numeric_limits<INDEX>::max()-1;
  constexpr static INDEX primal_infeasible = std::numeric_limits<INDEX>::max()-2;

   at_most_one_cell_factor(const INDEX size) : vector<REAL>(size)
   {}

   REAL LowerBound() const {
     return std::min(0.0, *std::min_element(this->begin(), this->end()));
   }

   REAL EvaluatePrimal() const {
     //if(primal_ == no_primal_decision) { primal_ = no_primal_active; }

     if(primal_ < size()) {
       return (*this)[primal_];
     } else if(primal_ == no_primal_active || no_primal_decision) {
       return 0.0;
     } else {
       return std::numeric_limits<REAL>::infinity();
     }
   }

   void init_primal() { primal_ = no_primal_decision; }
   template<class ARCHIVE> void serialize_primal(ARCHIVE& ar) { ar(primal_); }
   template<class ARCHIVE> void serialize_dual(ARCHIVE& ar) { ar( *static_cast<vector<REAL>*>(this) ); }

  template<typename SAT_SOLVER>
  void construct_sat_clauses(SAT_SOLVER& s) const
  {
    auto var = create_sat_variables(s, this->size());
    add_at_most_one_constraint_sat(s,var.begin(), var.end());
  }

  template<typename VEC>
  void reduce_sat(VEC& assumptions, const REAL th, sat_var begin) const
  {
    for(INDEX i=0; i<this->size(); ++i) {
      if((*this)[i] > th) { 
        assumptions.push_back(-to_literal(begin+i)); 
      }
    }
  }

  template<typename SAT_SOLVER>
  void convert_primal(SAT_SOLVER& s, sat_var first)
  {
    for(INDEX i=0; i<this->size(); ++i) {
      if(lglderef(s, to_literal(first+i)) == true) {
      //if(s.get_model()[first+i] == CMSat::l_True) {
        primal_ = i;
        return;
      }
    }
    primal_ = no_primal_active;
  }

private:
   INDEX primal_;
};

// connecting first entry of detection factor with one entry in at_most_one_cell_factor
// left factor is detection factor, right one at_most_one_cell_factor
class at_most_one_cell_message {
public:
  at_most_one_cell_message(const INDEX at_most_one_cell_factor_index) : at_most_one_cell_factor_index_(at_most_one_cell_factor_index) {}

  template<typename G>
  void RepamLeft(G& l, const REAL msg, const INDEX msg_dim)
  {
    assert(msg_dim == 0);
    l[0] += msg;
    assert(!std::isnan(l[0]));
  }

  template<typename G>
  void RepamRight(G& r, const REAL msg, const INDEX msg_dim)
  {
    assert(msg_dim == 0);
    r[at_most_one_cell_factor_index_] += msg;
  }

  template<typename LEFT_FACTOR, typename RIGHT_FACTOR>
  void ComputeRightFromLeftPrimal(const LEFT_FACTOR& l, RIGHT_FACTOR& r)
  {
    if(l.incoming_edge_ < l.no_incoming_edges() || l.outgoing_edge_ < l.no_outgoing_edges()) {
      if(r.primal_ == at_most_one_cell_factor::no_primal_decision || r.primal_ == at_most_one_cell_factor_index_) {
        r.primal_ = at_most_one_cell_factor_index_;
      } else {
        r.primal_ = at_most_one_cell_factor::primal_infeasible; // more than two detections active
      }
    } else {
      if(r.primal_ == r.no_primal_decision) {
        r.primal_ = r.no_primal_active;
      }
    }
  }

  template<typename LEFT_FACTOR, typename MSG_ARRAY, typename ITERATOR>
  static void SendMessagesToRight(const LEFT_FACTOR& l, MSG_ARRAY msg_begin, MSG_ARRAY msg_end, ITERATOR omegaIt)
  {
    const REAL delta = l.cost_of_detection();
    assert(!std::isnan(delta));

    for(; msg_begin!=msg_end; ++msg_begin, ++omegaIt) {
      (*msg_begin)[0] -= (*omegaIt)*delta;
    } 
  }


  template<typename RIGHT_FACTOR, typename MSG>
  void ReceiveRestrictedMessageFromRight(RIGHT_FACTOR& r, MSG& msg) 
  {
    assert(at_most_one_cell_factor_index_ < r.size());
    //std::cout << "r.primal = " << r.primal_ << " = " << r.no_primal_active << "\n";
    if(r.primal_ == r.no_primal_active || r.primal_ == r.no_primal_decision) { // no element chosen yet
      //make_right_factor_uniform(r, msg);
      const REAL cur_detection_cost = r[at_most_one_cell_factor_index_];
      r[at_most_one_cell_factor_index_] = std::numeric_limits<REAL>::infinity();
      const REAL rest_cost = std::min(0.0, *std::min_element(r.begin(), r.end()));
      r[at_most_one_cell_factor_index_] = cur_detection_cost;

      msg[0] -= (cur_detection_cost - rest_cost);
    } else if(r.primal_ < r.size()) { // one element already chosen
      if(r.primal_ == at_most_one_cell_factor_index_) {
        const REAL val_prev = msg.GetLeftFactor()->GetFactor()->operator[](0);
        msg[0] -= -100000;//std::numeric_limits<REAL>::infinity();
        const REAL val = msg.GetLeftFactor()->GetFactor()->operator[](0);
        //assert(msg.GetLeftFactor()->GetFactor()->operator[](0) < -10000);
      } else {
        //assert(r.primal_ < r.size());
        const REAL test_val_prev = msg.GetLeftFactor()->GetFactor()->operator[](0);
        msg[0] -= 100000;//std::numeric_limits<REAL>::infinity();
        const REAL test_val = msg.GetLeftFactor()->GetFactor()->operator[](0);
        //assert(test_val > 10000);
        //assert(msg.GetLeftFactor()->GetFactor()->operator[](0) > 10000);
      }
    } else {
      assert(r.primal_ == r.primal_infeasible);
      const REAL test_val_prev = msg.GetLeftFactor()->GetFactor()->operator[](0);
      msg[0] -= 100000;//std::numeric_limits<REAL>::infinity();
      const REAL test_val = msg.GetLeftFactor()->GetFactor()->operator[](0);
      //assert(test_val > 10000);
      //assert(msg.GetLeftFactor()->GetFactor()->operator[](0) > 10000); 
    }
  }
  /*
  template<typename RIGHT_FACTOR, typename MSG_ARRAY, typename ITERATOR>
  static void SendMessagesToLeft(const RIGHT_FACTOR& rightFactor, MSG_ARRAY msg_begin, MSG_ARRAY msg_end, ITERATOR omegaIt)
  {
    REAL omega = 0.0;
    for(auto it = msg_begin; it!=msg_end; ++it, ++omegaIt) {
      omega += *omegaIt;
    }
    assert(omega <= 1.0 + eps);

    INDEX c=0;
    for(auto it = msg_begin; it!=msg_end; ++it) ++c;
    assert(rightFactor.size() == c);

    for(; msg_begin!=msg_end; ++msg_begin) {
      (*msg_begin)[0] -= omega*rightFactor[ (*msg_begin).GetMessageOp().at_most_one_cell_factor_index_ ];
    }
  }
  */

  template<typename LEFT_FACTOR, typename G2>
  void send_message_to_right(const LEFT_FACTOR& l, G2& msg, const REAL omega = 1.0)
  {
    // make cost of detection same as cost of non-detection
    msg[0] -= omega*l.cost_of_detection();
  }
  template<typename RIGHT_FACTOR, typename G2>
  void send_message_to_left(RIGHT_FACTOR& r, G2& msg, const REAL omega = 1.0)
  {
    const REAL cur_detection_cost = r[at_most_one_cell_factor_index_];
    r[at_most_one_cell_factor_index_] = std::numeric_limits<REAL>::infinity();
    const REAL rest_cost = std::min(0.0, *std::min_element(r.begin(), r.end()));
    r[at_most_one_cell_factor_index_] = cur_detection_cost;

    msg[0] -= omega*(cur_detection_cost - rest_cost);
  }

  template<typename LEFT_FACTOR, typename RIGHT_FACTOR>
  bool CheckPrimalConsistency(const LEFT_FACTOR& l, const RIGHT_FACTOR& r) const
  {
    if(r.primal_ == r.primal_infeasible) {
      std::cout << "exclusion constraints not satisfied!\n";
      return false;
    }
    if(l.detection_active()) {
      return (r.primal_ == at_most_one_cell_factor_index_);
    } else {
      return (r.primal_ != at_most_one_cell_factor_index_);
    }
  }

  template<typename SAT_SOLVER, typename LEFT_FACTOR, typename RIGHT_FACTOR>
  void construct_sat_clauses(SAT_SOLVER& s, LEFT_FACTOR& l, RIGHT_FACTOR& r, sat_var left_begin, sat_var right_begin) const
  {
    auto right_var = right_begin+at_most_one_cell_factor_index_;
    make_sat_var_equal(s, to_literal(left_begin), to_literal(right_var));
  }
private:
   const INDEX at_most_one_cell_factor_index_;
};

// first entry: lower exit transition. second entry: upper transition messages (excluding exit transition)
// x_1 = 1 implies x_2 = 0
// equivalent to x_1 + x_2 <= 1
class exit_constraint_factor : public std::array<REAL,2> {
  template<exit_constraint_position POSITION> friend class exit_constraint_message;
public:

  constexpr static INDEX no_primal_decision = std::numeric_limits<INDEX>::max();
  constexpr static INDEX no_primal_active = std::numeric_limits<INDEX>::max()-1;
  constexpr static INDEX primal_infeasible = std::numeric_limits<INDEX>::max()-2; 

  using repam_type = std::array<REAL,2>;
  REAL LowerBound() const 
  {
    return std::min(0.0, std::min( (*this)[0], (*this)[1] )); 
  }
  REAL EvaluatePrimal() const 
  {
    if(primal_ < this->size()) {
      return (*this)[primal_];
    } else if(primal_ == no_primal_active) {
      return 0.0;
    } else {
      return std::numeric_limits<REAL>::infinity();
    }
  }

  void init_primal() { primal_ = no_primal_decision; }
  template<class ARCHIVE> void serialize_primal(ARCHIVE& ar) { ar( primal_ ); }
  template<class ARCHIVE> void serialize_dual(ARCHIVE& ar) { ar( *static_cast<repam_type*>(this) ); }

  template<typename SAT_SOLVER>
  void construct_sat_clauses(SAT_SOLVER& s) const
  {
    auto var = create_sat_variables(s, this->size());
    add_at_most_one_constraint_sat(s,var.begin(), var.end());
  }

  template<typename VEC>
  void reduce_sat(VEC& assumptions, const REAL th, sat_var begin) const
  {
    // do zrobienia: if reparametrization is < -th, then disallow no edge to be taken!
    for(INDEX i=0; i<this->size(); ++i) {
      if((*this)[i] > th) { 
        assumptions.push_back(-to_literal(begin+i)); 
      }
    }
  }

  template<typename SAT_SOLVER>
  void convert_primal(SAT_SOLVER& s, sat_var first) 
  {
    for(INDEX i=0; i<this->size(); ++i) {
      if(lglderef(s, to_literal(first+i)) == true) {
      //if(s.get_model()[first+i] == CMSat::l_True) {
        primal_ = i;
        return;
      }
    }
    primal_ = no_primal_active;
  }

private:
  INDEX primal_;
};

// left is detection_factor, right is exit_constraint_factor
template<exit_constraint_position POSITION>
class exit_constraint_message {
public:
  template<typename G>
  void RepamLeft(G& l, const REAL msg, const INDEX msg_dim)
  {
    assert(msg_dim == 0);
    if(POSITION == exit_constraint_position::lower) {
      const INDEX i = l.no_outgoing_edges()-1;
      l.outgoing(i) += msg;
    } else {
      for(INDEX i=0; i<l.no_outgoing_edges()-1; ++i) {
        l.outgoing(i) += msg;
      }
    }
  }

  template<typename G>
  void RepamRight(G& r, const REAL msg, const INDEX msg_dim)
  {
    assert(msg_dim == 0);
    if(POSITION == exit_constraint_position::lower) {
      r[0] += msg;
    } else {
      r[1] += msg;
    }
  }

  template<typename LEFT_FACTOR, typename RIGHT_FACTOR>
  void ComputeRightFromLeftPrimal(const LEFT_FACTOR& l, RIGHT_FACTOR& r)
  {
    //assert(false);
    if(POSITION == exit_constraint_position::lower) {
      if(l.outgoing_edge_ == l.no_outgoing_edges() - 1) {
        r.primal_ = 0;
      } 
    } else {

    }
  }

  template<typename LEFT_FACTOR, typename MSG>
  void send_message_to_right(const LEFT_FACTOR& l, MSG& msg, const REAL omega)
  { 
    REAL detection_incoming_cost = l[0];
    detection_incoming_cost += *std::min_element(l.incoming_begin(), l.incoming_end());
    REAL outgoing_cost = std::numeric_limits<REAL>::infinity();
    assert(l.no_outgoing_edges() > 0);
    if(l.no_outgoing_edges() > 1) {
      outgoing_cost = *std::min_element(l.outgoing_begin(), l.outgoing_end()-1);
    } 

    if(POSITION == exit_constraint_position::lower) {
      msg[0] -= omega*( l.outgoing( l.no_outgoing_edges()-1 ) - std::min( -detection_incoming_cost, outgoing_cost ) );
    } else {
      msg[0] -= omega*( outgoing_cost - std::min( -detection_incoming_cost, l.outgoing( l.no_outgoing_edges()-1 ) ) );
    }
  }


  template<typename LEFT_FACTOR, typename MSG_ARRAY, typename ITERATOR>
  static void SendMessagesToRight(const LEFT_FACTOR& l, MSG_ARRAY msg_begin, MSG_ARRAY msg_end, ITERATOR omegaIt)
  {
    REAL detection_incoming_cost = l[0];
    detection_incoming_cost += *std::min_element(l.incoming_begin(), l.incoming_end());
    REAL outgoing_cost = std::numeric_limits<REAL>::infinity();
    assert(l.no_outgoing_edges() > 0);
    if(l.no_outgoing_edges() > 1) {
      outgoing_cost = *std::min_element(l.outgoing_begin(), l.outgoing_end()-1);
    } 

    if(POSITION == exit_constraint_position::lower) {
      const REAL delta = l.outgoing( l.no_outgoing_edges()-1 ) - std::min( -detection_incoming_cost, outgoing_cost );
      assert(!std::isnan(delta));
      for(; msg_begin!=msg_end; ++msg_begin, ++omegaIt) { (*msg_begin)[0] -=  (*omegaIt)*delta; }
    } else {
      const REAL delta = outgoing_cost - std::min( -detection_incoming_cost, l.outgoing( l.no_outgoing_edges()-1 ) );
      assert(!std::isnan(delta));
      for(; msg_begin!=msg_end; ++msg_begin, ++omegaIt) { (*msg_begin)[0] -=  (*omegaIt)*delta; }
    }
  }

  template<typename RIGHT_FACTOR, typename MSG>
  void send_message_to_left(const RIGHT_FACTOR& r, MSG& msg, const REAL omega) 
  {
    assert(r.size() == 2);
    if(POSITION == exit_constraint_position::lower) {
      msg[0] -= omega*(r[0] - std::min(0.0, r[1]));
    } else {
      msg[0] -= omega*(r[1] - std::min(0.0, r[0]));
    }
  }

  template<typename LEFT_FACTOR, typename RIGHT_FACTOR>
  bool CheckPrimalConsistency(const LEFT_FACTOR& l, const RIGHT_FACTOR& r) const
  {
    if(POSITION == exit_constraint_position::lower) {
      if(r.primal_ == 0) {
        return l.outgoing_edge_ == l.no_outgoing_edges()-1;
      } else {
        return l.outgoing_edge_ < l.no_outgoing_edges()-1 || l.outgoing_edge_ == l.no_edge_taken;
      }
    } else {
      if(r.primal_ == 1) {
        return l.outgoing_edge_ < l.no_outgoing_edges()-1;
      } else {
        return l.outgoing_edge_ == l.no_outgoing_edges()-1 || l.outgoing_edge_ == l.no_edge_taken;
      }
    }
  }

  template<typename SAT_SOLVER, typename LEFT_FACTOR, typename RIGHT_FACTOR>
  void construct_sat_clauses(SAT_SOLVER& s, LEFT_FACTOR& l, RIGHT_FACTOR& r, sat_var left_begin, sat_var right_begin) const
  {
    if(POSITION == exit_constraint_position::lower) {
      make_sat_var_equal(s, to_literal(left_begin+l.size()-1), to_literal(right_begin));
    } else {
      // get indicator variable for sum of outgoing edges without last one
      std::vector<sat_var> outgoing_var(l.no_outgoing_edges()-1);
      for(INDEX i=0; i<outgoing_var.size(); ++i) {
        outgoing_var[i] = left_begin+1+l.no_incoming_edges()+i;
      }
      auto c = add_at_most_one_constraint_sat(s, outgoing_var.begin(), outgoing_var.end());
      make_sat_var_equal(s, to_literal(c), to_literal(right_begin+1));
    }
  }
};



// for conservation tracking

// detection can contain multiple cells
// layout of reparametrization: detection costs, appearance costs, disappearance costs, incoming edge costs, outgoing edge costs
class multiple_detection_factor {
public:
   multiple_detection_factor(const INDEX max_detections, const INDEX no_incoming_edges, const INDEX no_outgoing_edges, bool divides = false)
      : max_detections_(max_detections),
        no_incoming_edges_(no_incoming_edges),
        no_outgoing_edges_(no_outgoing_edges), // the last outgoing edge is the exit one (when exit is possible)
        divides_(divides)
  {
    pot_ = global_real_block_allocator.allocate(size());
    assert(pot_ != nullptr);
    std::fill(pot_, pot_ + size(), 0.0);
  }
  ~multiple_detection_factor() {
    global_real_block_allocator.deallocate(pot_,1);
  }
  multiple_detection_factor(const multiple_detection_factor& o) 
     : max_detections_(o.max_detections_),
     no_incoming_edges_(o.no_incoming_edges_),
     no_outgoing_edges_(o.no_outgoing_edges_)
  {
    pot_ = global_real_block_allocator.allocate(size());
    assert(pot_ != nullptr);
    for(INDEX i=0; i<size(); ++i) { pot_[i] = o.pot_[i]; }
  }
  void operator=(const multiple_detection_factor& o) {
     assert(max_detections_ == o.max_detections_);
     assert(no_incoming_edges_ == o.no_incoming_edges_);
     assert(no_outgoing_edges_ == o.no_outgoing_edges_);
     for(INDEX i=0; i<size(); ++i) { pot_[i] = o.pot_[i]; }
  }

  REAL operator[](const INDEX i) const {
    assert(i<size());
    return pot_[i];
  }
  REAL& operator[](const INDEX i) {
    assert(i<size());
    return pot_[i];
  }

  REAL detection(const INDEX i) const { return pot_[i]; }
  REAL& detection(const INDEX i) { return pot_[i]; }
  REAL appearance(const INDEX i) const { return pot_[max_detections_ + i]; }
  REAL& appearance(const INDEX i) { return pot_[max_detections_ + i]; }
  REAL disappearance(const INDEX i) const { return pot_[2*max_detections_ + i]; }
  REAL& disappearance(const INDEX i) { return pot_[2*max_detections_ + i]; }
    
  REAL incoming(const INDEX i) const {
    assert(i<no_incoming_edges_);
    return pot_[3*max_detections_+i];
  }
  REAL& incoming(const INDEX i) {
    assert(i<no_incoming_edges_);
    return pot_[3*max_detections_+i];
  }
  REAL outgoing(const INDEX i) const {
    assert(i<no_outgoing_edges_);
    return pot_[3*max_detections_+no_incoming_edges_+i];
  }
  REAL& outgoing(const INDEX i) {
    assert(i<no_outgoing_edges_);
    return pot_[3*max_detections_+no_incoming_edges_+i];
  }

  REAL* detection_begin() const { return pot_; }
  REAL* detection_end() const { return pot_ + max_detections_; }

  REAL* appearance_begin() const { return pot_ + max_detections_; }
  REAL* appearance_end() const { return pot_ + 2*max_detections_; }
  
  REAL* disappearance_begin() const { return pot_ + 2*max_detections_; }
  REAL* disappearance_end() const { return pot_ + 3*max_detections_; }

  REAL* incoming_begin() const { return pot_+3*max_detections_; };
  REAL* incoming_end() const { return pot_+3*max_detections_+no_incoming_edges_; };
  INDEX no_incoming_edges() const { return no_incoming_edges_; }

  REAL* outgoing_begin() const { return pot_+3*max_detections_+no_incoming_edges_; };
  REAL* outgoing_end() const { return pot_+size(); };
  INDEX no_outgoing_edges() const { return no_outgoing_edges_; }

  INDEX size() const { return 3*max_detections_ + no_incoming_edges_ + no_outgoing_edges_; }

  void init_primal() { incoming_edge_ = no_incoming_edges_; outgoing_edge_ = no_outgoing_edges_; }
  template<class ARCHIVE> void serialize_primal(ARCHIVE& ar) { ar( incoming_edge_, outgoing_edge_ ); }
  template<class ARCHIVE> void serialize_dual(ARCHIVE& ar) { ar( cereal::binary_data( pot_, sizeof(REAL)*size()) ); }

  void sort_incoming_cost(vector<REAL>& v) const
  {
     assert(v.size() == std::min(no_incoming_edges(), max_detections_));
     std::partial_sort_copy(incoming_begin(), incoming_end(), v.begin(), v.end());
     std::partial_sum(v.begin(), v.end(), v.begin());
  }

  void sort_outgoing_cost(vector<REAL>& v) const
  {
     assert(v.size() == std::min(no_outgoing_edges(), max_detections_));
     std::partial_sort_copy(outgoing_begin(), outgoing_end(), v.begin(), v.end());
     std::partial_sum(v.begin(), v.end(), v.begin()); 
  }

  REAL detection_cost(const vector<REAL>& smallest_incoming, const vector<REAL>& smallest_outgoing) const 
  {
     REAL min_detection_cost = 0.0; // for no detection

     assert(max_detections_ <= no_incoming_edges() || no_incoming_edges() == 0);
     assert(max_detections_ <= no_outgoing_edges() || no_outgoing_edges() == 0);

     for(INDEX no_detections=1; no_detections<=std::min(smallest_incoming.size(), smallest_outgoing.size()); ++no_detections) {
        REAL detection_cost = pot_[no_detections-1] + smallest_incoming[no_detections-1] + smallest_outgoing[no_detections-1];
        min_detection_cost = std::min(detection_cost, min_detection_cost); 
     }

     return min_detection_cost; 
  }

  REAL appearance_cost(const vector<REAL>& smallest_outgoing) const 
  {
    REAL min_appearance_cost = std::numeric_limits<REAL>::infinity();
    for(INDEX no_detections=1; no_detections<=max_detections_; ++no_detections) {
      REAL appearance_cost = appearance(no_detections-1) + outgoing(no_detections-1);
      min_appearance_cost = std::min(appearance_cost, min_appearance_cost);
    }
    return min_appearance_cost; 
  }

  REAL disappearance_cost(const vector<REAL>& smallest_incoming) const 
  {
    REAL min_disappearance_cost = std::numeric_limits<REAL>::infinity();
    for(INDEX no_detections=1; no_detections<=max_detections_; ++no_detections) {
      REAL disappearance_cost = disappearance(no_detections-1) + outgoing(no_detections-1);
      min_disappearance_cost = std::min(disappearance_cost, min_disappearance_cost);
    }
    return min_disappearance_cost;
  } 

  REAL LowerBound() const
  {
     vector<REAL> smallest_incoming( std::min(no_incoming_edges(), max_detections_) );
     sort_incoming_cost(smallest_incoming);
     vector<REAL> smallest_outgoing( std::min(no_outgoing_edges(), max_detections_) );
     sort_outgoing_cost(smallest_outgoing);

     const REAL min_detection_cost = detection_cost(smallest_incoming, smallest_outgoing);
     const REAL min_appearance_cost = appearance_cost(smallest_outgoing);
     const REAL min_disappearance_cost = disappearance_cost(smallest_incoming);

     return std::min({min_detection_cost, min_appearance_cost, min_disappearance_cost});
  }

  // conditioned on outgoing edge = 0/1, compute cost
  std::array<REAL,2> outgoing_cost(const INDEX outgoing_edge_index) const
  {
     const REAL current_arc_cost = outgoing(outgoing_edge_index);

     vector<REAL> smallest_incoming( std::min(no_incoming_edges(), max_detections_) );
     sort_incoming_cost(smallest_incoming);

     vector<REAL> smallest_outgoing( std::min(no_outgoing_edges(), max_detections_+1) );
     sort_outgoing_cost(smallest_outgoing);

     assert(max_detections_ <= no_incoming_edges() || no_incoming_edges() == 0);
     assert(max_detections_ <= no_outgoing_edges() || no_outgoing_edges() == 0);

     auto min_detection_cost = cost01(detection_begin(), smallest_incoming, smallest_outgoing, current_arc_cost);
     min_detection_cost[0] = std::min(min_detection_cost[0], 0.0);

     // disappearance cost
     const REAL min_disappearance_cost = disappearance_cost(smallest_incoming); // here outgoing edge is 0 by default

     // appearance cost
     auto min_appearance_cost = cost01(appearance_begin(), smallest_outgoing, current_arc_cost);

     const REAL cost_0 = std::min({min_detection_cost[0], min_disappearance_cost, min_appearance_cost[0]});
     const REAL cost_1 = std::min(min_detection_cost[1], min_appearance_cost[1]);
     return {cost_0, cost_1};

  }

  // compute cost when smallest_ is held to 0/1 for entry with cur_cost
  template<typename VEC1, typename VEC2, typename VEC3>
  std::array<REAL,2> cost01(const VEC1& cost1, const VEC2& cost2, const VEC3& smallest_, const REAL cur_cost) const
  {
     // current arc = 0
     REAL min_cost_0 = cost1[0] + cost2[0] + cur_cost==smallest_[0] ? smallest_[1]-cur_cost : smallest_[0];
     for(INDEX no_detections=2; no_detections<smallest_.size(); ++no_detections) { // do zrobienia: possibly this is false
       REAL cost_0 = cost1[no_detections-1] + cost2[no_detections-1];
       if(cur_cost <= smallest_[no_detections-1] - smallest_[no_detections-2]) {
         cost_0 += smallest_[no_detections] - cur_cost;
       } else {
         cost_0 += smallest_[no_detections-1];
       }
       min_cost_0 = std::min(min_cost_0, cost_0); 
     }

     // current arc = 1
     REAL min_cost_1 = cost1[0] + cost2[0] + cur_cost;
     for(INDEX no_detections=2; no_detections<=smallest_.size(); ++no_detections) {
        REAL cost_1 = cost1[no_detections-1] + cost2[no_detections-1];
        if(cur_cost > smallest_[no_detections-1] - smallest_[no_detections-2]) {
          cost_1 += cur_cost + smallest_[no_detections-2];
        } else {
          cost_1 += smallest_[no_detections-1];
        }
        min_cost_1 = std::min(min_cost_1, cost_1); 
     }

     return {min_cost_0, min_cost_1}; 
  }

  template<typename VEC1, typename VEC2>
  std::array<REAL,2> cost01(const VEC1& cost, const VEC2& smallest_, const REAL cur_cost) const
  {
    REAL cost_0 = cost[0] + smallest_[0] == cur_cost ? smallest_[1] - cur_cost : smallest_[0];
    for(INDEX no_detections=2; no_detections<smallest_.size(); ++no_detections) { // do zrobienia: possibly this is false
       REAL cost_0 = cost[no_detections-1] + smallest_[no_detections-1];
       if(cur_cost <= smallest_[no_detections-1] - smallest_[no_detections-2]) {
         cost_0 += smallest_[no_detections] - cur_cost;
       } else {
         cost_0 += smallest_[no_detections-1];
       }
    }

    REAL cost_1 = cost[0] + cur_cost;
    for(INDEX no_detections=2; no_detections<=smallest_.size(); ++no_detections) {
      REAL cost_1 = (*this)[no_detections-1];
      if(cur_cost > smallest_[no_detections-1] - smallest_[no_detections-2]) {
        cost_1 += cur_cost + smallest_[no_detections-2];
      } else {
        cost_1 += smallest_[no_detections-1];
      }
    }
    return {cost_0, cost_1};
  }

  std::array<REAL,2> incoming_cost(const INDEX incoming_edge_index) const
  {
     const REAL current_arc_cost = incoming(incoming_edge_index);

     vector<REAL> smallest_incoming( std::min(no_incoming_edges(), max_detections_) );
     sort_incoming_cost(smallest_incoming);

     vector<REAL> smallest_outgoing( std::min(no_outgoing_edges(), max_detections_+1) );
     sort_outgoing_cost(smallest_outgoing);

     assert(max_detections_ <= no_incoming_edges() || no_incoming_edges() == 0);
     assert(max_detections_ <= no_outgoing_edges() || no_outgoing_edges() == 0);

     auto min_detection_cost = cost01(detection_begin(), smallest_outgoing, smallest_incoming, current_arc_cost);
     min_detection_cost[0] = std::min(min_detection_cost[0], 0.0);

     // appearance cost
     const REAL min_appearance_cost = appearance_cost(smallest_outgoing); // here incoming edge is 0 by default

     // disappearance cost
     auto min_disappearance_cost = cost01(disappearance_begin(), smallest_incoming, current_arc_cost);

     const REAL cost_0 = std::min({min_detection_cost[0], min_appearance_cost, min_disappearance_cost[0]});
     const REAL cost_1 = std::min(min_detection_cost[1], min_disappearance_cost[1]);
     return {cost_0, cost_1};
  }

  REAL EvaluatePrimal() const 
  {
     return std::numeric_limits<REAL>::infinity();
  }
  INDEX incoming_edge_, outgoing_edge_; // must be vector of incoming and outgoing edges 
private:
  const INDEX max_detections_, no_incoming_edges_, no_outgoing_edges_;
  REAL* pot_;
  bool divides_; 
};

// transition message between multiple_detection_factor
class transition_message_multiple {
public:
   transition_message_multiple(const bool split, const INDEX outgoing_edge_index, const INDEX incoming_edge_index) 
      : outgoing_edge_index_(outgoing_edge_index),
      incoming_edge_index_(incoming_edge_index),
      split_(split)
   {}


   template<typename RIGHT_FACTOR, typename MSG>
   void ReceiveMessageFromRight(const RIGHT_FACTOR& r, MSG& msg)
   {
     auto cost = r.incoming_cost(incoming_edge_index_);
     msg[0] -= cost[1] - cost[0];
   }

   template<typename LEFT_FACTOR, typename MSG>
   void ReceiveMessageFromLeft(const LEFT_FACTOR& l, MSG& msg)
   {
     auto cost = l.outgoing_cost(outgoing_edge_index_);
     msg[0] -= cost[1] - cost[0];
   }

   template<typename LEFT_FACTOR, typename MSG>
   void SendMessageToRight(const LEFT_FACTOR& l, MSG msg, const REAL omega)
   {
     auto cost = l.outgoing_cost(outgoing_edge_index_);
     msg[0] -= omega*(cost[1] - cost[0]);
   
   }
   template<typename RIGHT_FACTOR, typename MSG>
   void SendMessageToLeft(const RIGHT_FACTOR& r, MSG msg, const REAL omega)
   {
     auto cost = r.incoming_cost(incoming_edge_index_);
     msg[0] -= omega*(cost[1] - cost[0]);
   }

   // send messages from detection factor along incoming edges
   template<typename LEFT_FACTOR, typename MSG_ARRAY, typename ITERATOR>
   static void tmp_SendMessagesToRight(const LEFT_FACTOR& leftFactor, MSG_ARRAY msg_begin, MSG_ARRAY msg_end, ITERATOR omegaIt)
   {
     return;
      //std::cout << "send to right:";
      //for(INDEX i=0; i<leftFactor.size(); ++i) std::cout << leftFactor[i] << ",";
      //std::cout << "\n";

      REAL detection_incoming_cost = leftFactor[0];
      detection_incoming_cost += *std::min_element(leftFactor.incoming_begin(), leftFactor.incoming_end());

      REAL omega = 0.0;
      for(auto it = msg_begin; it!=msg_end; ++it, ++omegaIt) {
         omega += *omegaIt;
      }
      assert(omega <= 1.0 + eps);
      //omega = 1.0;

      // compute smallest and second smallest value over all outgoing edges
      REAL smallest_outgoing = std::numeric_limits<REAL>::infinity();
      REAL second_smallest_outgoing = std::numeric_limits<REAL>::infinity();
      for(auto it=leftFactor.outgoing_begin(); it!=leftFactor.outgoing_end(); ++it) {
         if(*it <= smallest_outgoing) {
            second_smallest_outgoing = smallest_outgoing;
            smallest_outgoing = *it;
         } else if(*it <= second_smallest_outgoing) {
            second_smallest_outgoing = *it;
         }
      }

      // do zrobienia: check whether no of messages is greater than no of outgoing edges
      INDEX c=0;
      for(auto it = msg_begin; it!=msg_end; ++it)  ++c;
      assert(c+1 >= std::distance(leftFactor.outgoing_begin(), leftFactor.outgoing_end()));

      for(; msg_begin!=msg_end; ++msg_begin) {
         if((*msg_begin).GetMessageOp().split_) {
            //(*msg_begin)[0] -= 0.5*omega*(detection_incoming_cost + leftFactor.outgoing((*msg_begin).GetMessageOp().outgoing_edge_index_));
            (*msg_begin)[0] -= 0.5*omega*( leftFactor.outgoing( (*msg_begin).GetMessageOp().outgoing_edge_index_ ) - std::min(-detection_incoming_cost, second_smallest_outgoing) );
         } else {
            //(*msg_begin)[0] -= omega*(std::min(detection_incoming_cost + leftFactor[(*msg_begin).GetMessageOp().outgoing_edge_index_], 0.0));
            //(*msg_begin)[0] -= omega*(detection_incoming_cost + leftFactor.outgoing((*msg_begin).GetMessageOp().outgoing_edge_index_));
            (*msg_begin)[0] -= omega*( leftFactor.outgoing( (*msg_begin).GetMessageOp().outgoing_edge_index_ ) - std::min(-detection_incoming_cost, second_smallest_outgoing) );
         }
      }
   }

  template<typename G>
  void RepamLeft(G& l, const REAL msg, const INDEX msg_dim)
  {
    assert(msg_dim == 0);
    l.outgoing(outgoing_edge_index_) += msg;
  }
  template<typename G>
  void RepamRight(G& r, const REAL msg, const INDEX msg_dim)
  {
    assert(msg_dim == 0);
    r.incoming(incoming_edge_index_) += msg;
  }

private:
  const INDEX outgoing_edge_index_;
  const INDEX incoming_edge_index_;
  const bool split_; // is split really needed? Whether it is a split transition can be found out in SendMessagesTo... by checking whether there is two outgoing_edge_index 
};



} // end namespace LP_MP

#endif // LP_MP_DETECTION_FACTOR_HXX
