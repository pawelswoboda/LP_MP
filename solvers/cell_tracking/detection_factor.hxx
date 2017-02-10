#ifndef LP_MP_DETECTION_FACTOR_HXX
#define LP_MP_DETECTION_FACTOR_HXX

#include "config.hxx"
#include "vector.hxx"
#include "cereal/types/array.hpp"

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

  detection_factor(const INDEX no_incoming_edges, const INDEX no_outgoing_edges, const REAL detection_cost, bool can_exit = true)
    : no_incoming_edges_(no_incoming_edges),
    no_outgoing_edges_(no_outgoing_edges+can_exit) // the last outgoing edge is the exit one (when exit is possible)
  {
    assert(can_exit == true); // true for mother machine only
    assert(no_outgoing_edges+can_exit == no_outgoing_edges + 1);
    pot_ = global_real_block_allocator.allocate(size());
    assert(pot_ != nullptr);
    std::fill(pot_+1, pot_ + size(), 0.0);
    pot_[0] = detection_cost;
  }
  ~detection_factor() {
    global_real_block_allocator.deallocate(pot_,1);
  }
  detection_factor(const detection_factor& o) 
    : no_incoming_edges_(o.no_incoming_edges_),
    no_outgoing_edges_(o.no_outgoing_edges_)
  {
    pot_ = global_real_block_allocator.allocate(size());
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

    REAL lb = pot_[0];
    if(no_incoming_edges_ > 0) lb += *std::min_element(incoming_begin(),incoming_end()); 
    if(no_outgoing_edges_ > 0) lb += *std::min_element(outgoing_begin(),outgoing_end()); 
    //std::cout << lb << "\n";
    return lb;
  }

  void MaximizePotentialAndComputePrimal() {
    if(incoming_edge_ == no_edge_taken) {
      outgoing_edge_ = no_edge_taken;
      return;
    }
    if(outgoing_edge_ == no_edge_taken) {
      incoming_edge_ = no_edge_taken;
      return;
    }

    INDEX incoming_cand = 0;
    REAL lb = pot_[0];
    if(no_incoming_edges() > 0) {
      incoming_cand = std::min_element(incoming_begin(), incoming_end()) - incoming_begin();
      lb += incoming(incoming_cand); 
    }
    INDEX outgoing_cand = 0;
    if(no_outgoing_edges() > 0) {
      outgoing_cand = std::min_element(outgoing_begin(), outgoing_end()) - outgoing_begin();
      lb += outgoing(outgoing_cand); 
    }
    // one edge already labelled
    if(incoming_edge_ < no_incoming_edges() && no_outgoing_edges() > 0) {
      assert(outgoing_edge_ == no_primal_decision);
      outgoing_edge_ = outgoing_cand;
      return;
    } else if(outgoing_edge_ < no_outgoing_edges() && no_incoming_edges() > 0) {
      assert(incoming_edge_ == no_primal_decision);
      incoming_edge_ = incoming_cand;
      return;
    }
    // no edge labelled yet
    assert(incoming_edge_ == no_primal_decision && outgoing_edge_ == no_primal_decision);
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
    if(no_outgoing_edges() > 0 && no_incoming_edges() > 0) {
      if(incoming_edge_ < no_incoming_edges() && outgoing_edge_ < no_outgoing_edges()) {
        return pot_[0] + incoming(incoming_edge_) + outgoing(outgoing_edge_);
      } else {
        return 0.0;
      }
    }
    if(no_incoming_edges() == 0) {
      if(outgoing_edge_ < no_outgoing_edges()) {
        return pot_[0] + outgoing(outgoing_edge_);
      } else {
        return 0.0;
      }
    }
    assert(no_outgoing_edges() > 0);
    assert(false);
    return std::numeric_limits<REAL>::infinity();
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
private:
  const INDEX no_incoming_edges_, no_outgoing_edges_;
  REAL* pot_;
  
  INDEX incoming_edge_, outgoing_edge_;
};

// detection can contain multiple cells
// layout of reparametrization: detection costs, appearance costs, disappearance costs, incoming edge costs, outgoing edge costs
class multiple_detection_factor {
public:
   multiple_detection_factor(const INDEX max_detections, const INDEX no_incoming_edges, const INDEX no_outgoing_edges)
      : max_detections_(max_detections),
        no_incoming_edges_(no_incoming_edges),
        no_outgoing_edges_(no_outgoing_edges) // the last outgoing edge is the exit one (when exit is possible)
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

  void sort_incoming_cost(vector& v) const
  {
     assert(v.size() == std::min(no_incoming_edges(), max_detections_));
     std::partial_sort_copy(incoming_begin(), incoming_end(), v.begin(), v.end());
     std::partial_sum(v.begin(), v.end(), v.begin());
  }

  void sort_outgoing_cost(vector& v) const
  {
     assert(v.size() == std::min(no_outgoing_edges(), max_detections_));
     std::partial_sort_copy(outgoing_begin(), outgoing_end(), v.begin(), v.end());
     std::partial_sum(v.begin(), v.end(), v.begin()); 
  }

  REAL detection_cost(const vector& smallest_incoming, const vector& smallest_outgoing) const 
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

  REAL appearance_cost(const vector& smallest_outgoing) const 
  {
    REAL min_appearance_cost = std::numeric_limits<REAL>::infinity();
    for(INDEX no_detections=1; no_detections<=max_detections_; ++no_detections) {
      REAL appearance_cost = appearance(no_detections-1) + outgoing(no_detections-1);
      min_appearance_cost = std::min(appearance_cost, min_appearance_cost);
    }
    return min_appearance_cost; 
  }

  REAL disappearance_cost(const vector& smallest_incoming) const 
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
     vector smallest_incoming( std::min(no_incoming_edges(), max_detections_) );
     sort_incoming_cost(smallest_incoming);
     vector smallest_outgoing( std::min(no_outgoing_edges(), max_detections_) );
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

     vector smallest_incoming( std::min(no_incoming_edges(), max_detections_) );
     sort_incoming_cost(smallest_incoming);

     vector smallest_outgoing( std::min(no_outgoing_edges(), max_detections_+1) );
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

     vector smallest_incoming( std::min(no_incoming_edges(), max_detections_) );
     sort_incoming_cost(smallest_incoming);

     vector smallest_outgoing( std::min(no_outgoing_edges(), max_detections_+1) );
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
  static void SendMessagesToLeft(const RIGHT_FACTOR& rightFactor, MSG_ARRAY msg_begin, MSG_ARRAY msg_end, ITERATOR omegaIt)
  {
    //std::cout << "send to left :";
    //for(INDEX i=0; i<rightFactor.size(); ++i) std::cout << rightFactor[i] << ",";
    //std::cout << "\n";
    REAL detection_outgoing_cost = rightFactor[0];
    if(rightFactor.no_outgoing_edges() != 0) {
      detection_outgoing_cost += *std::min_element(rightFactor.outgoing_begin(), rightFactor.outgoing_end());
    }

    REAL omega = 0.0;
    for(auto it = msg_begin; it!=msg_end; ++it, ++omegaIt) {
      omega += *omegaIt;
    }
    assert(omega <= 1.0 + eps);
    assert(omega > 0);

    // compute smallest and second smallest value over all incoming edges
    auto smallest_incoming = two_smallest_elements<REAL>(rightFactor.incoming_begin(), rightFactor.incoming_end());
    //std::cout << "smallest = " << smallest_incoming << ", second smallest = " << second_smallest_incoming << "\n"; 

    // do zrobienia: check whether no of messages is greater than no of incoming edges
    //INDEX c=0;
    //for(auto it = msg_begin; it!=msg_end; ++it)  ++c;
    //assert(c == std::distance(rightFactor.incoming_begin(), rightFactor.incoming_end()));

    for(; msg_begin!=msg_end; ++msg_begin) {
      //std::cout << (*msg_begin).GetMessageOp().incoming_edge_index_  << ",";
      (*msg_begin)[0] -= omega*( rightFactor.incoming( (*msg_begin).GetMessageOp().incoming_edge_index_ ) - std::min(-detection_outgoing_cost, smallest_incoming[1]) );
    }
    //std::cout << "\n";
  }
  // send messages from detection factor along incoming edges
  template<typename LEFT_FACTOR, typename MSG_ARRAY, typename ITERATOR>
  static void SendMessagesToRight(const LEFT_FACTOR& leftFactor, MSG_ARRAY msg_begin, MSG_ARRAY msg_end, ITERATOR omegaIt)
  {
    //std::cout << "send to right:";
    //for(INDEX i=0; i<leftFactor.size(); ++i) std::cout << leftFactor[i] << ",";
    //std::cout << "\n";

    REAL detection_incoming_cost = leftFactor[0];
    if(leftFactor.no_incoming_edges() > 0) {
      detection_incoming_cost += *std::min_element(leftFactor.incoming_begin(), leftFactor.incoming_end());
    }

    REAL omega = 0.0;
    for(auto it = msg_begin; it!=msg_end; ++it, ++omegaIt) {
      omega += *omegaIt;
    }
    assert(omega <= 1.0 + eps);
    //omega = 1.0;

    // compute smallest and second smallest value over all outgoing edges
    auto smallest_outgoing = two_smallest_elements<REAL>(leftFactor.outgoing_begin(), leftFactor.outgoing_end());

    // do zrobienia: check whether no of messages is greater than no of outgoing edges
    //INDEX c=0;
    //for(auto it = msg_begin; it!=msg_end; ++it)  ++c;
    //assert(c+1 >= std::distance(leftFactor.outgoing_begin(), leftFactor.outgoing_end()));

    for(; msg_begin!=msg_end; ++msg_begin) {
      if((*msg_begin).GetMessageOp().split_) {
        //(*msg_begin)[0] -= 0.5*omega*(detection_incoming_cost + leftFactor.outgoing((*msg_begin).GetMessageOp().outgoing_edge_index_));
        (*msg_begin)[0] -= 0.5*omega*( leftFactor.outgoing( (*msg_begin).GetMessageOp().outgoing_edge_index_ ) - std::min(-detection_incoming_cost, smallest_outgoing[1]) );
      } else {
        //(*msg_begin)[0] -= omega*(std::min(detection_incoming_cost + leftFactor[(*msg_begin).GetMessageOp().outgoing_edge_index_], 0.0));
        //(*msg_begin)[0] -= omega*(detection_incoming_cost + leftFactor.outgoing((*msg_begin).GetMessageOp().outgoing_edge_index_));
        (*msg_begin)[0] -= omega*( leftFactor.outgoing( (*msg_begin).GetMessageOp().outgoing_edge_index_ ) - std::min(-detection_incoming_cost, smallest_outgoing[1]) );
      }
    }
  }

  template<typename LEFT_FACTOR, typename RIGHT_FACTOR>
  bool ComputeRightFromLeftPrimal(const LEFT_FACTOR& l, RIGHT_FACTOR& r)
  {
    return false;
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
    return false;
    if(r.incoming_edge_ == incoming_edge_index_ && r.outgoing_edge_ != outgoing_edge_index_) {
      l.outgoing_edge_ = outgoing_edge_index_;
      return true;
    } else {
      return false;
    }
  }

  template<typename RIGHT_FACTOR, typename G2>
  void make_right_factor_uniform(RIGHT_FACTOR& r, G2& msg)
  {
    REAL detection_outgoing_cost = r[0];
    if(r.no_outgoing_edges() > 0) {
      detection_outgoing_cost += *std::min_element(r.outgoing_begin(), r.outgoing_end());
    }

    REAL min_incoming_val = std::numeric_limits<REAL>::infinity();
    INDEX c=0;
    for(auto it = r.incoming_begin(); it!=r.incoming_end(); ++it, ++c) {
      if(c != incoming_edge_index_) {
        min_incoming_val = std::min(min_incoming_val, *it);
      } 
    }
    msg[0] -= std::min(detection_outgoing_cost + r.incoming(incoming_edge_index_) , 0.0)
      - std::min(detection_outgoing_cost + min_incoming_val, 0.0); // or +- exchanged 
  }
  template<typename RIGHT_FACTOR, typename G2>
  void ReceiveMessageFromRight(RIGHT_FACTOR& r, G2& msg)
  { 
    make_right_factor_uniform(r,msg);
  }
  template<typename RIGHT_FACTOR, typename G2>
  void ReceiveRestrictedMessageFromRight(RIGHT_FACTOR& r, G2& msg)
  { 
    return;
    if(r.incoming_edge_ < r.no_incoming_edges()) { // incoming edge has already been picked
      if(r.incoming_edge_ == incoming_edge_index_) {
        msg[0] -= -std::numeric_limits<REAL>::infinity(); 
        assert(msg.GetLeftFactor()->GetFactor()->outgoing(outgoing_edge_index_) < -10000);
      } else {
        msg[0] -= std::numeric_limits<REAL>::infinity(); 
        assert(msg.GetLeftFactor()->GetFactor()->outgoing(outgoing_edge_index_) > 10000);
      }
      //assert(r.incoming_edge_ != incoming_edge_index_);
    } else { // no incoming edge has been picked yet.
      return;
      make_right_factor_uniform(r,msg); 
    }
  }

  template<typename LEFT_FACTOR, typename G2>
  void make_left_factor_uniform(LEFT_FACTOR& l, G2& msg)
  { 
    const REAL detection_incoming_cost = l[0] + *std::min_element(l.incoming_begin(), l.incoming_end());

    REAL min_outgoing_val = std::numeric_limits<REAL>::infinity();
    INDEX c=0;
    for(auto it = l.outgoing_begin(); it!=l.outgoing_end(); ++it, ++c) {
      if(c != outgoing_edge_index_) {
        min_outgoing_val = std::min(min_outgoing_val, *it);
      } 
    }
    msg[0] -= std::min(detection_incoming_cost + l.outgoing(outgoing_edge_index_),  0.0)
      - std::min(detection_incoming_cost + min_outgoing_val,  0.0); // or +- exchanged
  }

  template<typename LEFT_FACTOR, typename G2>
  void ReceiveMessageFromLeft(LEFT_FACTOR& l, G2& msg)
  {
    make_left_factor_uniform(l,msg);
  }
  template<typename LEFT_FACTOR, typename G2>
  void ReceiveRestrictedMessageFromLeft(LEFT_FACTOR& l, G2& msg)
  { 
    return;
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
      make_left_factor_uniform(l,msg); 
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
  const bool split_; // is split really needed? Whether it is a split transition can be found out in SendMessagesTo... by checking whether there is two outgoing_edge_indices
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
      if(leftFactor.no_incoming_edges() > 0) {
         detection_incoming_cost += *std::min_element(leftFactor.incoming_begin(), leftFactor.incoming_end());
      }

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

// multiple cell detection hypotheses can be mutually exclusive. 
// simplex x1 + ... + xn <= 1
// to account for overlapping detections: only one can be active
class at_most_one_cell_factor : public vector {

  friend class at_most_one_cell_message;

public:

  constexpr static INDEX no_primal_decision = std::numeric_limits<INDEX>::max();
  constexpr static INDEX no_primal_active = std::numeric_limits<INDEX>::max()-1;
  constexpr static INDEX primal_infeasible = std::numeric_limits<INDEX>::max()-2;

   at_most_one_cell_factor(const INDEX size) : vector(size)
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
       return 1000;
       assert(false);
       return std::numeric_limits<REAL>::infinity();
     }
   }

   void init_primal() { primal_ = no_primal_decision; }
   template<class ARCHIVE> void serialize_primal(ARCHIVE& ar) { ar(primal_); }
   template<class ARCHIVE> void serialize_dual(ARCHIVE& ar) { ar( *static_cast<vector*>(this) ); }
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
    } 
  }

  template<typename LEFT_FACTOR, typename MSG>
  void SendMessageToRight(const LEFT_FACTOR& l, MSG& msg, const REAL omega)
  { 
    make_left_factor_uniform(l, msg, omega);
  }

  template<typename LEFT_FACTOR, typename MSG_ARRAY, typename ITERATOR>
  static void SendMessagesToRight(const LEFT_FACTOR& l, MSG_ARRAY msg_begin, MSG_ARRAY msg_end, ITERATOR omegaIt)
  {
    std::cout << "kwas2\n";
    const REAL delta = l.cost_of_detection();
    assert(!std::isnan(delta));

    for(; msg_begin!=msg_end; ++msg_begin, ++omegaIt) {
      (*msg_begin)[0] -= (*omegaIt)*delta;
    } 
  }


  template<typename RIGHT_FACTOR, typename MSG>
  void ReceiveMessageFromRight(const RIGHT_FACTOR& r, MSG& msg) 
  {
    make_right_factor_uniform(r,msg); 
  }

  template<typename RIGHT_FACTOR, typename MSG>
  void ReceiveRestrictedMessageFromRight(const RIGHT_FACTOR& r, MSG& msg) 
  {
    assert(at_most_one_cell_factor_index_ < r.size());
    if(r.primal_ == r.size()) { // no element chosen yet
      //make_right_factor_uniform(r, msg);
    } else if(r.primal_ < r.size()) { // one element already chosen
      if(r.primal_ == at_most_one_cell_factor_index_) {
        msg[0] -= -100000;//std::numeric_limits<REAL>::infinity();
        assert(msg.GetLeftFactor()->GetFactor()->operator[](0) < -10000);
      } else {
        //assert(r.primal_ < r.size());
        const REAL test_val_prev = msg.GetLeftFactor()->GetFactor()->operator[](0);
        msg[0] -= 100000;//std::numeric_limits<REAL>::infinity();
        const REAL test_val = msg.GetLeftFactor()->GetFactor()->operator[](0);
        assert(test_val > 10000);
        assert(msg.GetLeftFactor()->GetFactor()->operator[](0) > 10000);
      }
    } else {
        msg[0] -= 100000;//std::numeric_limits<REAL>::infinity();
        const REAL test_val = msg.GetLeftFactor()->GetFactor()->operator[](0);
        assert(test_val > 10000);
        assert(msg.GetLeftFactor()->GetFactor()->operator[](0) > 10000); 
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
  void make_left_factor_uniform(const LEFT_FACTOR& l, G2& msg, const REAL omega = 1.0)
  {
    // make cost of detection same as cost of non-detection
    msg[0] -= omega*l.cost_of_detection();
  }
  template<typename RIGHT_FACTOR, typename G2>
  void make_right_factor_uniform(const RIGHT_FACTOR& r, G2& msg, const REAL omega = 1.0)
  {
    // compute cost without current factor,
    auto smallest = two_smallest_elements<REAL>(r.begin(), r.end());
    const REAL cur_detection_cost = r[at_most_one_cell_factor_index_];
    REAL rest_cost;
    if(smallest[0] == cur_detection_cost) {
      rest_cost = std::min(0.0, smallest[1]);
    } else {
      rest_cost = std::min(0.0, smallest[0]);
    }

    msg[0] -= omega*(cur_detection_cost - rest_cost);
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
  using repam_type = std::array<REAL,2>;
  REAL LowerBound() const 
  {
    return std::min(0.0, std::min( (*this)[0], (*this)[1] )); 
  }
  REAL EvaluatePrimal() const 
  {
    if(primal_[0] && !primal_[1]) {
      return (*this)[0];
    } else if(primal_[1] && !primal_[0]) {
      return (*this)[0];
    } else if(!primal_[0] && !primal_[1]) {
      return 0.0;
    } else {
      //assert(false);
      return std::numeric_limits<REAL>::infinity();
    }
  }

  void init_primal() { primal_[0] = true; primal_[1] = true; }
  template<class ARCHIVE> void serialize_primal(ARCHIVE& ar) { ar( primal_ ); }
  template<class ARCHIVE> void serialize_dual(ARCHIVE& ar) { ar( *static_cast<repam_type*>(this) ); }
private:
  std::array<bool,2> primal_;
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
        r.primal_[0] = true;
      } 
    } else {

    }
  }

  template<typename LEFT_FACTOR, typename MSG>
  void SendMessageToRight(const LEFT_FACTOR& l, MSG& msg, const REAL omega)
  { 
    REAL detection_incoming_cost = l[0];
    if(l.no_incoming_edges() > 0) {
      detection_incoming_cost += *std::min_element(l.incoming_begin(), l.incoming_end());
    }
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
    if(l.no_incoming_edges() > 0) {
      detection_incoming_cost += *std::min_element(l.incoming_begin(), l.incoming_end());
    }
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
  void ReceiveMessageFromRight(const RIGHT_FACTOR& r, MSG& msg) 
  {
    assert(r.size() == 2);
    if(POSITION == exit_constraint_position::lower) {
      msg[0] -= r[0] - std::min(0.0, r[1]);
    } else {
      msg[0] -= r[1] - std::min(0.0, r[0]);
    }
  }
};

// obsolete factors and messages //
/*
// models either a possible split or pure transition. Used so that messages can be sent in batch from outgoing edges in detection factor
// variant whether just simple transition factor or split one
class transition_factor {
public:
  transition_factor(const REAL cost) : pot_(cost) {}
  REAL LowerBound() const {
    return std::min(0.0, REAL());
  }
  REAL EvaluatePrimal(PrimalSolutionStorage::Element primal) const {
    assert(*primal != unknownState);
    return primal[0]*pot_;
  }
  static constexpr REAL size() { return 1; }
  operator REAL() { return pot_; }
private:
    REAL pot_;
};

// possible states: 000, 110, 101, 111, the first of which has always cost zero, the last three are held explicitly
class split_factor : public std::array<REAL,3> {
public:
  split_factor() {}
  REAL LowerBound() const {
    return std::min(*std::min_element(this->begin(), this->end()), 0.0);
  }
  REAL EvaluatePrimal(PrimalSolutionStorage::Element primal) const {
    if(std::accumulate(primal, primal + size(),0) > 1) {
      return std::numeric_limits<REAL>::infinity();
    } else {
      REAL cost = 0.0;
      for(auto it=this->begin(); it!=this->end(); ++it, ++primal) {
        cost += (*primal)*(*it);
      }
      return cost;
    }
  }
};


// Ternary factor on edges in tracking. Models when cell splits.
// left factor is detection factor, right one is transition factor
class outgoing_edge_message {
public:
  outgoing_edge_message(const INDEX outgoing_edge_index) : outgoing_edge_index_(outgoing_edge_index) {}

  template<typename G>
  void RepamLeft(G& l, const REAL msg, const INDEX msg_dim)
  {
    assert(msg_dim == 0);
    l[outgoing_edge_index_] += msg;
  }
  template<typename G>
  void RepamRight(G& r, const REAL msg, const INDEX msg_dim)
  {
    assert(msg_dim == 0);
    r[0] += msg; // do zrobienia: will not work in general, when split factor is on right
  }

  template<typename MSG_ARRAY, typename LEFT_FACTOR>
  static void make_left_factor_uniform_parallel(MSG_ARRAY msg_begin, MSG_ARRAY msg_end, const LEFT_FACTOR& l, const REAL omega)
  {
    // assert that number of messages is number of outgoing edges
    const REAL detection_incoming_cost = l[0] + *std::min_element(l.incoming_begin(), l.incoming_end());

    for(; msg_begin!=msg_end; ++msg_begin) {
      (*msg_begin)[0] -= omega*(std::min(detection_incoming_cost + l[(*msg_begin).outgoing_edge_index_], 0.0));
    }
  }
  template<typename RIGHT_FACTOR, typename G2>
  void make_right_factor_uniform(const RIGHT_FACTOR& r, G2& msg, const REAL omega = 1.0)
  {
    msg[0] -= omega*(*r);
  }

private:
  const INDEX outgoing_edge_index_;
};

// left factor is transition_factor, right one is cell_detection_factor
class incoming_edge_message {
public:
  incoming_edge_message(const INDEX incoming_edge_index) 
    : incoming_edge_index_(incoming_edge_index) 
  {}

  template<typename G>
  void RepamLeft(G& l, const REAL msg, const INDEX msg_dim)
  {
    assert(msg_dim == 0);
    *l += msg;
  }
  template<typename G>
  void RepamRight(G& r, const REAL msg, const INDEX msg_dim)
  {
    assert(msg_dim == 0);
    r[incoming_edge_index_] += msg;
  }

  template<typename LEFT_FACTOR, typename MSG>
  static void make_left_factor_uniform(const LEFT_FACTOR& l, MSG& msg, const REAL omega = 1.0) 
  {
    msg[0] -= omega*(*l);
  }
  template<typename MSG_ITERATOR, typename RIGHT_FACTOR>
  void make_right_factor_uniform_parallel(MSG_ITERATOR msg_begin, MSG_ITERATOR msg_end, const RIGHT_FACTOR& r, const REAL omega = 1.0)
  {
    // assert that number of messages is number of outgoing edges
    const REAL detection_outgoing_cost = r[0] + *std::min_element(r.outgoing_begin(), r.outgoing_end());

    for(; msg_begin!=msg_end; ++msg_begin) {
      (*msg_begin)[0] -= omega*(std::min(detection_outgoing_cost + r[(*msg_begin).incoming_edge_index_], 0.0));
    }
  }
private:
  const INDEX incoming_edge_index_;
};

// for cell divisions
// outgoing transition
// left is transition factor, right is division factor
class transition_split_message_out {
public:
  template<typename G>
  void RepamLeft(G& l, const REAL msg, const INDEX msg_dim)
  {
    assert(msg_dim == 0);
    l += msg;
  }
  template<typename G>
  void RepamRight(G& r, const REAL msg, const INDEX msg_dim)
  {
    assert(msg_dim == 0);
    assert(false);
    r[2] += msg;
  }

  template<typename LEFT_FACTOR, typename MSG>
  static void make_left_factor_uniform(const LEFT_FACTOR& l, MSG& msg, const REAL omega = 1.0) 
  {
    msg[0] -= omega*(*l);
  }
  template<typename RIGHT_FACTOR, typename MSG>
  static void make_right_factor_uniform(const RIGHT_FACTOR& r, MSG& msg, const REAL omega = 1.0) 
  {
    assert(false);
  }
};

// left is division factor, right is transition factor
// incoming transition
class transition_split_message_in {
public:
  transition_split_message_in(const char index) : index_(index) 
  {
    assert(index_ < 2);
  }

  template<typename G>
  void RepamLeft(G& l, const REAL msg, const INDEX msg_dim)
  {
    assert(msg_dim == 0);
    l[index_] += msg;
  }
  template<typename G>
  void RepamRight(G& r, const REAL msg, const INDEX msg_dim)
  {
    assert(msg_dim == 0);
    assert(false);
    r += msg;
  }

  template<typename LEFT_FACTOR, typename MSG>
  static void make_left_factor_uniform(const LEFT_FACTOR& l, MSG& msg, const REAL omega = 1.0) 
  {
    assert(false);
  }
  template<typename RIGHT_FACTOR, typename MSG>
  static void make_right_factor_uniform(const RIGHT_FACTOR& r, MSG& msg, const REAL omega = 1.0) 
  {
    msg[0] -= omega*(*r);
  }

private:
  const char index_; 
};
*/

} // end namespace LP_MP

#endif // LP_MP_DETECTION_FACTOR_HXX
