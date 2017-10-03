#ifndef LP_MP_CELL_TRACKING_CONSTRUCTOR_HXX
#define LP_MP_CELL_TRACKING_CONSTRUCTOR_HXX

#include <vector>
#include "solver.hxx"
#include "pegtl/parse.hh"
#include "parse_rules.h"

// do zrobienia: do not include detection hypotheses that have no incoming and no outgoing edges into the LP

namespace LP_MP {

  struct exclusion_item { INDEX timestep; INDEX hypothesis_id;
    INDEX operator[](const INDEX i) const { if(i==0) return timestep; else return hypothesis_id;}
  };
  bool operator<=(const exclusion_item a, const exclusion_item b) 
  {
    if(a.timestep != b.timestep) {
      return a.timestep < b.timestep;
    } else {
      return a.hypothesis_id < b.hypothesis_id;
    } 
  }
  bool operator==(const exclusion_item a, const exclusion_item b) 
  {
    return a.timestep == b.timestep && a.hypothesis_id == b.hypothesis_id;
  }
  bool operator!=(const exclusion_item a, const exclusion_item b) 
  {
    return !(a == b);
  }
  static constexpr exclusion_item exclusion_item_delimiter = {std::numeric_limits<INDEX>::max(), std::numeric_limits<INDEX>::max()};

  // temporary structure which counts how many incoming and outgoing edges are already used by messages for building the model
  struct cell_tracking_transition_count {
    //transition_count(CONSTRUCTOR& c) 
    //{
    //  current_transition_no.resize( c.detection_factors_.size() );
    //  current_division_no.resize( c.detection_factors_.size() );
    //  for(INDEX i=0; i<c.detection_factors_.size(); ++i) {
    //    current_transition_no[i].resize( c.detection_factors_[i].size(), {0,0});
    //    current_division_no[i].resize( c.detection_factors_[i].size(), {0,0});
    //  }
    //}

    struct edge_type {
      std::array<INDEX,2> transition;
      std::array<INDEX,2> division;
    };
    std::vector<std::vector<edge_type>> edges;

    void add_detection_hypothesis(
        const INDEX timestep, const INDEX hypothesis_id,
        const INDEX no_incoming_transition_edges, const INDEX no_incoming_division_edges,
        const INDEX no_outgoing_transition_edges, const INDEX no_outgoing_division_edges
        )
    {
      if(timestep >= edges.size()) {
        edges.resize(timestep+1);
      }
      if(hypothesis_id >= edges[timestep].size()) {
        edges[timestep].resize(hypothesis_id + 1, {0,0});
      }
      assert(edges[timestep][hypothesis_id].transition[0] == 0 && edges[timestep][hypothesis_id].transition[1] == 0);
      assert(edges[timestep][hypothesis_id].division[0] == 0 && edges[timestep][hypothesis_id].division[1] == 0);

      edges[timestep][hypothesis_id].division[0] = no_incoming_transition_edges;
      edges[timestep][hypothesis_id].division[1] = no_outgoing_transition_edges;
    }
    INDEX next_incoming_transition_edge(const INDEX timestep, const INDEX hypothesis_id)
    {
      assert(timestep < edges.size() && hypothesis_id < edges[timestep].size());
      return edges[timestep][hypothesis_id].transition[0]++; 
    }
    INDEX next_incoming_division_edge(const INDEX timestep, const INDEX hypothesis_id)
    {
      assert(timestep < edges.size() && hypothesis_id < edges[timestep].size());
      return edges[timestep][hypothesis_id].division[0]++; 
    }
    INDEX next_outgoing_transition_edge(const INDEX timestep, const INDEX hypothesis_id)
    {
      assert(timestep < edges.size() && hypothesis_id < edges[timestep].size());
      return edges[timestep][hypothesis_id].transition[1]++; 
    }
    INDEX next_outgoing_division_edge(const INDEX timestep, const INDEX hypothesis_id)
    {
      assert(timestep < edges.size() && hypothesis_id < edges[timestep].size());
      return edges[timestep][hypothesis_id].division[1]++; 
    }
  };

template<typename DETECTION_FACTOR_CONTAINER, typename AT_MOST_ONE_CELL_FACTOR_CONTAINER, typename AT_MOST_ONE_CELL_MESSAGE_CONTAINER>
class basic_cell_tracking_constructor {
public:
  using detection_factor_container = DETECTION_FACTOR_CONTAINER;
  using at_most_one_cell_factor_container = AT_MOST_ONE_CELL_FACTOR_CONTAINER;
  using exclusion_factor = AT_MOST_ONE_CELL_FACTOR_CONTAINER;

  using CONSTRUCTOR = basic_cell_tracking_constructor<DETECTION_FACTOR_CONTAINER, AT_MOST_ONE_CELL_FACTOR_CONTAINER, AT_MOST_ONE_CELL_MESSAGE_CONTAINER>;

  template<typename SOLVER>
  basic_cell_tracking_constructor(SOLVER& solver) 
  : lp_(&solver.GetLP()) 
  {}

  //transition_count init_transition_counter() { return transition_count(*this); }

  void set_number_of_timesteps(const INDEX t)
  {
    assert(detection_factors_.size() == 0);
    detection_factors_.resize(t);
  }

  template<typename LP_TYPE>
  DETECTION_FACTOR_CONTAINER* 
  add_detection_hypothesis(
      LP_TYPE& lp, 
      const INDEX timestep, const INDEX hypothesis_id, 
      const REAL detection_cost, const REAL appearance_cost, const REAL disappearance_cost, 
      const INDEX no_incoming_transition_edges, const INDEX no_incoming_division_edges, 
      const INDEX no_outgoing_transition_edges, const INDEX no_outgoing_division_edges
      )
  { 
    assert(timestep < detection_factors_.size());
    if(hypothesis_id >= detection_factors_[timestep].size()) {
      detection_factors_[timestep].resize(hypothesis_id+1, nullptr);
    }
    assert(detection_factors_[timestep][hypothesis_id] == nullptr);

    // possible we should everywhere remove the possiblity of an empty factor
    //if(detection_cost == 0.0 && no_incoming_transition_edges == 0 && no_outgoing_edges == 0 && appearance_cost == 0.0 && disappearance_cost == 0.0) { return nullptr; }

    auto* f = new DETECTION_FACTOR_CONTAINER(no_incoming_transition_edges, no_incoming_division_edges, no_outgoing_transition_edges, no_outgoing_division_edges, detection_cost, appearance_cost, disappearance_cost); 
    lp.AddFactor(f);
    detection_factors_[timestep][hypothesis_id] = f;
    if(hypothesis_id > 0) {
      assert( detection_factors_[timestep][hypothesis_id-1] != nullptr); // need not be generally true, but then factor relation must be done more robust.
    }
    //std::cout << "H: " << timestep << ", " << hypothesis_id <<  "," << no_incoming_edges << "," << no_outgoing_edges << ", " << detection_cost << ", " << appearance_cost << ", " << disappearance_cost << std::endl;

    tc_.add_detection_hypothesis(timestep, hypothesis_id, no_incoming_transition_edges, no_incoming_division_edges, no_outgoing_transition_edges, no_outgoing_division_edges);

    return f; 
  }

  template<typename LP_TYPE, typename ITERATOR>
  AT_MOST_ONE_CELL_FACTOR_CONTAINER* add_exclusion_constraint(LP_TYPE& lp, const INDEX timestep, ITERATOR cell_detections_begin, ITERATOR cell_detections_end)
  {
    std::cout << "do not use this function anymore\n";
    assert(false);
    INDEX size = 0;
    for(auto it=cell_detections_begin; it!=cell_detections_end; ++it) {
      if(detection_factors_[timestep][ *it ] != nullptr) {
        ++size; 
      }
    }
    if(size == 0) { return nullptr; }

    auto* e = new AT_MOST_ONE_CELL_FACTOR_CONTAINER(size);
    lp.AddFactor(e);
    //std::cout << "Exclusion for time " << timestep << ": ";
    INDEX msg_idx = 0;
    for(auto it=cell_detections_begin; it!=cell_detections_end; ++it) {
      //std::cout << cell_detections[i] << ", ";
      if(detection_factors_[timestep][ *it ] != nullptr) {
        auto* f = detection_factors_[ timestep ][ *it ];
        auto* m = new AT_MOST_ONE_CELL_MESSAGE_CONTAINER(msg_idx, f, e);
        ++msg_idx;
        lp.AddMessage(m); 
      }
    }
    //std::cout << std::endl;
    return e; 
  }

  template<typename LP_TYPE, typename ITERATOR>
  AT_MOST_ONE_CELL_FACTOR_CONTAINER* add_exclusion_constraint(LP_TYPE& lp, ITERATOR begin, ITERATOR end) // iterator points to std::array<INDEX,2>
  {
    //std::sort(begin,end, [](auto a, auto b) { if(a[0] != b[0]) return a[0] < b[0]; else return a[1] < b[1]; } );
    assert(std::distance(begin, end) > 1);
    INDEX size = 0;
    std::array<INDEX,2> min_detection_factor = {std::numeric_limits<INDEX>::max(), std::numeric_limits<INDEX>::max()};
    std::array<INDEX,2> max_detection_factor = {0,0};

    for(auto it=begin; it!=end; ++it) {
      const INDEX timestep = (*it)[0];
      const INDEX hypothesis_id = (*it)[1];
      if(timestep < min_detection_factor[0] || (timestep == min_detection_factor[0] && hypothesis_id <= min_detection_factor[1])) {
        min_detection_factor[0] = timestep;
        min_detection_factor[1] = hypothesis_id;
      }
      if(timestep > max_detection_factor[0] || (timestep == max_detection_factor[0] && hypothesis_id >= max_detection_factor[1])) {
        max_detection_factor[0] = timestep;
        max_detection_factor[1] = hypothesis_id;
      }

      assert(detection_factors_[ timestep ][ hypothesis_id ] != nullptr); // for now, this need not hold true
      if(detection_factors_[ timestep ][ hypothesis_id ] != nullptr) {
        ++size; 
      }
    }

    auto* f_first = detection_factors_[ min_detection_factor[0] ][ min_detection_factor[1] ];
    auto* f_last  = detection_factors_[ max_detection_factor[0] ][ max_detection_factor[1] ];
    assert(f_first != nullptr);
    assert(f_last != nullptr);
    assert(f_first != f_last);

    assert(size > 1); // for now, this need not hold true
    if(size <= 1) { return nullptr; }

    auto* e = new AT_MOST_ONE_CELL_FACTOR_CONTAINER(size);
    lp.AddFactor(e);
    INDEX msg_idx = 0;
    for(auto it=begin; it!=end; ++it) {
      const INDEX timestep = (*it)[0];
      const INDEX hypothesis_id = (*it)[1];
      if(detection_factors_[timestep][hypothesis_id] != nullptr) {

        // only valid for sorted conflict set
        auto* f = detection_factors_[ timestep ][ hypothesis_id ];
        auto* m = new AT_MOST_ONE_CELL_MESSAGE_CONTAINER(msg_idx, f, e);

        ++msg_idx;
        lp.AddMessage(m); 
        if(f != f_last) {
          lp.ForwardPassFactorRelation(f,e);
        } else {
          lp.ForwardPassFactorRelation(e,f);
        }
        if(f != f_first) {
          lp.BackwardPassFactorRelation(f,e);
        } else {
          lp.BackwardPassFactorRelation(e,f);
        }
      }
    }

    //std::cout << std::endl;
    return e; 
  }

  template<typename LP_TYPE>
  void order_factors(LP_TYPE& lp) const
  {
    for(INDEX t=0; t<detection_factors_.size(); ++t) {
      assert(detection_factors_[t].size() > 0);
      for(INDEX i=0; i<detection_factors_[t].size(); ++i) {
        assert(detection_factors_[t][i] != nullptr);
      }
      if(t > 0) {
        lp.AddFactorRelation(detection_factors_[t-1].back(), detection_factors_[t][0]);
      }
      for(INDEX i=1; i<detection_factors_[t].size(); ++i) {
        lp.AddFactorRelation(detection_factors_[t][i-1], detection_factors_[t][i]);
      }
    }
  }

  // add violated exclusion constraints
  INDEX Tighten(const INDEX max_constraints_to_add)
  {
    assert(false); // not currently in use
    std::vector<std::tuple<INDEX,REAL>> exclusion_candidates; // offset into exclusions_ and guaranteed dual increase 

    // to do: parallelize
    for(INDEX i=0; i<exclusions_.size();) {
      auto item_begin = exclusions_.begin() + i;
      auto item_end = item_begin;
      while((*item_end) != exclusion_item_delimiter) {
        ++item_end; 
      }
      // check how large guaranteed dual increase is when we would add exclusion constraint
      REAL sum_detection_costs = 0.0;
      REAL smallest_detection_cost = std::numeric_limits<REAL>::infinity();
      for(auto it=item_begin; it!=item_end; ++it) {
        const REAL lb = detection_factors_[it->timestep][it->hypothesis_id]->LowerBound();
        assert(lb <= 0.0);
        sum_detection_costs += lb;
        smallest_detection_cost = std::min(lb, smallest_detection_cost);
        if(sum_detection_costs + eps < smallest_detection_cost) {
          exclusion_candidates.push_back( std::make_tuple(i, smallest_detection_cost - sum_detection_costs) );
        } 
      }
      i += 1 + std::distance(item_begin, item_end);
    }

    // possibly not necessary anymore
    std::sort(exclusion_candidates.begin(), exclusion_candidates.end(), [](auto a, auto b) { return std::get<1>(a) > std::get<1>(b); });

    if(exclusion_candidates.size() > max_constraints_to_add) {
      exclusion_candidates.resize(max_constraints_to_add);
    }
    for(INDEX i=0; i<exclusion_candidates.size(); ++i) {
      const INDEX idx = std::get<0>(exclusion_candidates[i]);
      auto item_begin = exclusions_.begin() + i;
      auto item_end = item_begin;
      while(*(item_end) != exclusion_item_delimiter) {
        ++item_end; 
      }
      add_exclusion_constraint(*lp_, item_begin, item_end); 
      i += 1 + std::distance(item_begin, item_end);
    }

    constexpr exclusion_item removal_mark = {std::numeric_limits<INDEX>::max(),0};
    // remove exclusions that were added
    for(INDEX i=0; i<exclusion_candidates.size(); ++i) {
      const INDEX idx = std::get<0>(exclusion_candidates[i]);
      auto it = exclusions_.begin() + idx;
      for(; (*it)!=exclusion_item_delimiter; ++it) {
        *it = removal_mark;
      }
      *it = removal_mark; 
    }
    auto it = std::remove(exclusions_.begin(), exclusions_.end(), removal_mark);
    exclusions_.resize(std::distance(exclusions_.begin(), it));

    if(verbosity >= 2) {
      std::cout << "added " << exclusion_candidates.size() << " exclusion factors, " << exclusions_.size() << " exclusions remain\n";
    }

    return exclusion_candidates.size();
  }

protected:
  using detection_factors_storage = std::vector<std::vector<DETECTION_FACTOR_CONTAINER*>>;
  detection_factors_storage detection_factors_;

  cell_tracking_transition_count tc_;

  using exclusion_factor_storage = std::vector<exclusion_item>;
  exclusion_factor_storage exclusions_; // hold all exclusions in a single array. delimiter is std::numeric_limits<INDEX>::max(). First entry in segment is timestep, followed by hypothesis ids

  LP* lp_;
};


template<typename BASIC_CELL_TRACKING_CONSTRUCTOR, typename TRANSITION_MESSAGE_CONTAINER>
class cell_tracking_constructor : public BASIC_CELL_TRACKING_CONSTRUCTOR {
public:
  using BASIC_CELL_TRACKING_CONSTRUCTOR::BASIC_CELL_TRACKING_CONSTRUCTOR;
  using transition_message_container = TRANSITION_MESSAGE_CONTAINER;

  //void set_number_of_timesteps(const INDEX t)
  //{
  //  assert(detection_factors_.size() == 0);
  //  detection_factors_.resize(t);
  //}

  template<typename LP_TYPE>
  void add_cell_transition(LP_TYPE& lp, const INDEX timestep_prev, const INDEX prev_cell, const INDEX timestep_next, const INDEX next_cell, const REAL cost) 
  {
    assert(timestep_prev + 1 == timestep_next);

    auto* out_cell_factor = this->detection_factors_[timestep_prev][prev_cell];
    const INDEX outgoing_edge_index  = this->tc_.next_outgoing_transition_edge(timestep_prev, prev_cell);//current_transition_no[timestep_prev][prev_cell][1];
    out_cell_factor->GetFactor()->set_outgoing_transition_cost(outgoing_edge_index, 0.5*cost);
    //tc.current_transition_no[timestep_prev][prev_cell][1]++;

    auto* in_cell_factor = this->detection_factors_[timestep_next][next_cell];
    const INDEX incoming_edge_index = this->tc_.next_incoming_transition_edge(timestep_next, next_cell);//current_transition_no[timestep_next][next_cell][0];
    //tc.current_transition_no[timestep_next][next_cell][0]++;
    in_cell_factor->GetFactor()->set_incoming_transition_cost(incoming_edge_index, 0.5*cost);
    auto* m = new TRANSITION_MESSAGE_CONTAINER(out_cell_factor, in_cell_factor, false, outgoing_edge_index, incoming_edge_index);
    lp.AddMessage(m);

    //std::cout << "MA: " << timestep << " " << prev_cell << ", " << next_cell << " " << cost << std::endl;
  }

  template<typename LP_TYPE>
  void add_cell_division(LP_TYPE& lp, const INDEX timestep_prev, const INDEX prev_cell, const INDEX timestep_next_1, const INDEX next_cell_1, const INDEX timestep_next_2, const INDEX next_cell_2, const REAL cost) 
  {
    auto* out_cell_factor = this->detection_factors_[timestep_prev][prev_cell];
    const INDEX outgoing_edge_index  = this->tc_.next_outgoing_division_edge(timestep_prev, prev_cell);//tc.current_transition_no[timestep_prev][prev_cell][1] + tc.current_division_no[timestep_prev][prev_cell][1];
    out_cell_factor->GetFactor()->set_outgoing_division_cost(outgoing_edge_index, 1.0/3.0*cost);
    //tc.current_division_no[timestep_prev][prev_cell][1]++;

    auto* in_cell_factor_1 = this->detection_factors_[timestep_next_1][next_cell_1];
    const INDEX incoming_edge_index_1 = this->tc_.next_incoming_division_edge(timestep_next_1, next_cell_1);//tc.current_transition_no[timestep_next_1][next_cell_1][0] + tc.current_division_no[timestep_next_1][next_cell_1][0];
    in_cell_factor_1->GetFactor()->set_incoming_division_cost(incoming_edge_index_1, 1.0/3.0*cost);
    //tc.current_division_no[timestep_next_1][next_cell_1][0]++;
    
    auto* in_cell_factor_2 = this->detection_factors_[timestep_next_2][next_cell_2];
    const INDEX incoming_edge_index_2 = this->tc_.next_incoming_division_edge(timestep_next_2, next_cell_2);//tc.current_transition_no[timestep_next_2][next_cell_2][0] + tc.current_division_no[timestep_next_2][next_cell_2][0];
    in_cell_factor_2->GetFactor()->set_incoming_division_cost(incoming_edge_index_2, 1.0/3.0*cost);
    //tc.current_division_no[timestep_next_2][next_cell_2][0]++;
    
    auto* m1 = new TRANSITION_MESSAGE_CONTAINER(out_cell_factor, in_cell_factor_1, true, outgoing_edge_index, incoming_edge_index_1);
    lp.AddMessage(m1);
    
    auto* m2 = new TRANSITION_MESSAGE_CONTAINER(out_cell_factor, in_cell_factor_2, true, outgoing_edge_index, incoming_edge_index_2);
    lp.AddMessage(m2);

    //std::cout << "DA: " << timestep << " " << prev_cell << ", " << next_cell_1 << " " << next_cell_2 << " " << cost << std::endl;
  }
};

template<typename BASIC_CELL_TRACKING_CONSTRUCTOR, typename MAPPING_EDGE_FACTOR_CONTAINER, typename INCOMING_EDGE_MESSAGE_CONTAINER, typename OUTGOING_EDGE_MESSAGE_CONTAINER>
class cell_tracking_constructor_duplicate_edges : public BASIC_CELL_TRACKING_CONSTRUCTOR {
public:
  using CONSTRUCTOR = cell_tracking_constructor_duplicate_edges<BASIC_CELL_TRACKING_CONSTRUCTOR, MAPPING_EDGE_FACTOR_CONTAINER, INCOMING_EDGE_MESSAGE_CONTAINER, OUTGOING_EDGE_MESSAGE_CONTAINER>;

  using BASIC_CELL_TRACKING_CONSTRUCTOR::BASIC_CELL_TRACKING_CONSTRUCTOR;

  template<typename LP_TYPE>
  void add_cell_transition(LP_TYPE& lp, const INDEX timestep_prev, const INDEX prev_cell, const INDEX timestep_next, const INDEX next_cell, const REAL cost) 
  {
    assert(timestep_prev + 1 == timestep_next);

    auto* out_cell_factor = this->detection_factors_[timestep_prev][prev_cell];
    const INDEX outgoing_edge_index  = this->tc_.next_incoming_transition_edge(timestep_prev, prev_cell);//tc.current_transition_no[timestep_prev][prev_cell][1];
    //assert( out_cell_factor->GetFactor()->outgoing[outgoing_edge_index] == 0.0 );
    out_cell_factor->GetFactor()->set_outgoing_transition_cost(outgoing_edge_index, 0.5*cost);
    //out_cell_factor->GetFactor()->outgoing[outgoing_edge_index] = cost;
    //tc.current_transition_no[timestep_prev][prev_cell][1]++;

    auto* in_cell_factor = this->detection_factors_[timestep_next][next_cell];
    const INDEX incoming_edge_index = this->tc_.next_incoming_transition_edge(timestep_next, next_cell);//tc.current_transition_no[timestep_next][next_cell][0];
    //assert( in_cell_factor->GetFactor()->set_incoming[incoming_edge_index] == 0.0 );
    //tc.current_transition_no[timestep_next][next_cell][0]++;
    in_cell_factor->GetFactor()->set_incoming_transition_cost(incoming_edge_index, 0.5*cost);

    auto* f = new MAPPING_EDGE_FACTOR_CONTAINER(0.0);
    lp.AddFactor(f);

    auto* m_outgoing = new OUTGOING_EDGE_MESSAGE_CONTAINER(f, out_cell_factor, outgoing_edge_index, false);
    lp.AddMessage(m_outgoing);
    auto* m_incoming = new INCOMING_EDGE_MESSAGE_CONTAINER(f, in_cell_factor, incoming_edge_index);
    lp.AddMessage(m_incoming);
    lp.AddFactorRelation(out_cell_factor, f);
    lp.AddFactorRelation(f, in_cell_factor);

    //std::cout << "MA: " << timestep << " " << prev_cell << ", " << next_cell << " " << cost << std::endl;
  }

  template<typename LP_TYPE>
  void add_cell_division(LP_TYPE& lp, const INDEX timestep_prev, const INDEX prev_cell, const INDEX timestep_next_1, const INDEX next_cell_1, const INDEX timestep_next_2, const INDEX next_cell_2, const REAL cost) 
  {
    // this doees not work for cell tracking without minimal distance. Then the edge indices need to be summed up with the transition edge indices
    auto* out_cell_factor = this->detection_factors_[timestep_prev][prev_cell];
    const INDEX outgoing_edge_index  = this->tc_.next_outgoing_division_edge(timestep_prev, prev_cell);//tc.current_transition_no[timestep_prev][prev_cell][1] + tc.current_division_no[timestep_prev][prev_cell][1];
    out_cell_factor->GetFactor()->set_outgoing_division_cost(outgoing_edge_index, 1.0/3.0*cost);
    //tc.current_division_no[timestep_prev][prev_cell][1]++;

    auto* in_cell_factor_1 = this->detection_factors_[timestep_next_1][next_cell_1];
    const INDEX incoming_edge_index_1 = this->tc_.next_incoming_division_edge(timestep_next_1, next_cell_1);//tc.current_transition_no[timestep_next_1][next_cell_1][0] + tc.current_division_no[timestep_next_1][next_cell_1][0];
    in_cell_factor_1->GetFactor()->set_incoming_division_cost(incoming_edge_index_1, 1.0/3.0*cost);
    //tc.current_division_no[timestep_next_1][next_cell_1][0]++;
    
    auto* in_cell_factor_2 = this->detection_factors_[timestep_next_2][next_cell_2];
    const INDEX incoming_edge_index_2 = this->tc_.next_incoming_division_edge(timestep_next_2, next_cell_2);//tc.current_transition_no[timestep_next_2][next_cell_2][0] + tc.current_division_no[timestep_next_2][next_cell_2][0];
    in_cell_factor_2->GetFactor()->set_incoming_division_cost(incoming_edge_index_2, 1.0/3.0*cost);
    //tc.current_division_no[timestep_next_2][next_cell_2][0]++;
    
    auto* f_1 = new MAPPING_EDGE_FACTOR_CONTAINER(0.0);
    lp.AddFactor(f_1);
    auto* f_2 = new MAPPING_EDGE_FACTOR_CONTAINER(0.0);
    lp.AddFactor(f_2);

    auto* m_outgoing_1 = new OUTGOING_EDGE_MESSAGE_CONTAINER(f_1, out_cell_factor, outgoing_edge_index, true);
    auto* m_outgoing_2 = new OUTGOING_EDGE_MESSAGE_CONTAINER(f_2, out_cell_factor, outgoing_edge_index, true);

    auto* m_incoming_1 = new INCOMING_EDGE_MESSAGE_CONTAINER(f_1, in_cell_factor_1, incoming_edge_index_1);
    auto* m_incoming_2 = new INCOMING_EDGE_MESSAGE_CONTAINER(f_2, in_cell_factor_2, incoming_edge_index_2);

    lp.AddMessage(m_outgoing_1);
    lp.AddMessage(m_outgoing_2);
    lp.AddMessage(m_incoming_1);
    lp.AddMessage(m_incoming_2);

    lp.AddFactorRelation(out_cell_factor, f_1);
    lp.AddFactorRelation(out_cell_factor, f_2);
    lp.AddFactorRelation(f_1, in_cell_factor_1);
    lp.AddFactorRelation(f_2, in_cell_factor_2);

    //std::cout << "DA: " << timestep << " " << prev_cell << ", " << next_cell_1 << " " << next_cell_2 << " " << cost << std::endl;
  }

};

template<typename CELL_TRACKING_CONSTRUCTOR>
class cell_tracking_with_division_distance_constructor : public CELL_TRACKING_CONSTRUCTOR {
public:
  using detection_factor_container = typename CELL_TRACKING_CONSTRUCTOR::detection_factor_container;
  using transition_message_container = typename CELL_TRACKING_CONSTRUCTOR::transition_message_container;
  using CELL_TRACKING_CONSTRUCTOR::CELL_TRACKING_CONSTRUCTOR;

  void set_division_distance(const INDEX d) 
  {
    division_distance_ = d;
  }

  template<typename LP_TYPE>
  detection_factor_container* 
  add_detection_hypothesis(
      LP_TYPE& lp, 
      const INDEX timestep, const INDEX hypothesis_id, 
      const REAL detection_cost, const REAL appearance_cost, const REAL disappearance_cost, 
      const INDEX no_incoming_transition_edges, const INDEX no_incoming_division_edges, 
      const INDEX no_outgoing_transition_edges, const INDEX no_outgoing_division_edges
      )
  { 
    assert(division_distance_ >= 2);
    assert(timestep < this->detection_factors_.size());
    if(hypothesis_id >= this->detection_factors_[timestep].size()) {
      this->detection_factors_[timestep].resize(hypothesis_id+1, nullptr);
    }
    assert(this->detection_factors_[timestep][hypothesis_id] == nullptr);

    auto* f = new detection_factor_container(no_incoming_transition_edges, no_incoming_division_edges, no_outgoing_transition_edges, no_outgoing_division_edges, detection_cost, appearance_cost, disappearance_cost, division_distance_);
    lp.AddFactor(f);
    this->detection_factors_[timestep][hypothesis_id] = f;
    if(hypothesis_id > 0) {
      assert( this->detection_factors_[timestep][hypothesis_id-1] != nullptr); // need not be generally true, but then factor relation must be done more robust.
    }
    //std::cout << "H: " << timestep << ", " << hypothesis_id <<  "," << no_incoming_transition_edges << "," << no_incoming_division_edges << "," << no_outgoing_transition_edges << ", " << no_outgoing_division_edges << "," << detection_cost << ", " << appearance_cost << ", " << disappearance_cost << std::endl;

    this->tc_.add_detection_hypothesis(timestep, hypothesis_id, no_incoming_transition_edges, no_incoming_division_edges, no_outgoing_transition_edges, no_outgoing_division_edges);
    return f; 
  }
  
  template<typename LP_TYPE>
  void add_cell_division(LP_TYPE& lp, const INDEX timestep_prev, const INDEX prev_cell, const INDEX timestep_next_1, const INDEX next_cell_1, const INDEX timestep_next_2, const INDEX next_cell_2, const REAL cost) 
  {
    auto* out_cell_factor = this->detection_factors_[timestep_prev][prev_cell];
    const INDEX outgoing_edge_index  = this->tc_.next_outgoing_division_edge(timestep_prev, prev_cell);//tc.current_division_no[timestep_prev][prev_cell][1];
    out_cell_factor->GetFactor()->set_outgoing_division_cost(outgoing_edge_index, 1.0/3.0*cost);
    //tc.current_division_no[timestep_prev][prev_cell][1]++;

    auto* in_cell_factor_1 = this->detection_factors_[timestep_next_1][next_cell_1];
    const INDEX incoming_edge_index_1 = this->tc_.next_incoming_division_edge(timestep_next_1, next_cell_1);//tc.current_division_no[timestep_next_1][next_cell_1][0];
    in_cell_factor_1->GetFactor()->set_incoming_division_cost(incoming_edge_index_1, 1.0/3.0*cost);
    //tc.current_division_no[timestep_next_1][next_cell_1][0]++;
    
    auto* in_cell_factor_2 = this->detection_factors_[timestep_next_2][next_cell_2];
    const INDEX incoming_edge_index_2 = this->tc_.next_incoming_division_edge(timestep_next_2, next_cell_2);//tc.current_division_no[timestep_next_2][next_cell_2][0];
    in_cell_factor_2->GetFactor()->set_incoming_division_cost(incoming_edge_index_2, 1.0/3.0*cost);
    //tc.current_division_no[timestep_next_2][next_cell_2][0]++;
    
    auto* m1 = new transition_message_container(out_cell_factor, in_cell_factor_1, true, outgoing_edge_index, incoming_edge_index_1);
    lp.AddMessage(m1);
    
    auto* m2 = new transition_message_container(out_cell_factor, in_cell_factor_2, true, outgoing_edge_index, incoming_edge_index_2);
    lp.AddMessage(m2);

    //std::cout << "DA: " << timestep << " " << prev_cell << ", " << next_cell_1 << " " << next_cell_2 << " " << cost << std::endl;
  }

private:
  INDEX division_distance_ = 0;
};

// this constructor first builds an ordinary cell tracking problem, solves it, and afterwards converts it into a cell tracking with division distance problem and solves it when it is preprocessed
// version with no duplicate edges
template<
typename CELL_TRACKING_CONSTRUCTOR,
typename CELL_TRACKING_DIVISION_DISTANCE_CONSTRUCTOR >
class cell_tracking_division_distance_conversion_constructor
{
public:
  using detection_factor_container = typename CELL_TRACKING_CONSTRUCTOR::detection_factor_container;
  using at_most_one_cell_factor_container = typename CELL_TRACKING_CONSTRUCTOR::at_most_one_cell_factor_container;
  using exclusion_factor = typename CELL_TRACKING_CONSTRUCTOR::exclusion_factor;

  cell_tracking_division_distance_conversion_constructor()
    : cdc_()
  {}

  ~cell_tracking_division_distance_conversion_constructor()
  {
    if(ctc_dd_ != nullptr) {
      delete ctc_dd_;
    }
  }

  template<typename LP_TYPE>
  detection_factor_container* 
  add_detection_hypothesis(
      LP_TYPE& lp, 
      const INDEX timestep, const INDEX hypothesis_id, 
      const REAL detection_cost, const REAL appearance_cost, const REAL disappearance_cost, 
      const INDEX no_incoming_transition_edges, const INDEX no_incoming_division_edges, 
      const INDEX no_outgoing_transition_edges, const INDEX no_outgoing_division_edges
      )
  { 
    if(timestep >= no_transition_edges_.size()) {
      no_transition_edges_.resize(timestep+1);
    }
    if(hypothesis_id >= no_transition_edges_[timestep].size()) {
      no_transition_edges_[timestep].resize(hypothesis_id + 1, {0,0});
    }
    no_transition_edges_[timestep][hypothesis_id][0] = no_incoming_transition_edges;
    no_transition_edges_[timestep][hypothesis_id][1] = no_outgoing_transition_edges;

    cdc_.add_detection_hypothesis(
        lp, timestep, hypothesis_id, detection_cost, appearance_cost, disappearance_cost, 
        no_incoming_transition_edges, no_incoming_division_edges,
        no_outgoing_transition_edges, no_outgoing_division_edges
        );
  }

  template<typename LP_TYPE>
  void add_cell_transition(LP_TYPE& lp, const INDEX timestep_prev, const INDEX prev_cell, const INDEX timestep_next, const INDEX next_cell, const REAL cost) 
  {
    cdc_.add_cell_transition.push_back(lp, timestep_prev, prev_cell, timestep_next, next_cell, cost);
    transition_edges_.push_back({timestep_prev, prev_cell, timestep_next, next_cell});

  }

  template<typename LP_TYPE>
  void add_cell_division(LP_TYPE& lp, const INDEX timestep_prev, const INDEX prev_cell, const INDEX timestep_next_1, const INDEX next_cell_1, const INDEX timestep_next_2, const INDEX next_cell_2, const REAL cost) 
  {
    cdc_.add_cell_division(lp, timestep_prev, prev_cell, timestep_next_1, next_cell_1, timestep_next_2, next_cell_2, cost);
    assert(timestep_prev + 1 == timestep_next_1);
    assert(timestep_prev + 1 == timestep_next_2);
    division_edges_.push_back({timestep_prev, prev_cell, timestep_next_1, next_cell_1, timestep_next_2, next_cell_2});
  }

  template<typename ITERATOR>
  void register_exclusion_constraint(ITERATOR cell_detections_begin, ITERATOR cell_detections_end)
  {
    assert(std::distance(cell_detections_begin, cell_detections_end) > 1);
    for(auto it=cell_detections_begin; it!=cell_detections_end; ++it) {
      exclusions_.push_back({(*it)[0], (*it)[1]});
    }
    auto begin = exclusions_.end() - std::distance(cell_detections_begin, cell_detections_end);
    auto end = exclusions_.end();
    std::sort(begin, end);
    exclusions_.push_back(exclusion_item_delimiter);
  }

  template<typename LP_TYPE, typename ITERATOR>
  at_most_one_cell_factor_container* add_exclusion_constraint(LP_TYPE& lp, ITERATOR begin, ITERATOR end) // iterator points to std::array<INDEX,2>
  {
    auto* f = cdc_.add_exclusion_constraint(lp, begin, end);
    exclusion_factors_.push_back(f);
    register_exclusion_constraint(begin, end);
  }

  template<typename LP>
  void convert(LP* lp)
  {
    auto* ctc_dd = new CELL_TRACKING_DIVISION_DISTANCE_CONSTRUCTOR();
    // read all factors from ordinary cell tracking problem, and add equivalent ones to problem with division distance.
    // cost of new factors is reparametrized cost of old ones

    // detection factors
    for(INDEX t=0; t<cdc_.detection_factors_.size(); ++t) {
      for(INDEX i=0; i<cdc_.detection_factors_[i].size(); ++i) {
        auto* fp = cdc_.detection_factors_[i][t];
        auto& f = *(fp->GetFactor());
        const REAL detection_cost = f.detection_cost();
        const REAL appearance_cost = f.appearance_cost();
        const REAL disappearance_cost = f.disappearance_cost();
        const INDEX no_incoming_transition = no_transition_edges_[i][t][0];
        const INDEX no_outgoing_transition = no_transition_edges_[i][t][1];
        const INDEX no_incoming_division = f.no_incoming_edges() - 1 - no_incoming_transition;
        const INDEX no_outgoing_division = f.no_outgoing_edges() - 1 - no_outgoing_transition;

        auto* fp_dd = ctc_dd_->add_detection_hypothesis( 
            lp, t, 
            detection_cost, appearance_cost, disappearance_cost,
            no_incoming_transition, no_incoming_division,
            no_outgoing_transition, no_outgoing_division
            );
        auto& f_dd = *(fp_dd->GetFactor());

        for(INDEX incoming_edge = 0; incoming_edge<no_incoming_transition; ++incoming_edge) {
          const REAL cost = f.incoming[incoming_edge];
          f_dd.set_incoming_transition_cost(incoming_edge, cost); 
        }
        for(INDEX incoming_edge = 0; incoming_edge<no_incoming_division; ++incoming_edge) {
          const REAL cost = f.incoming[ incoming_edge + no_transition_edges_[t][i] ] ;
          f_dd.set_incoming_division_cost(incoming_edge, cost); 
        }

        for(INDEX outgoing_edge = 0; outgoing_edge<no_outgoing_transition; ++outgoing_edge) {
          const REAL cost = f.outgoing[outgoing_edge];
          f_dd.set_outgoing_transition_cost(outgoing_edge, cost); 
        }
        for(INDEX outgoing_edge = 0; outgoing_edge<no_outgoing_division; ++outgoing_edge) {
          const REAL cost = f.outgoing[ outgoing_edge + no_transition_edges_[t][i] ] ;
          f_dd.set_outgoing_division_cost(outgoing_edge, cost); 
        } 
      } 
    }

    // link detections via edges
    /*
    {
      std::vector<std::vector<std::array<INDEX,2>>> edge_counter(cdc_.detection_factors_.size());
      for(INDEX i=0; i<edge_counter.size(); ++i) {
        edge_counter[i].resize( cdc_.detection_factors_[i].size(), {0,0} );
      }
      for(const auto t : transition_edges_) {
        auto* f_prev = ctc_dd->detection_factors_[t[0]][t[1]];
        auto* f_next = ctc_dd_->detection_factors_[t[2]][t[3]];
        const INDEX outgoing_edge_index = edge_counter[t[0]][t[1]][1]++;
        const INDEX incoming_edge_index = edge_counter[t[2]][t[3]][0]++;

        auto* m = new transition_message_container(f_prev, f_next, false, outgoing_edge_index, incoming_edge_index);
        lp.AddMessage(m);
      }
      for(const auto t : division_edges) {
        auto* f_prev = ctc_dd_->detection_factors_[t[0]][t[1]];
        auto* f_next_1 = ctc_dd->detection_factors_[t[2]][t[3]];
        auto* f_next_2 = ctc_dd->detection_factors_[t[4]][t[5]];
        const INDEX outgoing_edge_index = edge_counter[t[0]][t[1]][1]++;
        const INDEX incoming_edge_index_1 = edge_counter[t[2]][t[3]][0]++;
        const INDEX incoming_edge_index_2 = edge_counter[t[4]][t[5]][0]++;

        auto* m1 = new transition_message_container(f_prev, f_next_1, true, outgoing_edge_index, incoming_edge_index_1);
        lp.AddMessage(m1);

        auto* m2 = new transition_message_container(f_prev, f_next_2, true, outgoing_edge_index, incoming_edge_index_2);
        lp.AddMessage(m2); 
      }
    }
    */

    // exclusion factors
    /*
    {
      auto ef_it = exclusion_factors_.begin();
      for(INDEX i=0; i<cdc_.exclusions_.size();) {
        auto item_begin = cdc_.exclusions_.begin() + i;
        auto item_end = item_begin;
        while((*item_end) != exclusion_item_delimiter) {
          ++item_end; 
        }
        auto* f = ctc_dd_->add_exclusion_constraint(*lp_, item_begin, item_end); 
        i += 1 + std::distance(item_begin, item_end);

        auto* f_old = *ef_it;
        ++ef_it;
        for(INDEX j=0; j<f_old.size(); ++j) {
          f[j] = f_old[j];
        } 
      }
      assert(ef_it == exclusions_factors_.end());
    }
    */
    
    // free up space taken up by conversion information
    std::swap(no_transition_edges_,decltype(no_transition_edges_){});
    std::swap(transition_edges_,decltype(transition_edges_){});
    std::swap(division_edges_,decltype(division_edges_){});
    std::swap(exclusion_factors_,decltype(exclusion_factors_){});
    std::swap(exclusions_,decltype(exclusions_){});

  }
protected:
  // store the number of incoming/outgoing transition edges here
  std::vector<std::vector<std::array<INDEX,2>>> no_transition_edges_; 
  std::vector<std::array<INDEX,4>> transition_edges_;
  std::vector<std::array<INDEX,6>> division_edges_;

  std::vector<exclusion_factor*> exclusion_factors_;
  std::vector<exclusion_item> exclusions_;

  CELL_TRACKING_CONSTRUCTOR cdc_;
  CELL_TRACKING_DIVISION_DISTANCE_CONSTRUCTOR* ctc_dd_ = nullptr;
};

// we assume that cells are ordered by heigth and lower cells cannot exit if higher ones do not do so.
template<typename CELL_TRACKING_CONSTRUCTOR,
  typename EXIT_CONSTRAINT_FACTOR,
  typename EXIT_CONSTRAINT_LOWER_MESSAGE,
  typename EXIT_CONSTRAINT_UPPER_MESSAGE
  >
class cell_tracking_mother_machine_constructor : public CELL_TRACKING_CONSTRUCTOR {
public:
  using CELL_TRACKING_CONSTRUCTOR::CELL_TRACKING_CONSTRUCTOR;
  template<typename LP_TYPE>
  void add_exit_constraint(LP_TYPE& lp, const INDEX timestep, const INDEX lower_cell_detection, const INDEX upper_cell_detection)
  {
    if(this->detection_factors_[timestep][lower_cell_detection] == nullptr || this->detection_factors_[timestep][upper_cell_detection] == nullptr) { 
      return;
    }
    //std::cout << "t = " << timestep << ", lower cell = " << lower_cell_detection << ", upper cell = " << upper_cell_detection;
    if(this->detection_factors_[timestep][upper_cell_detection]->GetFactor()->no_outgoing_edges() > 1) { // only if there are outgoing edges which are not exit ones do we need exit constraints
      //std::cout << " add factor";
      auto* e = new EXIT_CONSTRAINT_FACTOR();
      lp.AddFactor(e);
      auto* m_lower = new EXIT_CONSTRAINT_LOWER_MESSAGE(this->detection_factors_[timestep][lower_cell_detection], e);
      lp.AddMessage(m_lower);
      auto* m_upper = new EXIT_CONSTRAINT_UPPER_MESSAGE(this->detection_factors_[timestep][upper_cell_detection], e);
      lp.AddMessage(m_upper); 
    }
    //std::cout << "\n";
  } 
};

namespace cell_tracking_parser_mother_machine {

   // import basic parsers
   using Parsing::opt_whitespace;
   using Parsing::mand_whitespace;
   using Parsing::opt_invisible;
   using Parsing::mand_invisible;
   using Parsing::positive_integer;
   using Parsing::real_number;

   struct my_eolf : pegtl::sor<pegtl::eol, pegtl::eof> {}; // do zrobienia: revert to usual pegtl::eolf
   struct comment_line : pegtl::seq< opt_whitespace, pegtl::sor< pegtl::seq< pegtl::string<'#'>, pegtl::until<my_eolf> >, my_eolf> > {};
   
   // t = number
   struct timestep : pegtl::seq< positive_integer > {};
   struct timestep_line : pegtl::seq< pegtl::string<'t'>, opt_whitespace, pegtl::string<'='>, opt_whitespace, timestep, opt_whitespace, my_eolf > {};

   // H hypothesis_number hypothesis_id {detection cost} {exit cost} (upper pos, lower pos) // upper < lower, inverted!
   struct cell_detection_hypothesis : pegtl::seq< positive_integer, opt_whitespace, positive_integer, opt_whitespace, real_number, opt_whitespace, real_number, opt_whitespace, pegtl::string<'('>, opt_whitespace, positive_integer, opt_whitespace, pegtl::string<','>, opt_whitespace, positive_integer, opt_whitespace, pegtl::string<')'> > {};
   struct cell_detection_hypothesis_line : pegtl::seq< pegtl::string<'H'>, opt_whitespace, cell_detection_hypothesis, opt_whitespace, my_eolf > {};
   
   // EC hyp_no1 + hyp_no2 ... <= 1
   struct exclusion : pegtl::seq< positive_integer, pegtl::star< pegtl::seq< opt_whitespace, pegtl::string<'+'>, opt_whitespace, positive_integer> > > {};
   struct exclusion_line : pegtl::seq< pegtl::string<'E','C'>, opt_whitespace, exclusion, opt_whitespace, pegtl::string<'<','='>, opt_whitespace, pegtl::string<'1'>, opt_whitespace, my_eolf> {};
   // timestep 1, first hypothesis id, timestep 1, second hypothesis id, cost
   struct mapping : pegtl::seq< positive_integer, mand_whitespace, positive_integer, mand_whitespace, positive_integer, mand_whitespace, positive_integer, mand_whitespace, real_number > {}; 
   struct mapping_line : pegtl::seq< pegtl::string<'M','A'>, opt_whitespace, mapping, opt_whitespace, my_eolf> {};
   // timestep 1, left hypothesis id, timestep 2, two right hypothesis ids, cost
   struct division : pegtl::seq< positive_integer, mand_whitespace, positive_integer, mand_whitespace, positive_integer, mand_whitespace, positive_integer, mand_whitespace, positive_integer, mand_whitespace, real_number > {}; 
   struct division_line : pegtl::seq< pegtl::string<'D','A'>, opt_whitespace, division, opt_whitespace, my_eolf> {};

   struct grammar : pegtl::seq< pegtl::star<pegtl::sor<
                    comment_line,
                    timestep_line, 
                    cell_detection_hypothesis_line,
                    exclusion_line,
                    mapping_line,
                    division_line
                    >>> {};

   struct input {
     // we must go through all cell mappings (cell transitions) to count the number of outgoing and incoming messages. Only then can we allocate cell detection hypothesis, as they need to know these numbers.
     std::vector<std::vector<std::tuple<REAL,REAL,INDEX,INDEX,REAL,REAL>>> cell_detection_stat_; // detection cost, exit cost,  # incoming edges, # outgoing edges, ( lower_boundary, upper_boundary ) <= last two only for mother machine
     std::vector<std::tuple<INDEX,INDEX,INDEX,REAL>> transitions_; // timestep, outgoing, incoming cell, cost
     std::vector<std::tuple<INDEX,INDEX,INDEX,INDEX,REAL>> divisions_; // timestep, outgoing, incoming1, incoming2, cost
     std::vector<std::vector<std::vector<INDEX>>> exclusions_; // timestep, {hyp_1, ..., hyp_n}
     bool eof = false;
   };

   template< typename Rule >
     struct action
     : pegtl::nothing< Rule > {};

   template<> struct action< pegtl::eof > {
     template<typename INPUT>
     static void apply(const INPUT & in, input& i)
     {
       //std::cout << "kwaskwas\n";
       i.eof = true;
     }
   };

   template<> struct action< timestep > {
     template<typename INPUT>
     static void apply(const INPUT & in, input& i)
     {
       const INDEX timestep = std::stoul(in.string());
       assert(i.cell_detection_stat_.size() == timestep);
       i.cell_detection_stat_.push_back({});
       i.exclusions_.push_back({});
     }
   };

   template<> struct action< cell_detection_hypothesis > {
     template<typename INPUT>
     static void apply(const INPUT & in, input& i)
     {
       std::stringstream s(in.string());
       INDEX number; s >> number;
       assert(number == i.cell_detection_stat_.back().size());
       INDEX hypothesis_id; s >> hypothesis_id; // this number is not used
       REAL detection_cost; s >> detection_cost;
       REAL exit_cost; s >> exit_cost;
       char opening_bracket; s >> opening_bracket; assert(opening_bracket == '(');
       REAL upper_boundary; s >> upper_boundary;
       char comma; s >> comma; assert(comma == ',');
       REAL lower_boundary; s >> lower_boundary;
       char closing_bracket; s >> closing_bracket; assert(closing_bracket == ')');
       assert(upper_boundary < lower_boundary );

       i.cell_detection_stat_.back().push_back(std::make_tuple(detection_cost,exit_cost,0,0, upper_boundary, lower_boundary));
       //std::cout << "H: " << number << ", " << detection_cost << std::endl;
     }
   };

   template<> struct action< exclusion > {
     template<typename INPUT>
     static void apply(const INPUT & in, input& i)
     {
       std::stringstream s(in.string());
       std::vector<INDEX> idx;
       INDEX number; s >> number; idx.push_back(number);
       //std::cout << "E: " << number << ", ";
       char plus;
       while(s >> plus) {
         assert(plus == '+');
         INDEX number; s >> number; idx.push_back(number); 
         //std::cout << number << ", ";
       }
       i.exclusions_.back().push_back(idx);
       //std::cout << " <= 1; " << i.exclusions_.size() << "\n";
     }
   };

   template<> struct action< mapping > {
     template<typename INPUT>
     static void apply(const INPUT & in, input& i)
     {
       std::stringstream s(in.string());
       INDEX timestep_prev; s >> timestep_prev;
       INDEX cell_prev; s >> cell_prev;
       INDEX timestep_next; s >> timestep_next;
       INDEX cell_next; s >> cell_next;
       REAL cost; s >> cost;

       assert(timestep_prev+1 == timestep_next);
       assert(timestep_next < i.cell_detection_stat_.size());
       assert(cell_prev < i.cell_detection_stat_[timestep_prev].size());
       assert(cell_next < i.cell_detection_stat_[timestep_next].size());

       std::get<3>(i.cell_detection_stat_[timestep_prev][cell_prev])++;
       std::get<2>(i.cell_detection_stat_[timestep_next][cell_next])++;

       i.transitions_.push_back( std::make_tuple(timestep_prev, cell_prev, cell_next, cost));

       //std::cout << "mapping: t = " << timestep_prev << ", h1 = " << cell_prev << ", h2 = " << cell_next << ", cost = " << cost << "\n";
     }
   };

   template<> struct action< division > {
     template<typename INPUT>
     static void apply(const INPUT & in, input& i)
     {
       std::stringstream s(in.string());
       INDEX timestep_prev; s >> timestep_prev;
       INDEX cell_prev; s >> cell_prev;
       INDEX timestep_next; s >> timestep_next;
       INDEX cell_next_1; s >> cell_next_1;
       INDEX cell_next_2; s >> cell_next_2;
       REAL cost; s >> cost;

       assert(timestep_prev+1 == timestep_next);
       assert(cell_next_1 != cell_next_2);
       assert(timestep_next < i.cell_detection_stat_.size());
       assert(cell_prev < i.cell_detection_stat_[timestep_prev].size());
       assert(cell_next_1 < i.cell_detection_stat_[timestep_next].size());
       assert(cell_next_2 < i.cell_detection_stat_[timestep_next].size());

       std::get<3>(i.cell_detection_stat_[timestep_prev][cell_prev])++;
       std::get<2>(i.cell_detection_stat_[timestep_next][cell_next_1])++;
       std::get<2>(i.cell_detection_stat_[timestep_next][cell_next_2])++;

       i.divisions_.push_back( std::make_tuple(timestep_prev, cell_prev, cell_next_1, cell_next_2, cost));
       //std::cout << "division: t = " << timestep_prev << ", h1 = " << cell_prev << ", h1,1 = " << cell_next_1 << ", h2,2 = " << cell_next_2 << ", cost = " << cost << "\n";
     }
   };


   bool read_input(const std::string& filename, input& i)
   {
      if(verbosity >= 1) {
        std::cout << "parsing " << filename << "\n";
      }
      pegtl::file_parser problem(filename);
      const bool success = problem.parse< grammar, action >(i); 
      //assert(i.eof);
      return success;
   }

   template<typename SOLVER>
   void construct_tracking_problem(input& i, SOLVER& s)
   {
      auto& cell_tracking_constructor = s.template GetProblemConstructor<0>();
      auto& lp = s.GetLP();

      //std::cout << "exclusion constraints disabled\n";
      cell_tracking_constructor.set_number_of_timesteps( i.cell_detection_stat_.size() );
      for(INDEX t=0; t<i.cell_detection_stat_.size(); ++t) {
        for(INDEX n=0; n<i.cell_detection_stat_[t].size(); ++n) {
          const REAL detection_cost = std::get<0>(i.cell_detection_stat_[t][n]);
          const REAL exit_cost = std::get<1>(i.cell_detection_stat_[t][n]);
          const INDEX no_incoming_edges = std::get<2>(i.cell_detection_stat_[t][n]);
          const INDEX no_outgoing_edges = std::get<3>(i.cell_detection_stat_[t][n]);

          //assert(t != i.cell_detection_stat_.size()-1 || exit_cost == 0.0);
          cell_tracking_constructor.add_detection_hypothesis( lp, t, n, detection_cost, no_incoming_edges, no_outgoing_edges, 0.0, exit_cost );
        }
        for(const auto& detections : i.exclusions_[t]) {
          cell_tracking_constructor.add_exclusion_constraint(lp, t, detections ); 
        }
        /* check whether all cell hypotheses are covered by exclusion constraints (must be valid for mother machine)
        std::vector<INDEX> all_indices;
        for(const auto& detections : i.exclusions_[t]) {
          all_indices.insert(all_indices.end(), detections.begin(), detections.end());
        }
        std::sort(all_indices.begin(), all_indices.end());
        all_indices.erase( std::unique( all_indices.begin(), all_indices.end()), all_indices.end() );
        assert(all_indices.size() == i.cell_detection_stat_[t].size());
        */
      }

      //auto tc = cell_tracking_constructor.init_transition_counter();
      for(auto& t : i.transitions_) {
        const INDEX timestep = std::get<0>(t);
        const INDEX prev_cell = std::get<1>(t);
        const INDEX next_cell = std::get<2>(t);
        const REAL cost = std::get<3>(t);
        cell_tracking_constructor.add_cell_transition( lp, timestep, prev_cell, timestep + 1, next_cell, cost);
      }
      for(auto& t : i.divisions_) {
        const INDEX timestep = std::get<0>(t);
        const INDEX prev_cell = std::get<1>(t);
        const INDEX next_cell_1 = std::get<2>(t);
        const INDEX next_cell_2 = std::get<3>(t);
        const REAL cost = std::get<4>(t);
        cell_tracking_constructor.add_cell_division( lp, timestep, prev_cell, timestep+1, next_cell_1, timestep+1, next_cell_2, cost);
      }
      //for(INDEX t=0; t<i.cell_detection_stat_.size(); ++t) {
      //  for(INDEX j=0; j<i.cell_detection_stat_[t].size(); ++j) {
      //    assert( std::get<2>(i.cell_detection_stat_[t][j]) == tc[t][j][0] );
      //    assert( std::get<3>(i.cell_detection_stat_[t][j]) == tc[t][j][1] );
      //  }
      //}
   }

   // we assume problem with single constructor
   template<typename SOLVER>
   bool ParseProblem(const std::string& filename, SOLVER& s)
   {
     input i;
     const bool read_suc = read_input(filename, i);
     construct_tracking_problem(i, s);
     return read_suc;
   }

   template<typename SOLVER>
   bool ParseProblemMotherMachine(const std::string& filename, SOLVER& s)
   {
     input i;
     const bool read_suc = read_input(filename, i);
     construct_tracking_problem(i, s);

      auto& cell_tracking_constructor = s.template GetProblemConstructor<0>();
      // add exit constraints based on positions of cell detection hypotheses
      // coordinates are in format (upper,lower) with upper < lower (hence coordinates are inverted!)
      //std::cout << "Exit constraints disabled\n";
      for(INDEX t=0; t<i.cell_detection_stat_.size(); ++t) {
        // possibly better: if there exists cell strictly between i and j, then no exit constraint needs to be added
        for(INDEX d1=0; d1<i.cell_detection_stat_[t].size(); ++d1) {
          for(INDEX d2=d1+1; d2<i.cell_detection_stat_[t].size(); ++d2) {
            auto upper_bound_d1 = std::get<4>(i.cell_detection_stat_[t][d1]);
            auto lower_bound_d1 = std::get<5>(i.cell_detection_stat_[t][d1]);
            auto upper_bound_d2 = std::get<4>(i.cell_detection_stat_[t][d2]);
            auto lower_bound_d2 = std::get<5>(i.cell_detection_stat_[t][d2]);
            assert(upper_bound_d1 < lower_bound_d1);
            assert(upper_bound_d2 < lower_bound_d2);
            //std::cout << "t=" << t << ", H1=" << d1 << "(" << lower_bound_d1 << "," << upper_bound_d1 << "), H2=" << d2 << "(" << lower_bound_d2 << "," << upper_bound_d2 << ")\n";
            if(upper_bound_d1 > lower_bound_d2) {
              cell_tracking_constructor.add_exit_constraint(s.GetLP(), t,d1,d2);
            }
            if(upper_bound_d2 > lower_bound_d1) {
              cell_tracking_constructor.add_exit_constraint(s.GetLP(), t,d2,d1);
            }
          }
        } 
      }

      cell_tracking_constructor.order_factors(s.GetLP());

      INDEX max_out_degree = 0;
      INDEX max_in_degree = 0;
      INDEX sum_out_degree = 0;
      INDEX sum_in_degree = 0;
      INDEX no_factors = 0;
      for(auto& d : i.cell_detection_stat_) {
        for(auto& f : d) {
          max_out_degree = std::max(max_out_degree, std::get<2>(f));
          max_in_degree = std::max(max_in_degree, std::get<3>(f));

          sum_out_degree += std::get<2>(f);
          sum_in_degree += std::get<3>(f);
          if(std::get<2>(f) > 0 || std::get<3>(f) > 0) no_factors++;
        }
      }
      if(verbosity >= 2) {
        std::cout << "maximum out-degree of detection factors = " << max_out_degree << "\n";
        std::cout << "maximum in-degree of detection factors = " << max_in_degree << "\n";

        std::cout << "average out-degree of detection factors = " << REAL(sum_out_degree)/REAL(no_factors) << "\n";
        std::cout << "average in-degree of detection factors = " << REAL(sum_in_degree)/REAL(no_factors) << "\n";
      }

      return read_suc;
   }
} // end namespace cell_tracking_parser_mother_machine


namespace cell_tracking_parser_2d {

   // import basic parsers
   using Parsing::opt_whitespace;
   using Parsing::mand_whitespace;
   using Parsing::opt_invisible;
   using Parsing::mand_invisible;
   using Parsing::positive_integer;
   using Parsing::real_number;

   struct my_eolf : pegtl::sor<pegtl::eol, pegtl::eof> {}; // do zrobienia: revert to usual pegtl::eolf
   struct comment_line : pegtl::seq< opt_whitespace, pegtl::sor< pegtl::seq< pegtl::string<'#'>, pegtl::until<my_eolf> >, my_eolf> > {};

   struct division_distance : pegtl::seq< positive_integer> {};
   struct division_distance_line : pegtl::seq< opt_whitespace, pegtl::string<'D','I','V','I','S','I','O','N',' ','D','I','S','T','A','N','C','E'>, opt_whitespace, pegtl::string<'='>, opt_whitespace, division_distance, opt_whitespace, my_eolf > {};
   
   // H timestep hypothesis_id {detection cost} (x_coord, lower coord) 
   struct cell_detection_hypothesis : pegtl::seq< positive_integer, mand_whitespace, positive_integer, mand_whitespace, real_number, opt_whitespace, pegtl::string<'('>, opt_whitespace, real_number, opt_whitespace, pegtl::string<','>, opt_whitespace, real_number, opt_whitespace, pegtl::string<')'> > {};
   struct cell_detection_hypothesis_line : pegtl::seq< pegtl::string<'H'>, opt_whitespace, cell_detection_hypothesis, opt_whitespace, my_eolf > {};
   
   // appearance cost: timestep hypothesis_id cost
   struct appearance : pegtl::seq<positive_integer, mand_whitespace, positive_integer, mand_whitespace, real_number > {};
   struct appearance_line : pegtl::seq< pegtl::string<'A','P','P'>, mand_whitespace, appearance, opt_whitespace, my_eolf > {};
   // disappearance cost: timestep hypothesis_id cost
   struct disappearance : pegtl::seq<positive_integer, mand_whitespace, positive_integer, mand_whitespace, real_number > {};
   struct disappearance_line : pegtl::seq< pegtl::string<'D','I','S','A','P','P'>, mand_whitespace, disappearance, opt_whitespace, my_eolf > {};
   // CONFSET timestep_1 hyp_no1 + timestep hyp_no2 ... <= 1
   struct conflict : pegtl::seq< positive_integer, mand_whitespace, positive_integer, pegtl::star< pegtl::seq< opt_whitespace, pegtl::string<'+'>, opt_whitespace, positive_integer, mand_whitespace, positive_integer> > > {};
   struct conflict_line : pegtl::seq< pegtl::string<'C','O','N','F','S','E','T'>, opt_whitespace, conflict, opt_whitespace, pegtl::string<'<','='>, opt_whitespace, pegtl::string<'1'>, opt_whitespace, my_eolf> {};
   // timestep 1, first hypothesis id, timestep 1, second hypothesis id, cost
   struct mapping : pegtl::seq< positive_integer, mand_whitespace, positive_integer, mand_whitespace, positive_integer, mand_whitespace, positive_integer, mand_whitespace, real_number > {}; 
   struct mapping_line : pegtl::seq< pegtl::string<'M','O','V','E'>, opt_whitespace, mapping, opt_whitespace, my_eolf> {};
   // timestep 1, left hypothesis id, right timestep 1, right hypothesid 1, right timestep 2, right hypothesis id 2, cost
   struct division : pegtl::seq< positive_integer, mand_whitespace, positive_integer, mand_whitespace, positive_integer, mand_whitespace, positive_integer, mand_whitespace, positive_integer, mand_whitespace, positive_integer, mand_whitespace, real_number > {}; 
   struct division_line : pegtl::seq< pegtl::string<'D','I','V'>, opt_whitespace, division, opt_whitespace, my_eolf> {};

   struct grammar : pegtl::seq< pegtl::star<pegtl::sor<
                    comment_line,
                    cell_detection_hypothesis_line,
                    appearance_line,
                    disappearance_line,
                    conflict_line,
                    mapping_line,
                    division_line
                    >>> {};

   struct grammar_division_distance : 
                    pegtl::seq< 
                    pegtl::star< comment_line >,
                    division_distance_line,
                    
                    pegtl::star<pegtl::sor<
                    comment_line,
                    cell_detection_hypothesis_line,
                    appearance_line,
                    disappearance_line,
                    conflict_line,
                    mapping_line,
                    division_line
                    >>> {};


   struct input {
     struct detection_factor_stat {
       REAL detection_cost = 0.0, appearance_cost = 0.0, disappearance_cost = 0.0;
       INDEX no_incoming_transition_edges = 0, no_incoming_division_edges = 0, no_outgoing_transition_edges = 0, no_outgoing_division_edges = 0;
     };
     // we must go through all cell mappings (cell transitions and divisions) to count the number of outgoing and incoming messages. Only then can we allocate cell detection hypothesis, as they need to know these numbers.
     std::vector<std::vector<detection_factor_stat>> cell_detection_stat; // detection cost, appearance cost, disappearance cost,  # incoming edges, # outgoing edges
     std::vector<std::tuple<INDEX,INDEX,INDEX,INDEX,REAL>> mappings; // timestep outgoing, outgoing cell, timestep incoming, incoming cell, cost
     std::vector<std::tuple<INDEX,INDEX,INDEX,INDEX,INDEX,INDEX,REAL>> divisions; // timestep outgoing, outgoing cell, incoming timestep 1, incoming cell 1, incoming timestep 2, incoming cell 2, cost
     std::vector<std::vector<std::array<INDEX,2>>> conflicts; // timestep, {hyp_1, ..., hyp_n}
     INDEX division_distance = 1;
     bool eof = false;
   };

   template< typename Rule >
     struct action
     : pegtl::nothing< Rule > {};

   template<> struct action< pegtl::eof > {
     template<typename INPUT>
     static void apply(const INPUT & in, input& i)
     {
       //std::cout << "kwaskwas\n";
       i.eof = true;
     }
   };

   template<> struct action< division_distance > {
     template<typename INPUT>
     static void apply(const INPUT & in, input& i)
     {
       assert(i.division_distance == 1);
       std::stringstream s(in.string()); 
       s >> i.division_distance; 
       assert(i.division_distance > 0);
     }
   };

   template<> struct action< cell_detection_hypothesis > {
     template<typename INPUT>
     static void apply(const INPUT & in, input& i)
     {
       std::stringstream s(in.string());
       INDEX timestep; s >> timestep;
       INDEX hypothesis_id; s >> hypothesis_id; 
       REAL detection_cost; s >> detection_cost;
       char opening_bracket; s >> opening_bracket; assert(opening_bracket == '(');
       REAL x; s >> x;
       char comma; s >> comma; assert(comma == ',');
       REAL y; s >> y;
       char closing_bracket; s >> closing_bracket; assert(closing_bracket == ')');

       if(i.cell_detection_stat.size() <= timestep) {
         i.cell_detection_stat.resize(timestep+1);
       }
       if(i.cell_detection_stat[timestep].size() <= hypothesis_id) {
         i.cell_detection_stat[timestep].resize(hypothesis_id+1);
       }
       i.cell_detection_stat[timestep][hypothesis_id].detection_cost = detection_cost;
       //std::cout << "cell deteciton hypothesis " << timestep << "," << hypothesis_id << "\n";
     }
   };

   template<> struct action< appearance > {
     template<typename INPUT>
       static void apply(const INPUT & in, input& i)
       {
         std::stringstream s(in.string());
         INDEX timestep; s >> timestep;
         INDEX hypothesis_id; s >> hypothesis_id;
         REAL cost; s >> cost;
         i.cell_detection_stat[timestep][hypothesis_id].appearance_cost = cost;
         //std::cout << "APP " << timestep << " " << hypothesis_id << " " << cost << "\n";
       }
   };

   template<> struct action< disappearance > {
     template<typename INPUT>
       static void apply(const INPUT & in, input& i)
       {
         std::stringstream s(in.string());
         INDEX timestep; s >> timestep;
         INDEX hypothesis_id; s >> hypothesis_id;
         REAL cost; s >> cost;
         i.cell_detection_stat[timestep][hypothesis_id].disappearance_cost = cost;
         //std::cout << "DISAPP " << timestep << " " << hypothesis_id << " " << cost << "\n";
       }
   };

   template<> struct action< conflict > {
     template<typename INPUT>
     static void apply(const INPUT & in, input& i)
     {
       std::stringstream s(in.string());
       std::vector<std::array<INDEX,2>> idx; // timestep, hypothesis id
       INDEX timestep; s >> timestep; 
       INDEX hypothesis_id; s >> hypothesis_id; 
       idx.push_back({timestep, hypothesis_id});
       char plus;
       while(s >> plus) {
         assert(plus == '+');
         INDEX timestep; s >> timestep; 
         INDEX hypothesis_id; s >> hypothesis_id; 
         idx.push_back({timestep, hypothesis_id});
       }
       i.conflicts.push_back(idx);
       //std::cout << "CONFSET ";
       //for(INDEX i=0; i<idx.size(); ++i) {
       //  std::cout << idx[i][0] << " " << idx[i][1] << " + ";
       //}
       //std::cout << " <= 1\n";
     }
   };

   template<> struct action< mapping > {
     template<typename INPUT>
     static void apply(const INPUT & in, input& i)
     {
       std::stringstream s(in.string());
       INDEX timestep_prev; s >> timestep_prev;
       INDEX cell_prev; s >> cell_prev;
       INDEX timestep_next; s >> timestep_next;
       INDEX cell_next; s >> cell_next;
       REAL cost; s >> cost;

       assert(timestep_prev+1 == timestep_next);
       assert(timestep_next < i.cell_detection_stat.size());
       assert(cell_prev < i.cell_detection_stat[timestep_prev].size());
       assert(cell_next < i.cell_detection_stat[timestep_next].size());

       i.cell_detection_stat[timestep_prev][cell_prev].no_outgoing_transition_edges++;
       i.cell_detection_stat[timestep_next][cell_next].no_incoming_transition_edges++;

       i.mappings.push_back( std::make_tuple(timestep_prev, cell_prev, timestep_next, cell_next, cost));

       //std::cout << "mapping: t = " << timestep_prev << ", h1 = " << cell_prev << ", h2 = " << cell_next << ", cost = " << cost << "\n";
     }
   };

   template<> struct action< division > {
     template<typename INPUT>
     static void apply(const INPUT & in, input& i)
     {
       std::stringstream s(in.string());
       INDEX timestep_prev; s >> timestep_prev;
       INDEX cell_prev; s >> cell_prev;
       INDEX timestep_next_1; s >> timestep_next_1;
       INDEX cell_next_1; s >> cell_next_1;
       INDEX timestep_next_2; s >> timestep_next_2;
       INDEX cell_next_2; s >> cell_next_2;
       REAL cost; s >> cost;

       assert(timestep_prev+1 == timestep_next_1);
       assert(timestep_prev+1 == timestep_next_2);
       assert(cell_next_1 != cell_next_2);
       assert(timestep_next_1 < i.cell_detection_stat.size());
       assert(timestep_next_2 < i.cell_detection_stat.size());
       assert(cell_prev < i.cell_detection_stat[timestep_prev].size());
       assert(cell_next_1 < i.cell_detection_stat[timestep_next_1].size());
       assert(cell_next_2 < i.cell_detection_stat[timestep_next_2].size());

       i.cell_detection_stat[timestep_prev][cell_prev].no_outgoing_division_edges++;
       i.cell_detection_stat[timestep_next_1][cell_next_1].no_incoming_division_edges++;
       i.cell_detection_stat[timestep_next_2][cell_next_2].no_incoming_division_edges++;

       i.divisions.push_back( std::make_tuple(timestep_prev, cell_prev, timestep_next_1, cell_next_1, timestep_next_2, cell_next_2, cost));
       //std::cout << "division: t = " << timestep_prev << ", h1 = " << cell_prev << ", h1,1 = " << cell_next_1 << ", h2,2 = " << cell_next_2 << ", cost = " << cost << "\n";
     }
   };


   template<typename GRAMMAR>
   bool read_input(const std::string& filename, input& i)
   {
      if(verbosity >= 1) {
        std::cout << "parsing " << filename << "\n";
      }
      pegtl::file_parser problem(filename);
      const bool success = problem.parse< GRAMMAR, action >(i); 
      //assert(i.eof);
      return success;
   }

   template<typename SOLVER>
   void construct_tracking_problem(input& i, SOLVER& s)
   {
      auto& cell_tracking_constructor = s.template GetProblemConstructor<0>();
      auto& lp = s.GetLP();

      //std::cout << "exclusion constraints disabled\n";
      cell_tracking_constructor.set_number_of_timesteps( i.cell_detection_stat.size() );
      for(INDEX t=0; t<i.cell_detection_stat.size(); ++t) {
        for(INDEX n=0; n<i.cell_detection_stat[t].size(); ++n) {
          const REAL detection_cost = i.cell_detection_stat[t][n].detection_cost;
          const REAL appearance_cost = i.cell_detection_stat[t][n].appearance_cost;
          const REAL disappearance_cost = i.cell_detection_stat[t][n].disappearance_cost;
          const INDEX no_incoming_transition_edges = i.cell_detection_stat[t][n].no_incoming_transition_edges;
          const INDEX no_outgoing_transition_edges = i.cell_detection_stat[t][n].no_outgoing_transition_edges;
          const INDEX no_incoming_division_edges = i.cell_detection_stat[t][n].no_incoming_division_edges;
          const INDEX no_outgoing_division_edges = i.cell_detection_stat[t][n].no_outgoing_division_edges;

          cell_tracking_constructor.add_detection_hypothesis( lp, t, n, detection_cost, appearance_cost, disappearance_cost, no_incoming_transition_edges, no_incoming_division_edges, no_outgoing_transition_edges, no_outgoing_division_edges);
        }
      }
      for(const auto& conflict_set : i.conflicts) {
        cell_tracking_constructor.add_exclusion_constraint(lp, conflict_set.begin(), conflict_set.end()); 
        //cell_tracking_constructor.register_exclusion_constraint(conflict_set.begin(), conflict_set.end() ); 
      }

      //auto tc = cell_tracking_constructor.init_transition_counter();
      // the order is important! Possibly change and treat transition and division edges separately
      for(auto& t : i.mappings) {
        const INDEX timestep_prev = std::get<0>(t);
        const INDEX prev_cell = std::get<1>(t);
        const INDEX timestep_next = std::get<2>(t);
        const INDEX next_cell = std::get<3>(t);
        const REAL cost = std::get<4>(t);
        cell_tracking_constructor.add_cell_transition( lp, timestep_prev, prev_cell, timestep_next, next_cell, cost);
      }
      for(auto& t : i.divisions) {
        const INDEX timestep_prev = std::get<0>(t);
        const INDEX prev_cell = std::get<1>(t);
        const INDEX timestep_next_1 = std::get<2>(t);
        const INDEX next_cell_1 = std::get<3>(t);
        const INDEX timestep_next_2 = std::get<4>(t);
        const INDEX next_cell_2 = std::get<5>(t);
        const REAL cost = std::get<6>(t);
        cell_tracking_constructor.add_cell_division( lp, timestep_prev, prev_cell, timestep_next_1, next_cell_1, timestep_next_2, next_cell_2, cost);
      }
      //for(INDEX t=0; t<i.cell_detection_stat.size(); ++t) {
      //  for(INDEX j=0; j<i.cell_detection_stat[t].size(); ++j) {
      //    assert( i.cell_detection_stat[t][j].no_outgoing_transition_edges + i.cell_detection_stat[t][j].no_outgoing_division_edges == tc.current_transition_no[t][j][1] + tc.current_division_no[t][j][1] );
      //    assert( i.cell_detection_stat[t][j].no_incoming_transition_edges + i.cell_detection_stat[t][j].no_incoming_division_edges == tc.current_transition_no[t][j][0] + tc.current_division_no[t][j][0] );
      //  }
      //}

      cell_tracking_constructor.order_factors(lp);
   }

   // we assume problem with single constructor
   template<typename SOLVER>
   bool ParseProblem(const std::string& filename, SOLVER& s)
   {
     input i;
     const bool read_suc = read_input<grammar>(filename, i);
     construct_tracking_problem(i, s);
     return read_suc;
   }

   // we assume problem with single constructor
   template<typename SOLVER>
   bool parse_problem_with_division_distance(const std::string& filename, SOLVER& s)
   {
     input i;
     const bool read_suc = read_input<grammar_division_distance>(filename, i);
     auto& cell_tracking_constructor = s.template GetProblemConstructor<0>();
     cell_tracking_constructor.set_division_distance(i.division_distance);
     construct_tracking_problem(i, s);
     return read_suc;
   }

} // end namespace cell_tracking_parser_2d




} // end namespace LP_MP

#endif // LP_MP_CELL_TRACKING_CONSTRUCTOR_HXX

