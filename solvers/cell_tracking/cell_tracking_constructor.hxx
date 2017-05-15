#ifndef LP_MP_CELL_TRACKING_CONSTRUCTOR_HXX
#define LP_MP_CELL_TRACKING_CONSTRUCTOR_HXX

#include <vector>
#include "solver.hxx"
#include "pegtl/parse.hh"
#include "parse_rules.h"

// do zrobienia: do not include detection hypotheses that have no incoming and no outgoing edges into the LP

namespace LP_MP {

template<typename DETECTION_FACTOR_CONTAINER, typename AT_MOST_ONE_CELL_FACTOR_CONTAINER,
         typename TRANSITION_MESSAGE_CONTAINER, typename AT_MOST_ONE_CELL_MESSAGE_CONTAINER>
class cell_tracking_constructor {
public:
  using CONSTRUCTOR = cell_tracking_constructor<DETECTION_FACTOR_CONTAINER, AT_MOST_ONE_CELL_FACTOR_CONTAINER, TRANSITION_MESSAGE_CONTAINER, AT_MOST_ONE_CELL_MESSAGE_CONTAINER>;
  template<typename SOLVER>
  cell_tracking_constructor(SOLVER& solver) {}

  // temporary structure which counts how many incoming and outgoing edges are already used by messages for building the model
  struct transition_count {
    transition_count(CONSTRUCTOR& c) 
    {
      current_transition_no.resize( c.detection_factors_.size() );
      for(INDEX i=0; i<c.detection_factors_.size(); ++i) {
        current_transition_no[i].resize( c.detection_factors_[i].size(), {0,0});
      }
    }
    std::vector<std::array<INDEX,2>>& operator[](const INDEX i) { return current_transition_no[i]; }
    std::vector<std::vector<std::array<INDEX,2>>> current_transition_no; 
  };
  transition_count init_transition_counter() { return transition_count(*this); }

  void set_number_of_timesteps(const INDEX t)
  {
    assert(detection_factors_.size() == 0);
    detection_factors_.resize(t);
  }

  template<typename LP_TYPE>
  DETECTION_FACTOR_CONTAINER* add_detection_hypothesis(LP_TYPE& lp, const INDEX timestep, const INDEX number, const REAL detection_cost, const INDEX no_incoming_edges, const INDEX no_outgoing_edges, const REAL exit_cost)
  {

    assert(timestep < detection_factors_.size());
    if(number >= detection_factors_[timestep].size()) {
      detection_factors_[timestep].resize(number+1, nullptr);
    }
    assert(detection_factors_[timestep][number] == nullptr);

    if(detection_cost == 0.0 && no_incoming_edges == 0 && no_outgoing_edges == 0 && exit_cost == 0.0) { return nullptr; }

    auto* f = new DETECTION_FACTOR_CONTAINER(std::max(no_incoming_edges,INDEX(1)), no_outgoing_edges, detection_cost, true);
    f->GetFactor()->outgoing( f->GetFactor()->no_outgoing_edges()-1 ) = exit_cost;
    lp.AddFactor(f);
    detection_factors_[timestep][number] = f;
    //std::cout << "Added ";
    //std::cout << "H: " << timestep << ", " << number << ", " << detection_cost << std::endl;
    return f; 
  }

  template<typename LP_TYPE>
  AT_MOST_ONE_CELL_FACTOR_CONTAINER* add_exclusion_constraint(LP_TYPE& lp, const INDEX timestep, const std::vector<INDEX>& cell_detections)
  {
    INDEX size = 0;
    for(INDEX i=0; i<cell_detections.size(); ++i) {
      if(detection_factors_[timestep][ cell_detections[i] ] != nullptr) {
        ++size; 
      }
    }
    if(size == 0) { return nullptr; }

    auto* e = new AT_MOST_ONE_CELL_FACTOR_CONTAINER(size);
    lp.AddFactor(e);
    //std::cout << "Added Excusion for time " << timestep << ": ";
    INDEX msg_idx = 0;
    for(INDEX i=0; i<cell_detections.size(); ++i) {
      //std::cout << cell_detections[i] << ", ";
      if(detection_factors_[timestep][ cell_detections[i] ] != nullptr) {
        auto* m = new AT_MOST_ONE_CELL_MESSAGE_CONTAINER(msg_idx, detection_factors_[ timestep ][ cell_detections[i] ], e);
        ++msg_idx;
        lp.AddMessage(m); 
      }
    }
    //std::cout << std::endl;
    return e; 
  }

  template<typename LP_TYPE>
  void add_cell_transition(LP_TYPE& lp, const INDEX timestep, const INDEX prev_cell, const INDEX next_cell, const REAL cost, transition_count& tc) 
  {
    auto* out_cell_factor = detection_factors_[timestep][prev_cell];
    const INDEX outgoing_edge_index  = tc[timestep][prev_cell][1];
    assert( out_cell_factor->GetFactor()->outgoing(outgoing_edge_index) == 0.0 );
    out_cell_factor->GetFactor()->outgoing(outgoing_edge_index) = 0.5*cost;
    //out_cell_factor->GetFactor()->outgoing(outgoing_edge_index) = cost;
    tc[timestep][prev_cell][1]++;

    auto* in_cell_factor = detection_factors_[timestep+1][next_cell];
    const INDEX incoming_edge_index = tc[timestep+1][next_cell][0];
    assert( in_cell_factor->GetFactor()->incoming(incoming_edge_index) == 0.0 );
    tc[timestep+1][next_cell][0]++;
    in_cell_factor->GetFactor()->incoming(incoming_edge_index) = 0.5*cost;
    auto* m = new TRANSITION_MESSAGE_CONTAINER(out_cell_factor, in_cell_factor, false, outgoing_edge_index, incoming_edge_index);
    lp.AddMessage(m);

    //std::cout << "MA: " << timestep << " " << prev_cell << ", " << next_cell << " " << cost << std::endl;
  }

  template<typename LP_TYPE>
  void add_cell_division(LP_TYPE& lp, const INDEX timestep, const INDEX prev_cell, const INDEX next_cell_1, const INDEX next_cell_2, const REAL cost, transition_count& tc) 
  {
    auto* out_cell_factor = detection_factors_[timestep][prev_cell];
    const INDEX outgoing_edge_index  = tc[timestep][prev_cell][1];
    assert( out_cell_factor->GetFactor()->outgoing(outgoing_edge_index) == 0.0 );
    out_cell_factor->GetFactor()->outgoing(outgoing_edge_index) = 1.0/3.0*cost;
    //out_cell_factor->GetFactor()->outgoing(outgoing_edge_index) = cost;
    tc[timestep][prev_cell][1]++;

    auto* in_cell_factor_1 = detection_factors_[timestep+1][next_cell_1];
    const INDEX incoming_edge_index_1 = tc[timestep+1][next_cell_1][0];
    assert( in_cell_factor_1->GetFactor()->incoming(incoming_edge_index_1) == 0.0 );
    in_cell_factor_1->GetFactor()->incoming(incoming_edge_index_1) = 1.0/3.0*cost;
    tc[timestep+1][next_cell_1][0]++;
    
    auto* in_cell_factor_2 = detection_factors_[timestep+1][next_cell_2];
    const INDEX incoming_edge_index_2 = tc[timestep+1][next_cell_2][0];
    assert( in_cell_factor_2->GetFactor()->incoming(incoming_edge_index_2) == 0.0 );
    in_cell_factor_2->GetFactor()->incoming(incoming_edge_index_2) = 1.0/3.0*cost;
    tc[timestep+1][next_cell_2][0]++;
    
    auto* m1 = new TRANSITION_MESSAGE_CONTAINER(out_cell_factor, in_cell_factor_1, true, outgoing_edge_index, incoming_edge_index_1);
    lp.AddMessage(m1);
    
    auto* m2 = new TRANSITION_MESSAGE_CONTAINER(out_cell_factor, in_cell_factor_2, true, outgoing_edge_index, incoming_edge_index_2);
    lp.AddMessage(m2);

    //std::cout << "DA: " << timestep << " " << prev_cell << ", " << next_cell_1 << " " << next_cell_2 << " " << cost << std::endl;
  }

  // output lineage tree
  template<typename STREAM>
  void WritePrimal(STREAM& s, PrimalSolutionStorage& primal) const 
  {
    return;
    for(INDEX t=0; t<detection_factors_.size(); ++t) {
      std::cout << "t=" << t << "\n";
      for(INDEX i=0; i<detection_factors_[t].size(); ++i) {
        auto* f = detection_factors_[t][i]->GetFactor();
        if(f->outgoing_edge() < f->no_outgoing_edges()) {
          std::cout << i << " ";
        }
      }
      std::cout << "\n";
    }
  }

protected:
  using detection_factors_storage = std::vector<std::vector<DETECTION_FACTOR_CONTAINER*>>;
  detection_factors_storage detection_factors_;
  //std::map<std::tuple<INDEX,INDEX,INDEX>, TRANSITION_FACTOR_CONTAINER*> transition_factors_;
  //std::map<std::tuple<INDEX,INDEX,INDEX,INDEX>, TRANSITION_FACTOR_CONTAINER*> division_factors_;

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
      std::cout << "parsing " << filename << "\n";
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
          cell_tracking_constructor.add_detection_hypothesis( lp, t, n, detection_cost, no_incoming_edges, no_outgoing_edges, exit_cost );
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

      auto tc = cell_tracking_constructor.init_transition_counter();
      for(auto& t : i.transitions_) {
        const INDEX timestep = std::get<0>(t);
        const INDEX prev_cell = std::get<1>(t);
        const INDEX next_cell = std::get<2>(t);
        const REAL cost = std::get<3>(t);
        cell_tracking_constructor.add_cell_transition( lp, timestep, prev_cell, next_cell, cost, tc );
      }
      for(auto& t : i.divisions_) {
        const INDEX timestep = std::get<0>(t);
        const INDEX prev_cell = std::get<1>(t);
        const INDEX next_cell_1 = std::get<2>(t);
        const INDEX next_cell_2 = std::get<3>(t);
        const REAL cost = std::get<4>(t);
        cell_tracking_constructor.add_cell_division( lp, timestep, prev_cell, next_cell_1, next_cell_2, cost, tc );
      }
      for(INDEX t=0; t<i.cell_detection_stat_.size(); ++t) {
        for(INDEX j=0; j<i.cell_detection_stat_[t].size(); ++j) {
          assert( std::get<2>(i.cell_detection_stat_[t][j]) == tc[t][j][0] );
          assert( std::get<3>(i.cell_detection_stat_[t][j]) == tc[t][j][1] );
        }
      }
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
    std::cout << "maximum out-degree of detection factors = " << max_out_degree << "\n";
    std::cout << "maximum in-degree of detection factors = " << max_in_degree << "\n";

    std::cout << "average out-degree of detection factors = " << REAL(sum_out_degree)/REAL(no_factors) << "\n";
    std::cout << "average in-degree of detection factors = " << REAL(sum_in_degree)/REAL(no_factors) << "\n";

      return read_suc;
   }
} // end namespace cell_tracking_parser_mother_machine


/*
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
   
   // H timestep hypothesis_id {detection cost} (x_coord, lower coord) 
   struct cell_detection_hypothesis : pegtl::seq< positive_integer, mand_whitespace, positive_integer, mand_whitespace, real_number, opt_whitespace, pegtl::string<'('>, opt_whitespace, real_number, opt_whitespace, pegtl::string<','>, opt_whitespace, real_number, opt_whitespace, pegtl::string<')'> > {};
   struct cell_detection_hypothesis_line : pegtl::seq< pegtl::string<'H'>, opt_whitespace, cell_detection_hypothesis, opt_whitespace, my_eolf > {};
   
   // appearance cost: timestep hypothesis_id cost
   struct appearance_cost : pegtl::seq< pegtl::string<'A','P','P'>, mand_whitespace, positive_integer
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
      std::cout << "parsing " << filename << "\n";
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
          cell_tracking_constructor.add_detection_hypothesis( lp, t, n, detection_cost, no_incoming_edges, no_outgoing_edges, exit_cost );
        }
        for(const auto& detections : i.exclusions_[t]) {
          cell_tracking_constructor.add_exclusion_constraint(lp, t, detections ); 
        }
        // check whether all cell hypotheses are covered by exclusion constraints (must be valid for mother machine)
        //std::vector<INDEX> all_indices;
        //for(const auto& detections : i.exclusions_[t]) {
        //  all_indices.insert(all_indices.end(), detections.begin(), detections.end());
        //}
        //std::sort(all_indices.begin(), all_indices.end());
        //all_indices.erase( std::unique( all_indices.begin(), all_indices.end()), all_indices.end() );
        //assert(all_indices.size() == i.cell_detection_stat_[t].size());
        
      }

      auto tc = cell_tracking_constructor.init_transition_counter();
      for(auto& t : i.transitions_) {
        const INDEX timestep = std::get<0>(t);
        const INDEX prev_cell = std::get<1>(t);
        const INDEX next_cell = std::get<2>(t);
        const REAL cost = std::get<3>(t);
        cell_tracking_constructor.add_cell_transition( lp, timestep, prev_cell, next_cell, cost, tc );
      }
      for(auto& t : i.divisions_) {
        const INDEX timestep = std::get<0>(t);
        const INDEX prev_cell = std::get<1>(t);
        const INDEX next_cell_1 = std::get<2>(t);
        const INDEX next_cell_2 = std::get<3>(t);
        const REAL cost = std::get<4>(t);
        cell_tracking_constructor.add_cell_division( lp, timestep, prev_cell, next_cell_1, next_cell_2, cost, tc );
      }
      for(INDEX t=0; t<i.cell_detection_stat_.size(); ++t) {
        for(INDEX j=0; j<i.cell_detection_stat_[t].size(); ++j) {
          assert( std::get<2>(i.cell_detection_stat_[t][j]) == tc[t][j][0] );
          assert( std::get<3>(i.cell_detection_stat_[t][j]) == tc[t][j][1] );
        }
      }
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
    std::cout << "maximum out-degree of detection factors = " << max_out_degree << "\n";
    std::cout << "maximum in-degree of detection factors = " << max_in_degree << "\n";

    std::cout << "average out-degree of detection factors = " << REAL(sum_out_degree)/REAL(no_factors) << "\n";
    std::cout << "average in-degree of detection factors = " << REAL(sum_in_degree)/REAL(no_factors) << "\n";

      return read_suc;
   }
} // end namespace cell_tracking_parser_2d
*/




} // end namespace LP_MP

#endif // LP_MP_CELL_TRACKING_CONSTRUCTOR_HXX

