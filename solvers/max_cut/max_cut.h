#ifndef LP_MP_MAX_CUT_H
#define LP_MP_MAX_CUT_H

#include "LP_MP.h"
#include "solver.hxx"
#include "factors_messages.hxx"
#include "max_cut_factors.hxx"
#include "factors/constant_factor.hxx"
#include "../multicut/multicut_constructor.hxx"

#include "parse_rules.h"

#include <iostream>
#include <vector>

namespace LP_MP {

struct FMC_CYCLE_MAX_CUT {
   constexpr static const char* name = "Max cut with cycle constraints";

   using max_cut_edge_factor_container = FactorContainer<max_cut_edge_factor, FMC_CYCLE_MAX_CUT, 0, true>;
   using max_cut_triplet_factor_container = FactorContainer<max_cut_triplet_factor, FMC_CYCLE_MAX_CUT, 1>;
   using ConstantFactorContainer = FactorContainer<ConstantFactor, FMC_CYCLE_MAX_CUT, 2>;
      
   using max_cut_edge_triplet_message_0_container = MessageContainer<max_cut_edge_triplet_message_0, 0, 1, variableMessageNumber, 1, FMC_CYCLE_MAX_CUT, 0>;
   using max_cut_edge_triplet_message_1_container = MessageContainer<max_cut_edge_triplet_message_1, 0, 1, variableMessageNumber, 1, FMC_CYCLE_MAX_CUT, 1>;
   using max_cut_edge_triplet_message_2_container = MessageContainer<max_cut_edge_triplet_message_2, 0, 1, variableMessageNumber, 1, FMC_CYCLE_MAX_CUT, 2>;

   using FactorList = meta::list< 
      max_cut_edge_factor_container,
      max_cut_triplet_factor_container,
      ConstantFactorContainer 
         >;

   using MessageList = meta::list<
      max_cut_edge_triplet_message_0_container,
      max_cut_edge_triplet_message_1_container,
      max_cut_edge_triplet_message_2_container
      >;

   using max_cut_c = MulticutConstructor<FMC_CYCLE_MAX_CUT,0,1,0,1,2,2,cut_type::maxcut>;
   using ProblemDecompositionList = meta::list<max_cut_c>;
};

// It would be nice to be able to derive from FMC_MULTICUT. This is not possible due to deviating FMCs. Possibly parametrize above FMC with template
struct FMC_ODD_BICYCLE_WHEEL_MAX_CUT {
   constexpr static const char* name = "Max cut with cycle and odd bicycle wheel constraints";

   using edge_factor_container = FactorContainer<max_cut_edge_factor, FMC_ODD_BICYCLE_WHEEL_MAX_CUT, 0>;
   using triplet_factor_container = FactorContainer<max_cut_triplet_factor, FMC_ODD_BICYCLE_WHEEL_MAX_CUT, 1>;
   using four_clique_factor_container = FactorContainer<max_cut_4_clique_factor, FMC_ODD_BICYCLE_WHEEL_MAX_CUT, 2>; 
   using five_clique_factor_container = FactorContainer<max_cut_5_clique_factor, FMC_ODD_BICYCLE_WHEEL_MAX_CUT, 3>;
   using constant_factor_container = FactorContainer<ConstantFactor, FMC_ODD_BICYCLE_WHEEL_MAX_CUT, 4>;
      
   using edge_triplet_message_0_container = MessageContainer<max_cut_edge_triplet_message_0, 0, 1, variableMessageNumber, 1, FMC_ODD_BICYCLE_WHEEL_MAX_CUT, 0>;
   using edge_triplet_message_1_container = MessageContainer<max_cut_edge_triplet_message_1, 0, 1, variableMessageNumber, 1, FMC_ODD_BICYCLE_WHEEL_MAX_CUT, 1>;
   using edge_triplet_message_2_container = MessageContainer<max_cut_edge_triplet_message_2, 0, 1, variableMessageNumber, 1, FMC_ODD_BICYCLE_WHEEL_MAX_CUT, 2>;

   using triplet_4_clique_message_012_container = MessageContainer<max_cut_triplet_4_clique_message_012, 1, 2, variableMessageNumber, 1, FMC_ODD_BICYCLE_WHEEL_MAX_CUT, 3>;
   using triplet_4_clique_message_013_container = MessageContainer<max_cut_triplet_4_clique_message_013, 1, 2, variableMessageNumber, 1, FMC_ODD_BICYCLE_WHEEL_MAX_CUT, 4>;
   using triplet_4_clique_message_023_container = MessageContainer<max_cut_triplet_4_clique_message_023, 1, 2, variableMessageNumber, 1, FMC_ODD_BICYCLE_WHEEL_MAX_CUT, 5>;
   using triplet_4_clique_message_123_container = MessageContainer<max_cut_triplet_4_clique_message_123, 1, 2, variableMessageNumber, 1, FMC_ODD_BICYCLE_WHEEL_MAX_CUT, 6>;

   using four_five_clique_message_0123_container = MessageContainer<max_cut_4_5_clique_message_0123, 2, 3, variableMessageNumber, 1, FMC_ODD_BICYCLE_WHEEL_MAX_CUT, 7>;
   using four_five_clique_message_0124_container = MessageContainer<max_cut_4_5_clique_message_0124, 2, 3, variableMessageNumber, 1, FMC_ODD_BICYCLE_WHEEL_MAX_CUT, 8>;
   using four_five_clique_message_0134_container = MessageContainer<max_cut_4_5_clique_message_0134, 2, 3, variableMessageNumber, 1, FMC_ODD_BICYCLE_WHEEL_MAX_CUT, 9>;
   using four_five_clique_message_0234_container = MessageContainer<max_cut_4_5_clique_message_0234, 2, 3, variableMessageNumber, 1, FMC_ODD_BICYCLE_WHEEL_MAX_CUT, 10>;
   using four_five_clique_message_1234_container = MessageContainer<max_cut_4_5_clique_message_1234, 2, 3, variableMessageNumber, 1, FMC_ODD_BICYCLE_WHEEL_MAX_CUT, 11>;

   using FactorList = meta::list< 
      edge_factor_container,
      triplet_factor_container,
      four_clique_factor_container,
      five_clique_factor_container,
      constant_factor_container 
         >;

   using MessageList = meta::list<
      edge_triplet_message_0_container,
      edge_triplet_message_1_container,
      edge_triplet_message_2_container, 

      triplet_4_clique_message_012_container,
      triplet_4_clique_message_013_container,
      triplet_4_clique_message_023_container,
      triplet_4_clique_message_123_container,

      four_five_clique_message_0123_container,
      four_five_clique_message_0124_container,
      four_five_clique_message_0134_container,
      four_five_clique_message_0234_container,
      four_five_clique_message_1234_container
      >;

   using max_cut_c = MulticutConstructor<FMC_ODD_BICYCLE_WHEEL_MAX_CUT,0,1,0,1,2,4,cut_type::maxcut>;
   using max_cut_ow = MulticutOddWheelConstructor<max_cut_c,2,3,4,5,6,cut_type::maxcut>;
   using max_cut_bow = multicut_odd_bicycle_wheel_constructor<max_cut_ow, 3, 7,8,9,10,11, cut_type::maxcut>;
   using ProblemDecompositionList = meta::list<max_cut_bow>;
};

namespace max_cut_dimacs_input {
   // the dimacs input actually specifies the qpbo graph. We must revert the graph building process to obtain the original max-cut problem

   using Parsing::opt_whitespace;
   using Parsing::mand_whitespace;
   using Parsing::positive_integer;
   using Parsing::real_number;

   struct init_line : pegtl::seq< opt_whitespace, pegtl::string<'p',' ','m','a','x'>, mand_whitespace, positive_integer, mand_whitespace, positive_integer, opt_whitespace > {};
   struct source_line : pegtl::seq< opt_whitespace, pegtl::string<'n',' ','1',' ','t'>, opt_whitespace > {};
   struct sink_line : pegtl::seq< opt_whitespace, pegtl::string<'n',' ','2',' ','s'>, opt_whitespace > {};
   struct pos_unary_potential_line : pegtl::seq< opt_whitespace, pegtl::string<'a',' ','1'>, mand_whitespace, positive_integer, mand_whitespace, positive_integer, opt_whitespace > {};
   struct neg_unary_potential_line : pegtl::seq< opt_whitespace, pegtl::string<'a'>, mand_whitespace, positive_integer, pegtl::string<'2'>, mand_whitespace, positive_integer, opt_whitespace > {};
   // actually, the next is two lines: the arc and the reverse arc. The input respects this convention
   //struct pairwise_potential_line :

      /*
   struct grammar : pegtl::must<
                    init_line, pegtl::eol,
                    numberOfVariables_line, pegtl::eol,
                    pegtl::star<edge_line, pegtl::eol>,
                    pegtl::opt<edge_line>,
                    pegtl::star<pegtl::sor<mand_whitespace, pegtl::eol>>,
                    pegtl::eof> {};
                    */

}

namespace max_cut_simple_text_format {

   using Parsing::opt_whitespace;
   using Parsing::mand_whitespace;
   using Parsing::positive_integer;
   using Parsing::real_number;

   // #nodes #edges
   struct init_line : pegtl::seq< opt_whitespace, positive_integer, opt_whitespace, positive_integer, opt_whitespace> {};
   // ${node 1} ${node 2} ${cost}
   struct edge_line : pegtl::seq< opt_whitespace, positive_integer, opt_whitespace, positive_integer, opt_whitespace, real_number, opt_whitespace> {};

   struct grammar : pegtl::must<
                    init_line, pegtl::eol,
                    pegtl::star<edge_line, pegtl::eolf>
                    > {};

   struct input {
      INDEX no_variables;
      std::vector<std::tuple<INDEX,INDEX,REAL>> edges;
   };

   template<typename Rule >
      struct action
      : pegtl::nothing< Rule > {};

   template<>
   struct action<init_line > {
      template<typename INPUT>
      static void apply(const INPUT & in, input& ipt)
      {
         std::istringstream iss(in.string());
         INDEX no_nodes; iss >> no_nodes;
         INDEX no_edges; iss >> no_edges;

         ipt.no_variables = no_nodes;
      }
   };

   template<>
   struct action<edge_line > {
      template<typename INPUT>
      static void apply(const INPUT & in, input& ipt)
      {
         std::istringstream iss(in.string());
         INDEX i; iss >> i;
         INDEX j; iss >> j;
         REAL cost; iss >> cost;
         assert(i <= j);
         assert(std::max(i,j) <= ipt.no_variables);

         ipt.edges.push_back(std::make_tuple(i,j,cost));
      }
   };

   template<typename CONSTRUCTOR>
   void construct_max_cut_problem(CONSTRUCTOR& c, const input& ipt)
   {
      for(const auto& edge : ipt.edges) {
         const INDEX i = std::get<0>(edge);
         const INDEX j = std::get<1>(edge);
         assert(i > 0 && j > 0 && i < j);
         const REAL cost = std::get<2>(edge);
         c.AddUnaryFactor(std::min(i,j)-1, std::max(i,j)-1, -cost); // max cut problems maximize objective
      }
   }

   template<typename CONSTRUCTOR>
   void construct_QUBO_problem(CONSTRUCTOR& c, const input& ipt)
   {
      // transform QUBO to max cut: add one additional node
      std::vector<REAL> additional_edges(ipt.no_variables, 0.0);
      REAL sum = 0.0;
      for(const auto& edge : ipt.edges) {
         const INDEX i = std::get<0>(edge);
         const INDEX j = std::get<1>(edge);
         assert(i > 0 && j > 0);
         const REAL cost = std::get<2>(edge);

         if(i == j) {
            additional_edges[i-1] += cost/2.0;
            sum += cost;
         } else {
            const REAL cost = std::get<2>(edge);
            c.AddUnaryFactor(std::min(i,j)-1, std::max(i,j)-1, -cost); // max cut problems maximize objective
            additional_edges[i-1] += cost;
            additional_edges[j-1] += cost;
            sum += 2*cost;
         }
      }
      for(INDEX i=0; i<ipt.no_variables; ++i) {
         c.AddUnaryFactor(i, ipt.no_variables, -additional_edges[i]);
      }
      c.AddToConstant(sum);
   } 

   template<typename SOLVER>
   bool ParseProblemMaxCut(const std::string filename, SOLVER& s)
   {
      input ipt;
      std::cout << "parsing " << filename << "\n";

      pegtl::file_parser problem(filename);

      auto& c = s.template GetProblemConstructor<0>();
      if(problem.parse< grammar, action >(ipt)) {
         construct_max_cut_problem(c, ipt);
         return true;
      } else { 
         return false;
      }
   } 

   template<typename SOLVER>
   bool ParseProblemQUBO(const std::string filename, SOLVER& s)
   {
      input ipt;
      std::cout << "parsing " << filename << "\n";

      pegtl::file_parser problem(filename);

      auto& c = s.template GetProblemConstructor<0>();
      if(problem.parse< grammar, action >(ipt)) {
         construct_QUBO_problem(c, ipt);
         return true;
      } else { 
         return false;
      }
   } 

}


} // end namespace LP_MP

#endif // LP_MP_MAX_CUT_H

