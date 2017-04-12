#ifndef LP_MP_TOMOGRAPHY_H
#define LP_MP_TOMOGRAPHY_H

#include <iostream>
#include <fstream>
#include "LP_MP.h"
#include "solver.hxx"
#include "factors/simplex_factor.hxx"
#include "messages/simplex_marginalization_message.hxx"
#include "problem_constructors/mrf_problem_construction.hxx"

#include "discrete_tomography_factor_counting.hxx"
#include "discrete_tomography_message_counting.hxx"
#include "discrete_tomography_message_counting_pairwise.hxx"
#include "discrete_tomography_tree_constructor.hxx"
#include "discrete_tomography_counting_naive.hxx"
#include "dt_sequential.hxx"

#include "parse_rules.h"

/*
  This file defines the FMC (factors,messages,...) for the tomography solver
*/

namespace LP_MP{

  struct FMC_DT {
    static constexpr char* name = "Discrete Tomography with MRF and Counting Factors";

    typedef FactorContainer<UnarySimplexFactor, FMC_DT, 0, false> UnaryFactor;
    typedef FactorContainer<PairwiseSimplexFactor, FMC_DT, 1> PairwiseFactor;
    typedef FactorContainer<DiscreteTomographyFactorCounting2, FMC_DT, 2> DiscreteTomographyCountingFactorContainer;
    //using dt_sequential_factor = FactorContainer<dt_sum_state_factor, FMC_DT, 3>;
    using dt_sequential_pairwise_factor = FactorContainer<dt_sum_state_pairwise_factor, FMC_DT, 3>;
   
    // do zrobienia: possibly try out MPLP as well
    typedef MessageContainer<UnaryPairwiseMessageLeft<MessageSendingType::SRMP>, 0, 1, variableMessageNumber, 1, FMC_DT, 0 > UnaryPairwiseMessageLeftContainer;
    typedef MessageContainer<UnaryPairwiseMessageRight<MessageSendingType::SRMP>, 0, 1, variableMessageNumber, 1, FMC_DT, 1 > UnaryPairwiseMessageRightContainer;
    
    typedef MessageContainer<DiscreteTomographyMessageCounting2, 2, 2, atMostTwoMessages, atMostTwoMessages, FMC_DT, 2>
      DiscreteTomographyCountingMessageLeft;
    typedef MessageContainer<DiscreteTomographyMessageCounting2, 2, 2, atMostOneMessage, atMostOneMessage, FMC_DT, 3>
      DiscreteTomographyCountingMessageRight;

    // connect counting factors with pairwise ones
    using pairwise_to_center_counting = MessageContainer<DiscreteTomographyMessageCountingPairwise2<CountingPairwiseMessageType::center>, 1, 2, variableMessageNumber, 1, FMC_DT, 4>;
    using pairwise_to_left_counting = MessageContainer<DiscreteTomographyMessageCountingPairwise2<CountingPairwiseMessageType::left>, 1, 2, variableMessageNumber, atMostOneMessage, FMC_DT, 5>;
    using pairwise_to_right_counting = MessageContainer<DiscreteTomographyMessageCountingPairwise2<CountingPairwiseMessageType::right>, 1, 2, variableMessageNumber, atMostOneMessage, FMC_DT, 6>;

    using left_pairwise_to_left_counting = MessageContainer<DiscreteTomographyMessageCountingPairwise3<Chirality::left,Chirality::left>, 1, 2, variableMessageNumber, atMostOneMessage, FMC_DT, 7>;
    using right_pairwise_to_left_counting = MessageContainer<DiscreteTomographyMessageCountingPairwise3<Chirality::left,Chirality::right>, 1, 2, variableMessageNumber, atMostOneMessage, FMC_DT, 8>;
    using left_pairwise_to_right_counting = MessageContainer<DiscreteTomographyMessageCountingPairwise3<Chirality::right,Chirality::left>, 1, 2, variableMessageNumber, atMostOneMessage, FMC_DT, 9>;
    using right_pairwise_to_right_counting = MessageContainer<DiscreteTomographyMessageCountingPairwise3<Chirality::right,Chirality::right>, 1, 2, variableMessageNumber, atMostOneMessage, FMC_DT, 10>;

    // connect sequential dt factors
    using dt_sequential_message = MessageContainer<dt_pairwise_message, 3, 3, atMostOneMessage, atMostOneMessage, FMC_DT, 11>;
    //using dt_sequential_message_left = MessageContainer<dt_sum_pairwise_message<Chirality::left>, 3,4, atMostOneMessage, atMostOneMessage, FMC_DT, 11>;
    //using dt_sequential_message_right = MessageContainer<dt_sum_pairwise_message<Chirality::right>, 3,4, atMostOneMessage, atMostOneMessage, FMC_DT, 12>;
    using dt_sequential_recursive_message_left = MessageContainer<dt_sum_counting_message<Chirality::left>, 3,2, atMostOneMessage, atMostOneMessage, FMC_DT, 12>;
    using dt_sequential_recursive_message_right = MessageContainer<dt_sum_counting_message<Chirality::right>, 3,2, atMostOneMessage, atMostOneMessage, FMC_DT, 13>;
    using dt_pairwise_pairwise_message = MessageContainer<dt_sum_pairwise_pairwise_message, 1, 3, variableMessageNumber, atMostOneMessage, FMC_DT, 14>;
    using dt_sum_unary_message_left = MessageContainer<dt_sum_unary_message<Chirality::left>, 0, 3, variableMessageNumber, 1, FMC_DT, 18>;
    using dt_sum_unary_message_right = MessageContainer<dt_sum_unary_message<Chirality::right>, 0, 3, variableMessageNumber, 1, FMC_DT, 19>;

    //using dt_unary_sum_message_container = MessageContainer<dt_unary_sum_message, 0, 3, variableMessageNumber, 1, FMC_DT, 19>;

    // tightening
    typedef FactorContainer<SimpleTighteningTernarySimplexFactor, FMC_DT, 4> EmptyTripletFactor;
    typedef MessageContainer<PairwiseTripletMessage<0,1,MessageSendingType::SRMP>, 1, 4, variableMessageNumber, 1, FMC_DT, 15> PairwiseTriplet12MessageContainer;
    typedef MessageContainer<PairwiseTripletMessage<0,2,MessageSendingType::SRMP>, 1, 4, variableMessageNumber, 1, FMC_DT, 16> PairwiseTriplet13MessageContainer;
    typedef MessageContainer<PairwiseTripletMessage<1,2,MessageSendingType::SRMP>, 1, 4, variableMessageNumber, 1, FMC_DT, 17> PairwiseTriplet23MessageContainer;

    using FactorList = meta::list< UnaryFactor, PairwiseFactor, DiscreteTomographyCountingFactorContainer, dt_sequential_pairwise_factor, EmptyTripletFactor >;
    using MessageList = meta::list<
       UnaryPairwiseMessageLeftContainer,
       UnaryPairwiseMessageRightContainer,
       DiscreteTomographyCountingMessageLeft,
       DiscreteTomographyCountingMessageRight,
       pairwise_to_center_counting,
       pairwise_to_left_counting,
       pairwise_to_right_counting,
       left_pairwise_to_left_counting,
       right_pairwise_to_left_counting,
       left_pairwise_to_right_counting,
       right_pairwise_to_right_counting,
       dt_sequential_message,
       dt_sequential_recursive_message_left,
       dt_sequential_recursive_message_right,
       dt_pairwise_pairwise_message,
       PairwiseTriplet12MessageContainer,
       PairwiseTriplet13MessageContainer,
       PairwiseTriplet23MessageContainer,
       dt_sum_unary_message_left,
       dt_sum_unary_message_right
       >;

    using mrf = StandardMrfConstructor<FMC_DT,0,1,0,1>;
    using tighteningMrf = TighteningMRFProblemConstructor<mrf,4,15,16,17>;
    using dt_recursive = DiscreteTomographyTreeConstructor<FMC_DT,0,
          DiscreteTomographyCountingFactorContainer,
          DiscreteTomographyCountingMessageLeft,
          DiscreteTomographyCountingMessageRight,
          pairwise_to_center_counting,
          pairwise_to_left_counting,
          pairwise_to_right_counting,
          left_pairwise_to_left_counting,
          right_pairwise_to_left_counting,
          left_pairwise_to_right_counting,
          right_pairwise_to_right_counting>;
    using dt_sequential = dt_sequential_constructor<FMC_DT, 0,
          dt_sequential_pairwise_factor, 
          dt_sequential_message,
          dt_pairwise_pairwise_message,
          dt_sum_unary_message_left,
          dt_sum_unary_message_right
          >;
    using dt_combined = dt_combined_constructor<FMC_DT, 0,
          dt_sequential,
          dt_recursive,
          dt_sequential_recursive_message_left,
          dt_sequential_recursive_message_right>;


    //using ProblemDecompositionList = meta::list<tighteningMrf,dt_sequential>;
    using ProblemDecompositionList = meta::list<tighteningMrf,dt_sequential>;
    //using ProblemDecompositionList = meta::list<tighteningMrf,dt_recursive>;
    //using ProblemDecompositionList = meta::list<tighteningMrf,dt_combined>;
	  
  };

  /*
  struct FMC_DT_NAIVE {
    static constexpr char* name = "Discrete Tomography, naive LP model";

    typedef FactorContainer<Simplex, ExplicitRepamStorage, FMC_DT_NAIVE, 0, false> UnaryFactor;
    typedef FactorContainer<Simplex, ExplicitRepamStorage, FMC_DT_NAIVE, 1> PairwiseFactor;
    typedef FactorContainer<DiscreteTomographyFactorCountingNaive, ExplicitRepamStorage, FMC_DT_NAIVE, 2> DiscreteTomographyCountingFactorContainer;
   
    typedef MessageContainer<LeftMargMessage, 0, 1, variableMessageNumber, 1, variableMessageSize, FMC_DT_NAIVE, 0 > UnaryPairwiseMessageLeft;
    typedef MessageContainer<RightMargMessage, 0, 1, variableMessageNumber, 1, variableMessageSize, FMC_DT_NAIVE, 1 > UnaryPairwiseMessageRight;
    
    typedef MessageContainer<DiscreteTomographyUnaryToFactorCountingNaiveMessage, 0, 2, variableMessageNumber, variableMessageNumber, variableMessageSize, FMC_DT_NAIVE, 2>
      DiscreteTomographyCountingMessage;

    using FactorList = meta::list< UnaryFactor, PairwiseFactor, DiscreteTomographyCountingFactorContainer >;
    using MessageList = meta::list<
      UnaryPairwiseMessageLeft,
      UnaryPairwiseMessageRight,
      DiscreteTomographyCountingMessage
      >;

    using mrf = StandardMrfConstructor<FMC_DT_NAIVE,0,1,0,1>;
    using dt = DiscreteTomographyNaiveConstructor<FMC_DT_NAIVE,0,0,2,2>;
    using ProblemDecompositionList = meta::list<mrf,dt>;
  };


  struct FMC_DT_COMBINED {
    static constexpr char* name = "Discrete Tomography, naive LP and tight model combined";

    typedef FactorContainer<Simplex, ExplicitRepamStorage, FMC_DT_COMBINED, 0, true> UnaryFactor;
    typedef FactorContainer<Simplex, ExplicitRepamStorage, FMC_DT_COMBINED, 1> PairwiseFactor;
    typedef FactorContainer<DiscreteTomographyFactorCounting, ExplicitRepamStorage, FMC_DT_COMBINED, 2> DiscreteTomographyCountingFactorContainer;
   
    typedef MessageContainer<LeftMargMessage, 0, 1, variableMessageNumber, 1, variableMessageSize, FMC_DT_COMBINED, 0 > UnaryPairwiseMessageLeft;
    typedef MessageContainer<RightMargMessage, 0, 1, variableMessageNumber, 1, variableMessageSize, FMC_DT_COMBINED, 1 > UnaryPairwiseMessageRight;
    
    typedef MessageContainer<DiscreteTomographyMessageCounting<DIRECTION::left>, 2, 2, atMostOneMessage, atMostOneMessage, variableMessageSize, FMC_DT_COMBINED, 2>
      DiscreteTomographyCountingMessageLeft;
    typedef MessageContainer<DiscreteTomographyMessageCounting<DIRECTION::right>, 2, 2, atMostOneMessage, atMostOneMessage, variableMessageSize, FMC_DT_COMBINED, 3>
      DiscreteTomographyCountingMessageRight;
    typedef MessageContainer<DiscreteTomographyMessageCountingPairwise, 1, 2, variableMessageNumber, 1, variableMessageSize, FMC_DT_COMBINED, 4>
      DiscreteTomographyCountingPairwiseMessageContainer;

   // tightening
   typedef FactorContainer<Simplex, ExplicitRepamStorage, FMC_DT_COMBINED, 3> EmptyTripletFactor;
   typedef MessageContainer<PairwiseTriplet12Message, 1, 3, variableMessageNumber, 1, variableMessageSize, FMC_DT_COMBINED, 5> PairwiseTriplet12MessageContainer;
   typedef MessageContainer<PairwiseTriplet13Message, 1, 3, variableMessageNumber, 1, variableMessageSize, FMC_DT_COMBINED, 6> PairwiseTriplet13MessageContainer;
   typedef MessageContainer<PairwiseTriplet23Message, 1, 3, variableMessageNumber, 1, variableMessageSize, FMC_DT_COMBINED, 7> PairwiseTriplet23MessageContainer;

   // naive factors
   typedef FactorContainer<DiscreteTomographyFactorCountingNaive, ExplicitRepamStorage, FMC_DT_COMBINED, 4> DiscreteTomographyCountingFactorNaiveContainer;
   typedef MessageContainer<DiscreteTomographyUnaryToFactorCountingNaiveMessage, 0, 4, variableMessageNumber, variableMessageNumber, variableMessageSize, FMC_DT_COMBINED, 8>
      DiscreteTomographyCountingNaiveMessage;

    using FactorList = meta::list< UnaryFactor, PairwiseFactor, DiscreteTomographyCountingFactorContainer, EmptyTripletFactor, DiscreteTomographyCountingFactorNaiveContainer >;
    using MessageList = meta::list<
      UnaryPairwiseMessageLeft,
      UnaryPairwiseMessageRight,
      DiscreteTomographyCountingMessageLeft,
      DiscreteTomographyCountingMessageRight,
      DiscreteTomographyCountingPairwiseMessageContainer,
      PairwiseTriplet12MessageContainer,
      PairwiseTriplet13MessageContainer,
      PairwiseTriplet23MessageContainer,
      DiscreteTomographyCountingNaiveMessage
      >;

    using mrf = StandardMrfConstructor<FMC_DT_COMBINED,0,1,0,1>;
    using tighteningMrf = TighteningMRFProblemConstructor<mrf,3,5,6,7>;
    using dt = DiscreteTomographyTreeConstructor<FMC_DT_COMBINED,0,2,2,3,4>;
    using dt_naive = DiscreteTomographyNaiveConstructor<FMC_DT_COMBINED,0,0,4,8>;
    using ProblemDecompositionList = meta::list<tighteningMrf,dt,dt_naive>;
	  
  };
  */

  namespace DiscreteTomographyTextInput {

     using Parsing::opt_whitespace;
     using Parsing::mand_whitespace;
     using Parsing::positive_integer;
     using Parsing::real_number;

    struct ProjectionPreamble : pegtl::string<'P','R','O','J','E','C','T','I','O','N','S'> {};
    struct ProjectionVector : pegtl::seq< pegtl::string<'('>, opt_whitespace, real_number, opt_whitespace, pegtl::star< pegtl::string<','>, opt_whitespace, real_number, opt_whitespace >, opt_whitespace, pegtl::string<')'> > {};
    struct first_positive_integer : positive_integer {};
    struct ProjectionLine : pegtl::seq<first_positive_integer,pegtl::plus<opt_whitespace,pegtl::string<'+'>,opt_whitespace,positive_integer>,opt_whitespace,pegtl::string<'='>,opt_whitespace,ProjectionVector> {};

   
    // projection grammar
    struct grammar : pegtl::seq<
      pegtl::until<ProjectionPreamble>,
      pegtl::star<pegtl::sor<pegtl::seq<opt_whitespace,pegtl::eol>,ProjectionLine>>,
      pegtl::eof> {};

    struct Projections {
       std::vector<std::vector<INDEX>> projectionVar;
       std::vector<std::vector<REAL>> projectionCost;
    };

    template<typename Rule>
      struct action : pegtl::nothing<Rule> {};
    
    template<>
       struct action<pegtl::string<'('>> {
       template<typename INPUT>
          static void apply(const INPUT& in, Projections& p)
          {
             p.projectionCost.push_back({});
          }
       };

    template<>
       struct action<real_number> {
       template<typename INPUT>
          static void apply(const INPUT& in, Projections& p)
          {
             p.projectionCost.back().push_back(std::stod(in.string()));
             //realStack.push(std::stod(in.string()));
          }
       };

    template<>
       struct action<first_positive_integer> {
       template<typename INPUT>
          static void apply(const INPUT& in, Projections& p)
          {
             p.projectionVar.push_back({});
             p.projectionVar.back().push_back(std::stoul(in.string()));
          }
       };

    template<>
       struct action<positive_integer> {
       template<typename INPUT>
          static void apply(const INPUT& in, Projections& p)
          {
             p.projectionVar.back().push_back(std::stoul(in.string()));
             //integerStack.push(std::stoul(in.string()));
          }
       };

    template<typename SOLVER>
    inline bool ParseProblem(const std::string& filename, SOLVER& pd) {
       std::cout << "parsing " << filename << "\n";
       pegtl::file_parser problem(filename);

       UaiMrfInput::MrfInput mrfInput;
       bool ret = problem.parse< UaiMrfInput::grammar, UaiMrfInput::action>(mrfInput);
       if(ret != true) {
          throw std::runtime_error("could not read mrf problem in uai format for discrete tomography");
       }
       UaiMrfInput::build_mrf(pd.template GetProblemConstructor<0>(), mrfInput);

       pd.template GetProblemConstructor<1>().SetNumberOfLabels(mrfInput.cardinality_[0]);

       Projections p;
       ret = problem.parse< grammar, action>(p);
       if(ret != true) {
          throw std::runtime_error("could not read projection constraints for discrete tomography");
       }
       assert(p.projectionVar.size() == p.projectionCost.size());
       for(INDEX i=0; i<p.projectionVar.size(); ++i) {
          pd.template GetProblemConstructor<1>().AddProjection(p.projectionVar[i], p.projectionCost[i]);
          //LP_tree t;
          //pd.template GetProblemConstructor<1>().AddProjection(p.projectionVar[i], p.projectionCost[i], &t);
          //t.compute_subgradient();
          //assert(false); 
       }
       return true;
    }

    template<typename SOLVER>
    inline bool ParseProblemCombined(const std::string& filename, SOLVER& pd) {
       std::cout << "parsing " << filename << "\n";
       pegtl::file_parser problem(filename);

       UaiMrfInput::MrfInput mrfInput;
       bool ret = problem.parse< UaiMrfInput::grammar, UaiMrfInput::action>(mrfInput);
       if(ret != true) {
          throw std::runtime_error("could not read mrf problem in uai format for discrete tomography");
       }
       UaiMrfInput::build_mrf(pd.template GetProblemConstructor<0>(), mrfInput);

       pd.template GetProblemConstructor<1>().SetNumberOfLabels(mrfInput.cardinality_[0]);
       pd.template GetProblemConstructor<2>().SetNumberOfLabels(mrfInput.cardinality_[0]);

       Projections p;
       ret = problem.parse< grammar, action>(p);
       if(ret != true) {
          throw std::runtime_error("could not read projection constraints for discrete tomography");
       }
       assert(p.projectionVar.size() == p.projectionCost.size());
       for(INDEX i=0; i<p.projectionVar.size(); ++i) {
          pd.template GetProblemConstructor<1>().AddProjection(p.projectionVar[i], p.projectionCost[i]);
       }
       for(INDEX i=0; i<p.projectionVar.size(); ++i) {
          pd.template GetProblemConstructor<2>().AddProjection(p.projectionVar[i], p.projectionCost[i]);
       }
       return true;
    }

  }

} // end namespace LP_MP
#endif // LP_MP_TOMOGRAPHY_H
