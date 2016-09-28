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

#include "parse_rules.h"

/*
  This file defines the FMC (factors,messages,...) for the tomography solver
*/

namespace LP_MP{

  typedef UnaryLoop<> UnaryLoopType;
  typedef PairwiseLoop<0> LeftLoopType;
  typedef PairwiseLoop<1> RightLoopType;

  //typedef SimplexMarginalizationMessage<UnaryLoopType,LeftLoopType,true,true,false,true> LeftMargMessage;
  //typedef SimplexMarginalizationMessage<UnaryLoopType,RightLoopType,true,true,false,true> RightMargMessage;
  typedef SimplexMarginalizationMessage<UnaryLoopType,LeftLoopType,true,false,false,true> LeftMargMessage;
  typedef SimplexMarginalizationMessage<UnaryLoopType,RightLoopType,true,false,false,true> RightMargMessage;

  typedef SimplexFactor<> Simplex;

  // messages to triplet factors for tightening
  typedef PairwiseTripletLoop<0,1> PairwiseTripletLoopType12;
  typedef PairwiseTripletLoop<0,2> PairwiseTripletLoopType13;
  typedef PairwiseTripletLoop<1,2> PairwiseTripletLoopType23;
  typedef SimplexMarginalizationMessage<UnaryLoopType,PairwiseTripletLoopType12,true,false,false,true> PairwiseTriplet12Message;
  typedef SimplexMarginalizationMessage<UnaryLoopType,PairwiseTripletLoopType13,true,false,false,true> PairwiseTriplet13Message;
  typedef SimplexMarginalizationMessage<UnaryLoopType,PairwiseTripletLoopType23,true,false,false,true> PairwiseTriplet23Message;


  struct FMC_DT {
    static constexpr char* name = "Discrete Tomography";

    typedef FactorContainer<Simplex, ExplicitRepamStorage, FMC_DT, 0, true, true> UnaryFactor;
    typedef FactorContainer<Simplex, ExplicitRepamStorage, FMC_DT, 1> PairwiseFactor;
    typedef FactorContainer<DiscreteTomographyFactorCounting, ExplicitRepamStorage, FMC_DT, 2> DiscreteTomographyCountingFactorContainer;
   
    typedef MessageContainer<LeftMargMessage, 0, 1, variableMessageNumber, 1, variableMessageSize, FMC_DT, 0 > UnaryPairwiseMessageLeft;
    typedef MessageContainer<RightMargMessage, 0, 1, variableMessageNumber, 1, variableMessageSize, FMC_DT, 1 > UnaryPairwiseMessageRight;
    
    typedef MessageContainer<DiscreteTomographyMessageCounting<DIRECTION::left>, 2, 2, atMostOneMessage, atMostOneMessage, variableMessageSize, FMC_DT, 2>
      DiscreteTomographyCountingMessageLeft;
    typedef MessageContainer<DiscreteTomographyMessageCounting<DIRECTION::right>, 2, 2, atMostOneMessage, atMostOneMessage, variableMessageSize, FMC_DT, 3>
      DiscreteTomographyCountingMessageRight;
    typedef MessageContainer<DiscreteTomographyMessageCountingPairwise, 1, 2, variableMessageNumber, 1, variableMessageSize, FMC_DT, 4>
      DiscreteTomographyCountingPairwiseMessageContainer;

   // tightening
   typedef FactorContainer<Simplex, ExplicitRepamStorage, FMC_DT, 3> EmptyTripletFactor;
   typedef MessageContainer<PairwiseTriplet12Message, 1, 3, variableMessageNumber, 1, variableMessageSize, FMC_DT, 5> PairwiseTriplet12MessageContainer;
   typedef MessageContainer<PairwiseTriplet13Message, 1, 3, variableMessageNumber, 1, variableMessageSize, FMC_DT, 6> PairwiseTriplet13MessageContainer;
   typedef MessageContainer<PairwiseTriplet23Message, 1, 3, variableMessageNumber, 1, variableMessageSize, FMC_DT, 7> PairwiseTriplet23MessageContainer;

    using FactorList = meta::list< UnaryFactor, PairwiseFactor, DiscreteTomographyCountingFactorContainer, EmptyTripletFactor >;
    using MessageList = meta::list<
      UnaryPairwiseMessageLeft,
      UnaryPairwiseMessageRight,
      DiscreteTomographyCountingMessageLeft,
      DiscreteTomographyCountingMessageRight,
      DiscreteTomographyCountingPairwiseMessageContainer,
      PairwiseTriplet12MessageContainer,
      PairwiseTriplet13MessageContainer,
      PairwiseTriplet23MessageContainer
      >;

    using mrf = StandardMrfConstructor<FMC_DT,0,1,0,1>;
    using tighteningMrf = TighteningMRFProblemConstructor<mrf,3,5,6,7>;
    using dt = DiscreteTomographyTreeConstructor<FMC_DT,0,2,2,3,4>;
    using ProblemDecompositionList = meta::list<tighteningMrf,dt>;
	  
  };

  struct FMC_DT_NAIVE {
    static constexpr char* name = "Discrete Tomography, naive LP model";

    typedef FactorContainer<Simplex, ExplicitRepamStorage, FMC_DT_NAIVE, 0, true, true> UnaryFactor;
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
          static void apply(const pegtl::action_input& in, Projections& p)
          {
             p.projectionCost.push_back({});
          }
       };

    template<>
       struct action<real_number> {
          static void apply(const pegtl::action_input& in, Projections& p)
          {
             p.projectionCost.back().push_back(std::stod(in.string()));
             //realStack.push(std::stod(in.string()));
          }
       };

    template<>
       struct action<first_positive_integer> {
          static void apply(const pegtl::action_input& in, Projections& p)
          {
             p.projectionVar.push_back({});
             p.projectionVar.back().push_back(std::stoul(in.string()));
          }
       };

    template<>
       struct action<positive_integer> {
          static void apply(const pegtl::action_input& in, Projections& p)
          {
             p.projectionVar.back().push_back(std::stoul(in.string()));
             //integerStack.push(std::stoul(in.string()));
          }
       };

    template<typename FMC>
    inline bool ParseProblem(const std::string& filename, Solver<FMC>& pd) {
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
       }
       return true;
    }

  }

} // end namespace LP_MP
#endif // LP_MP_TOMOGRAPHY_H
