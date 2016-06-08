#ifndef LP_MP_TOMOGRAPHY_H
#define LP_MP_TOMOGRAPHY_H

#include <iostream>
#include <fstream>
#include "LP_MP.h"
#include "factors/simplex_factor.hxx"
#include "messages/simplex_marginalization_message.hxx"
#include "problem_constructors/mrf_problem_construction.hxx"

#include "discrete_tomography_factor_counting.hxx"
#include "discrete_tomography_message_counting.hxx"
#include "discrete_tomography_message_counting_pairwise.hxx"
#include "discrete_tomography_tree_constructor.hxx"

#include "parse_rules.h"

/*
  This file defines the FMC (factors,messages,...) for the tomography solver
*/

namespace LP_MP{

  typedef UnaryLoop<> UnaryLoopType;
  typedef PairwiseLoop<0> LeftLoopType;
  typedef PairwiseLoop<1> RightLoopType;

  typedef SimplexMarginalizationMessage<UnaryLoopType,LeftLoopType,true,false,false,true> LeftMargMessage;
  typedef SimplexMarginalizationMessage<UnaryLoopType,RightLoopType,true,false,false,true> RightMargMessage;

  typedef SimplexFactor<> Simplex;


  struct FMC_DT {

    typedef FactorContainer<Simplex, ExplicitRepamStorage, FMC_DT, 0, true, true > UnaryFactor;
    typedef FactorContainer<Simplex, ExplicitRepamStorage, FMC_DT, 1, false, false > PairwiseFactor;
    typedef FactorContainer<DiscreteTomographyFactorCounting, ExplicitRepamStorage, FMC_DT, 2, false, false> DiscreteTomographyCountingFactorContainer;
   
    typedef MessageContainer<LeftMargMessage, 0, 1, variableMessageNumber, 1, variableMessageSize, FMC_DT, 0 > UnaryPairwiseMessageLeft;
    typedef MessageContainer<RightMargMessage, 0, 1, variableMessageNumber, 1, variableMessageSize, FMC_DT, 1 > UnaryPairwiseMessageRight;
    
    typedef MessageContainer<DiscreteTomographyMessageCounting<DIRECTION::left>, 2, 2, variableMessageNumber, variableMessageNumber, variableMessageSize, FMC_DT, 2>
      DiscreteTomographyCountingMessageLeft;
    typedef MessageContainer<DiscreteTomographyMessageCounting<DIRECTION::right>, 2, 2, variableMessageNumber, variableMessageNumber, variableMessageSize, FMC_DT, 3>
      DiscreteTomographyCountingMessageRight;
    typedef MessageContainer<DiscreteTomographyMessageCountingPairwise, 1, 2, variableMessageNumber, 1, variableMessageSize, FMC_DT, 4>
      DiscreteTomographyCountingPairwiseMessageContainer;

    using FactorList = meta::list< UnaryFactor, PairwiseFactor, DiscreteTomographyCountingFactorContainer >;
    using MessageList = meta::list<
      UnaryPairwiseMessageLeft,
      UnaryPairwiseMessageRight,
      DiscreteTomographyCountingMessageLeft,
      DiscreteTomographyCountingMessageRight,
      DiscreteTomographyCountingPairwiseMessageContainer>;

    using mrf = StandardMrfConstructor<FMC_DT,0,1,0,1>;
    using dt = DiscreteTomographyTreeConstructor<FMC_DT,0,2,3,4,5>;
    using ProblemDecompositionList = meta::list<mrf,dt>;
	  
  };

  namespace DiscreteTomographyTextInput {


    struct ProjectionPreamble : pegtl::string<'P','R','O','J','E','C','T','I','O','N','S'> {};
    struct ProjectionVector : pegtl::seq< pegtl::string<'('>, opt_whitespace, real_number, opt_whitespace, pegtl::star< pegtl::string<','>, opt_whitespace, real_number,opt_whitespace >, opt_whitespace, pegtl::string<')'> > {};
    struct ProjectionLine : pegtl::seq<opt_whitespace,pegtl::plus<positive_integer,opt_whitespace,pegtl::string<'+'>,opt_whitespace>,positive_integer,opt_whitespace,pegtl::string<'='>,opt_whitespace,ProjectionVector,opt_whitespace> {};

    struct grammar : pegtl::seq<
      // add mrf grammar here
      pegtl::star<opt_whitespace,pegtl::eol>,
      ProjectionPreamble,
      pegtl::star<opt_whitespace,pegtl::eol,ProjectionLine>,
      pegtl::eof> {};

    template<typename Rule>
      struct action : pegtl::nothing<Rule> {};

    template<>
    struct action<real_number> {
      static void apply(const pegtl::input& in, ProblemDecomposition<FMC_DT>& pd, std::stack<SIGNED_INDEX>& integerStack, std::stack<REAL>& realStack)
      {
	realStack.push(std::stod(in.string()));
      }
    };
    template<>
    struct action<positive_integer> {
      static void apply(const pegtl::input& in, ProblemDecomposition<FMC_DT>& pd, std::stack<SIGNED_INDEX>& integerStack, std::stack<REAL>& realStack)
      {
	integerStack.push(std::stoi(in.string()));
      }
    };

    template<>
    struct action<ProjectionLine> {
      static void apply(const pegtl::input& in, ProblemDecomposition<FMC_DT>& pd, std::stack<SIGNED_INDEX>& integerStack, std::stack<REAL>& realStack)
      {
	std::vector<INDEX> projectionVar;
	while(!integerStack.empty()) {
	  projectionVar.push_back(integerStack.top());
	  integerStack.pop();
	}
	std::reverse(projectionVar.begin(),projectionVar.end());

	std::vector<REAL> projectionCost;
	while(!realStack.empty()) {
	  projectionCost.push_back(realStack.top());
	  realStack.pop(); 
	}
	std::reverse(projectionCost.begin(), projectionCost.end());
	auto& dt = pd.template GetProblemConstructor<1>();
	//.AddProjection(projectionVar,projectionCost);
      }
    };

    bool ParseLiftedProblem(const std::string filename&, ProblemDecomposition<FMC_DT>& pd) {
      return true;
    }

  }
#endif // LP_MP_TOMOGRAPHY_H
