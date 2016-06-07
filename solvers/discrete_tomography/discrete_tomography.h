#ifndef LP_MP_TOMOGRAPHY_H
#define LP_MP_TOMOGRAPHY_H

#include "parse_rules.h"

/*
  This file defines the FMC (factors,messages,...) for the tomography solver
*/



namespace LP_MP{

  typedef UnaryLoop<> UnaryLoopType;
  typedef PairwiseLoop<0> LeftLoopType;
  typedef PairwiseLoop<1> RightLoopType;

  typedef MultiplexMargMessage<UnaryLoopType,LeftLoopType,true,false,false,true> LeftMargMessage;
  typedef MultiplexMargMessage<UnaryLoopType,RightLoopType,true,false,false,true> RightMargMessage;

  typedef MultiplexFactor<std::vector<REAL>, const_ones_array, const_one> Simplex;


  struct FMC_DT {

    typedef FactorContainer<Simplex, ExplicitRepamStorage, FMC_DT, 0, true, true > UnaryFactor;
    typedef FactorContainer<Simplex, ExplicitRepamStorage, FMC_DT, 1, false, false > PairwiseFactor;
    typedef FactorContainer<DiscreteTomographyCountingFactor, ExplicitRepamStorage, FMC_DT, 2, false, false> DiscreteTomographyCountingFactorContainer;
   
    typedef MessageContainer<LeftMargMessage, StandardMessageStorage, FMC_DT, 1 > UnaryPairwiseMessageLeft;
    typedef MessageContainer<RightMargMessage, StandardMessageStorage, FMC_DT, 2 > UnaryPairwiseMessageRight;
    typedef MessageContainer<DiscreteTomographyCountingMessage<Direction::Left>, FMC_DT, 3> DiscreteTomographyCountingMessageLeft;
    typedef MessageContainer<DiscreteTomographyCountingMessage<Direction::Right>, FMC_DT, 4> DiscreteTomographyCountingMessageRight;
    typedef MessageContainer<DiscreteTomographyCountingPairwiseMessage, FMC_DT, 5> DiscreteTomographyCountingCountingPairwiseMessageContainer;


    using FactorList = meta::list< UnaryFactor, PairwiseFactor, DiscreteTomographyCountingFactorContainer >;
    using MessageList = meta::list<
      MessageListItem< UnaryPairwiseMessageLeft,  0, 1, std::vector, FixedSizeMessageContainer<1>::type >,
      MessageListItem< UnaryPairwiseMessageRight, 0, 1, std::vector, FixedSizeMessageContainer<1>::type >,
      MessageListItem< DiscreteTomographyCountingMessageLeft, 2, 2,  std::vector, std::vector >,
      MessageListItem< DiscreteTomographyCountingMessageRight, 2, 2,  std::vector, std::vector >,
      MessageListItem< DiscreteTomographyCountingPairwiseMessageContainer, 1, 2,  std::vector, FixedSizeMessageContainer<1>::type >
      >;

    using mrf = MRFProblemConstructor<FMC_DT,0,1,0,1>;
    using dt = DiscreteTomographyConstructor<FMC_DT,mrf,2,3,4,5>;
    using ProblemDecompositionList = meta::list<mrf,dt>;
	  
  };

  namespace DiscreteTomographyTextInput {


    struct ProjectionPreamble : pegtl::string<'P','R','O','J','E','C','T','I','O','N','S'> {};
    struct ProjectionVector : pegtl::seq<pegtl::string<'('>,opt_whitespace, real_number, opt_whitespace, pegtl::star< pegtl::string<','>,opt_whitespace, real_number,opt_whitespace >, opt_whitespace, pegtl::string<')'>
      struct ProjectionLine : pegtl::seq<opt_whitespace,pegtl::plus<positive_integer,opt_whitespace,pegtl::string<'+'>,opt_whitespace>,positive_integer,opt_whitespace,pegtl::string<'='>,opt_whitespace,ProjectionsVector,opt_whitespace> {};

    struct grammar : pegtl::seq< ....
      pegtl::star<opt_whitespace,pegtl::eol>,
      ProjectionPreamble,
      pegtl::star<opt_whitespace,pegtl::eol,ProjectionLine>,
      pegtl::eof> {};

    struct action<real_number> {
      static void apply(const pegtl::input& in, ProblemDecomposition<FMC_DT>& pd, std::stack<SIGNED_INDEX>& integerStack, std::stack<REAL>& realStack)
      {
	realStack.push(std::stod(in.string()));
      }
    };
    struct action<positive_integer> {
      static void apply(const pegtl::input& in, ProblemDecomposition<FMC_DT>& pd, std::stack<SIGNED_INDEX>& integerStack, std::stack<REAL>& realStack)
      {
	integerStack.push(std::stoi(in.string()));
      }
    };
    
    struct action<ProjectionLine> {
      static void apply(const pegtl::input& in, ProblemDecomposition<FMC_DT>& pd, std::stack<SIGNED_INDEX>& integerStack, std::stack<REAL>& realStack)
      {
	std::vector<INDEX> projectionVar;
	while(!integerStack.empty) {
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
	pd.GetProblemConstructor<1>().AddProjection(projectionVar,projectionCost);
      }
    };

  }
#endif // LP_MP_TOMOGRAPHY_H
