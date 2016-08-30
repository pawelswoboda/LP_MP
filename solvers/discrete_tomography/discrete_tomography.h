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

  namespace DiscreteTomographyTextInput {

    struct init_line : pegtl::seq< opt_whitespace, pegtl::string<'M','A','R','K','O','V'>, opt_whitespace > {};
    struct numberOfVariables_line : pegtl::seq< opt_whitespace, positive_integer, opt_whitespace > {};
    // vector of integers denoting how many labels each variable has
    struct cardinality_line : pegtl::seq< opt_whitespace, positive_integer, pegtl::star< opt_whitespace, positive_integer>, opt_whitespace> {};
    struct numberOfCliques : pegtl::seq< opt_whitespace, positive_integer, opt_whitespace> {};
    // first is the number of variables in the clique, then the actual variables.
    struct cliqueScope_line : pegtl::seq< opt_whitespace, positive_integer, pegtl::plus< opt_whitespace, positive_integer>, opt_whitespace, pegtl::eol> {};
    struct cliqueScopes : pegtl::star<cliqueScope_line> {};
    // a function table is begun by number of entries and then a list of real numbers. Here we record all the values in the real stack
    struct functionTables : pegtl::seq<pegtl::star<pegtl::sor<mand_whitespace, pegtl::eol>>, real_number, pegtl::plus<pegtl::star<pegtl::sor<mand_whitespace, pegtl::eol>>, real_number>> {};

    struct ProjectionPreamble : pegtl::string<'P','R','O','J','E','C','T','I','O','N','S'> {};
    struct ProjectionVector : pegtl::seq< pegtl::string<'('>, opt_whitespace, real_number, opt_whitespace, pegtl::star< pegtl::string<','>, opt_whitespace, real_number, opt_whitespace >, opt_whitespace, pegtl::string<')'> > {};
    struct ProjectionLine : pegtl::seq<positive_integer,pegtl::plus<opt_whitespace,pegtl::string<'+'>,opt_whitespace,positive_integer>,opt_whitespace,pegtl::string<'='>,opt_whitespace,ProjectionVector> {};

   
    struct grammar : pegtl::seq<
      // mrf grammar
      init_line, pegtl::eol,
      numberOfVariables_line, pegtl::eol,
      cardinality_line, pegtl::eol,
      numberOfCliques, pegtl::eol,
      cliqueScopes,
      functionTables,
      pegtl::star<pegtl::sor<mand_whitespace,pegtl::eol>>,
      pegtl::star<opt_whitespace,pegtl::eol>,
      // projection grammar
      ProjectionPreamble,
      pegtl::star<pegtl::sor<pegtl::seq<opt_whitespace,pegtl::eol>,ProjectionLine>>,
      pegtl::eof> {};

    struct MRFInput {
      INDEX numberOfVariables_;
      INDEX numberOfCliques_;
      std::vector<INDEX> cardinality_;
      std::vector<std::vector<INDEX>> cliqueScope_;
      std::vector<std::vector<REAL>> functionTable_;
    };
    
    template<typename Rule>
      struct action : pegtl::nothing<Rule> {};


    
    template<>
       struct action<real_number> {
          static void apply(const pegtl::input& in, Solver<FMC_DT>& pd, std::stack<SIGNED_INDEX>& integerStack, std::stack<REAL>& realStack, MRFInput& mrfInput)
          {
             realStack.push(std::stod(in.string()));
          }
       };
    template<>
       struct action<positive_integer> {
          static void apply(const pegtl::input& in, Solver<FMC_DT>& pd, std::stack<SIGNED_INDEX>& integerStack, std::stack<REAL>& realStack, MRFInput& mrfInput)
          {
             integerStack.push(std::stoul(in.string()));
          }
       };


    template<>
       struct action< numberOfVariables_line > {
          static void apply(const pegtl::input & in, Solver<FMC_DT>&, std::stack<SIGNED_INDEX>& integer_stack, std::stack<REAL>& real_stack, MRFInput& mrfInput)
          {
             mrfInput.numberOfVariables_ = integer_stack.top();
             integer_stack.pop();
          }
       };
    template<>
       struct action< numberOfCliques > {
          static void apply(const pegtl::input & in, Solver<FMC_DT>&, std::stack<SIGNED_INDEX>& integer_stack, std::stack<REAL>& real_stack, MRFInput& mrfInput)
          {
             mrfInput.numberOfCliques_ = integer_stack.top();
             integer_stack.pop();
          }
       };
    template<>
       struct action< cardinality_line > {
          static void apply(const pegtl::input & in, Solver<FMC_DT>&, std::stack<SIGNED_INDEX>& integer_stack, std::stack<REAL>& real_stack, MRFInput& mrfInput)
          {
             assert(integer_stack.size() ==mrfInput.numberOfVariables_);
             while(!integer_stack.empty()) {
                mrfInput.cardinality_.push_back(integer_stack.top());
                integer_stack.pop();
             }
             std::reverse(mrfInput.cardinality_.begin(), mrfInput.cardinality_.end());
             assert(real_stack.empty());
             assert(integer_stack.empty());
          }
       };
    template<>
       struct action< cliqueScope_line > {
          static void apply(const pegtl::input & in, Solver<FMC_DT>&, std::stack<SIGNED_INDEX>& integer_stack, std::stack<REAL>& real_stack, MRFInput& mrfInput)
          {
             mrfInput.cliqueScope_.push_back(std::vector<INDEX>(0));
             while(integer_stack.size() > 1) {
                mrfInput.cliqueScope_.back().push_back(integer_stack.top());
                integer_stack.pop();
             }
             std::reverse(mrfInput.cliqueScope_.back().begin(), mrfInput.cliqueScope_.back().end());
             const INDEX cliqueSize = integer_stack.top();
             integer_stack.pop();
             assert(mrfInput.cliqueScope_.back().size() == cliqueSize);
          }
       };
    // reconstruct the function tables
    template<>
       struct action< functionTables > {
          static void apply(const pegtl::input & in, Solver<FMC_DT>& pd, std::stack<SIGNED_INDEX>& integer_stack, std::stack<REAL>& real_stack, MRFInput& mrfInput)
          {
             // first read in real stack into array
             std::vector<REAL> values;
             while(!real_stack.empty()) {
                values.push_back(real_stack.top());
                real_stack.pop();
             }
             std::reverse(values.begin(), values.end());

             // now iterate over function tables. first read in number of values of function table, then add them to the function table.
             INDEX functionTableIdx = 0; // at which index does current function table start?
             while( functionTableIdx<values.size() ) {
                // assert that this functionTableSize is actually a real variable
                const INDEX functionTableSize = values[functionTableIdx];
                assert(functionTableIdx+functionTableSize < values.size());
                const INDEX potentialNumber = mrfInput.functionTable_.size();
                const INDEX cardinality = mrfInput.cliqueScope_[potentialNumber].size();

                mrfInput.functionTable_.push_back(std::vector<REAL>(functionTableSize,0.0));
                if(cardinality == 1) { // unary factor
                   assert(functionTableSize == mrfInput.cardinality_[ mrfInput.cliqueScope_[potentialNumber][0] ]);
                   for(INDEX label=0; label<functionTableSize; ++label) {
                      mrfInput.functionTable_.back().operator[](label) = values[functionTableIdx+1+label];
                      //std::cout << values[functionTableIdx+1+label] << ", ";
                   }
                   //std::cout << "\n";
                } else if(cardinality == 2) { // pairwise factor
                   const INDEX var1 = mrfInput.cliqueScope_[potentialNumber][0];
                   const INDEX var2 = mrfInput.cliqueScope_[potentialNumber][1];
                   assert(functionTableSize == mrfInput.cardinality_[var1] * mrfInput.cardinality_[var2]);
                   for(INDEX label1=0; label1<mrfInput.cardinality_[var1]; ++label1) {
                      for(INDEX label2=0; label2<mrfInput.cardinality_[var2]; ++label2) {
                         // note: we must transpose the matrix, that we have read in
                         mrfInput.functionTable_.back().operator[](label1 + label2*mrfInput.cardinality_[var1]) = values[functionTableIdx + 1 + label2 + label1*mrfInput.cardinality_[var2]];
                      }
                   }
                } else {
                   assert(false);
                   throw std::runtime_error("Only unary and pairwise potentials supported now");
                }
                functionTableIdx += functionTableSize+1;
             }

             // construct mrf
             auto& mrf = pd.template GetProblemConstructor<0>();

             // first input the unaries, as pairwise potentials need them to be able to link to them
             for(INDEX i=0; i<mrfInput.numberOfCliques_; ++i) {
                if(mrfInput.cliqueScope_[i].size() == 1) {
                   const INDEX var = mrfInput.cliqueScope_[i][0];
                   assert(mrfInput.functionTable_[i].size() == mrfInput.cardinality_[var]);
                   mrf.AddUnaryFactor(var,mrfInput.functionTable_[i]);
                } else if(mrfInput.cliqueScope_[i].size() > 2) {
                   throw std::runtime_error("only pairwise models are accepted currently");
                }
             }
             // now the pairwise potentials. 
             for(INDEX i=0; i<mrfInput.numberOfCliques_; ++i) {
                if(mrfInput.cliqueScope_[i].size() == 2) {
                   const INDEX var1 = mrfInput.cliqueScope_[i][0];
                   const INDEX var2 = mrfInput.cliqueScope_[i][1];
                   assert(var1<var2);
                   assert(mrfInput.functionTable_[i].size() == mrfInput.cardinality_[var1]*mrfInput.cardinality_[var2]);
                   mrf.AddPairwiseFactor(var1,var2,mrfInput.functionTable_[i]); // or do we have to transpose the values?
                }
             }
             assert(integer_stack.empty());
             assert(real_stack.empty());
          }
       };

    template<>
       struct action<ProjectionPreamble> {
          static void apply(const pegtl::input& in, Solver<FMC_DT>& pd, std::stack<SIGNED_INDEX>& integerStack, std::stack<REAL>& realStack, MRFInput& mrfInput)
          {
             for(INDEX i=0; i<mrfInput.numberOfVariables_-1; ++i) {
                assert(mrfInput.cardinality_[i] == mrfInput.cardinality_[i+1]);
             }
             pd.template GetProblemConstructor<1>().SetNumberOfLabels(mrfInput.cardinality_[0]);
             assert(integerStack.empty());
             assert(realStack.empty());
          }
       };

    template<>
       struct action<ProjectionLine> {
          static void apply(const pegtl::input& in, Solver<FMC_DT>& pd, std::stack<SIGNED_INDEX>& integerStack, std::stack<REAL>& realStack, MRFInput& mrfInput)
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
             pd.template GetProblemConstructor<1>().AddProjection(projectionVar,projectionCost);
          }
       };

    bool ParseProblem(const std::string& filename, Solver<FMC_DT>& pd) {
       std::stack<SIGNED_INDEX> integerStack;
       std::stack<REAL> realStack;
       MRFInput mrfInput;

       pegtl::file_parser problem(filename);
       std::cout << "parsing " << filename << "\n";

       return problem.parse< grammar, action>(pd, integerStack, realStack, mrfInput);
    }

  }

} // end namespace LP_MP
#endif // LP_MP_TOMOGRAPHY_H
