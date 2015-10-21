#ifndef ASSIGNMENT_PROBLEM_CONSTRUCTION_HXX
#define ASSIGNMENT_PROBLEM_CONSTRUCTION_HXX

#include "LP_MP.h"
#include "factors_messages.hxx"
#include "MinCost/MinCost.h"
#include "problem_construction_helper.hxx"
#include "problem_decomposition.hxx"

#include <sstream>

namespace LP_MP {

// construct assignment problem in some prespecified factor message network given the number of the assignment factors and the equality message joining them
// do zrobienia: move into include/...
// do zrobienia: integrate with linear_assignment to obtain primal solution from reparametrization
// do zrobienia: add eval and primalBound
// do zrobienia: given a line from problem file, read in contents and call the right methods

template<class FMC, INDEX FACTOR_NO, INDEX MESSAGE_NO>
class AssignmentProblemConstructor {
   using AssignmentFactor = meta::at_c<typename FMC::factor_list, FACTOR_NO>;
   using AssignmentConstraintMessage = meta::at_c<meta::at_c<typename FMC::msg_list, MESSAGE_NO>, 0>;
public:
   AssignmentProblemConstructor() {}

   void ReadLine(const std::string& line) {
      char assignment_char;
      INDEX l, r;
      REAL cost;
      std::stringstream s{line};
      s >> assignment_char >> l >> r >> cost;
      if(s.fail()) throw std::runtime_error("incompatible line");
      if(assignment_char != 'a') throw std::runtime_error("assignment must start with letter 'a'");
      if(l >= leftGraph_.size()) { SetNumberLeftNodes(l); }
      if(r >= rightGraph_.size()) { SetNumberRightNodes(l); }
      AddAssignmentCost(l,r,cost);
   }

   AssignmentProblemConstructor& SetNumberNodes(const INDEX nodeNumber, std::vector<std::vector<INDEX> >& graph, std::vector<INDEX>& marginal)
   {
      assert(nodeNumber >= graph.size());
      graph.resize(nodeNumber+1, std::vector<INDEX>(0));
      marginal.resize(nodeNumber+1, 1);
      return *this;
   }
   AssignmentProblemConstructor& SetNumberLeftNodes(const INDEX nodeNumber) {
      SetNumberNodes(nodeNumber,leftGraph_, leftMarginal_);
      return *this;
   }
   AssignmentProblemConstructor& SetNumberRightNodes(const INDEX nodeNumber) {
      SetNumberNodes(nodeNumber,rightGraph_, rightMarginal_);
      return *this;
   }

   AssignmentProblemConstructor& SetLeftMarginal(const INDEX i, const INDEX val) {
      leftMarginal_[i] = val;
      return *this;
   }
   AssignmentProblemConstructor& SetRightMarginal(const INDEX i, const INDEX val) {
      rightMarginal_[i] = val;
      return *this;
   }

   // note: this function has no symmetric formulation, as we need the unaries for the assignments either way.
   AssignmentProblemConstructor& AddAssignmentCost(const INDEX leftNode, const INDEX rightNode, const REAL cost)
   {
      leftUnaryCost_.resize(leftGraph_.size());
      assert(leftNode < leftGraph_.size());
      assert(rightNode < rightGraph_.size());
      leftGraph_[leftNode].push_back(rightNode);
      assert(HasUniqueValues(leftGraph_[leftNode]));
      leftUnaryCost_[leftNode].push_back(0.5*cost);

      rightUnaryCost_.resize(rightGraph_.size());
      rightGraph_[rightNode].push_back(leftNode);
      assert(HasUniqueValues(rightGraph_[rightNode]));
      rightUnaryCost_[rightNode].push_back(0.5*cost);

      for(INDEX i=0; i<leftGraph_.size(); i++) {
         assert(HasUniqueValues(leftGraph_[i]));
      }
      for(INDEX i=0; i<rightGraph_.size(); i++) {
         assert(HasUniqueValues(rightGraph_[i]));
      }

      return *this;
   }

   void Construct(ProblemDecomposition<FMC>& pd) {}

   AssignmentProblemConstructor& Construct(LP* lp) {
      ConstructAssignmentFactors(lp, leftUnaryCost_, leftMarginal_, leftAssignmentFactor_);
      ConstructAssignmentFactors(lp, rightUnaryCost_, rightMarginal_, rightAssignmentFactor_);
      ConstructAssignmentConstraints(lp, leftGraph_, leftAssignmentFactor_, rightGraph_, rightAssignmentFactor_);
      return *this;
   }

   const std::vector<AssignmentFactor*>& GetLeftAssignmentFactors() { return leftAssignmentFactor_; }
   const std::vector<AssignmentFactor*>& GetRightAssignmentFactors() { return rightAssignmentFactor_; }

private:

   void ConstructAssignmentFactors(
         LP* const lp,
         const std::vector<std::vector<REAL> >& cost,
         const std::vector<INDEX>& marginal,
         std::vector<AssignmentFactor*>& assignmentFactor)
   {
      assignmentFactor.resize(cost.size());
      for(INDEX i=0; i<cost.size(); ++i) {
         std::vector<REAL> cur_cost = cost[i];
         cur_cost.push_back(0.0); // cost for non-assignment
         std::vector<INDEX> capacity(cur_cost.size(), marginal[i]);
         assignmentFactor[i] = AddMultiplexFactor<AssignmentFactor>(lp, cur_cost, capacity, marginal[i]);

         /*
         std::cout << *std::min_element(cost[i].begin(), cost[i].end()) << "\n";
         for(INDEX j=0; j<cur_cost.size(); ++j) {
            std::cout << cur_cost[j] << ", ";
         }
         std::cout << "\n\n";
         */

      }
   }

   void ConstructAssignmentConstraints(
         LP* const lp,
         const std::vector<std::vector<INDEX> >& graph1,
         const std::vector<AssignmentFactor*>& assignmentFactor1,
         const std::vector<std::vector<INDEX> >& graph2,
         const std::vector<AssignmentFactor*>& assignmentFactor2)
   {
      for(INDEX i_0=0; i_0<graph1.size(); i_0++) {
         for(INDEX i_0_index=0; i_0_index<graph1[i_0].size(); i_0_index++) {
            const INDEX i_1 = graph1[i_0][i_0_index];
            const INDEX i_1_index = find(graph2[i_1].begin(), graph2[i_1].end(), i_0) - graph2[i_1].begin();
            assert(i_1_index < graph2[i_1].size());
            assert(graph2[i_1][i_1_index] == i_0);
            LinkSingleVariableEqualityMessage<AssignmentConstraintMessage>(lp, assignmentFactor1[i_0], i_0_index, assignmentFactor2[i_1], i_1_index);
         }
      }
   }

private:
   std::vector<std::vector<INDEX> > leftGraph_, rightGraph_;

   std::vector<AssignmentFactor*> leftAssignmentFactor_;
   std::vector<AssignmentFactor*> rightAssignmentFactor_;

   // value to which {left|right}AssignmentFactor_ must sum
   std::vector<INDEX> leftMarginal_;
   std::vector<INDEX> rightMarginal_;

   std::vector<std::vector<REAL> > leftUnaryCost_;
   std::vector<std::vector<REAL> > rightUnaryCost_;

};




} // end namespace LP_MP

#endif // ASSIGNMENT_PROBLEM_CONSTRUCTION_HXX
