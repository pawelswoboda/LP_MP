#include "cosegmentation.h"

namespace LP_MP {


Cosegmentation& Cosegmentation::SetNumberLeftNodes(const INDEX nodeNumber)
{
   SetNumberNodes(nodeNumber,leftGraph_);
   return *this;
}
Cosegmentation& Cosegmentation::SetNumberRightNodes(const INDEX nodeNumber)
{
   SetNumberNodes(nodeNumber,rightGraph_);
   return *this;
}
Cosegmentation& Cosegmentation::SetNumberNodes(const INDEX nodeNumber, std::vector<std::vector<INDEX> >& graph)
{
   assert(nodeNumber >= graph.size());
   graph.resize(nodeNumber, std::vector<INDEX>(0));
   return *this;
}

// note: this function has no symmetric formulation, as we need the unaries for the assignments either way.
Cosegmentation& Cosegmentation::AddAssignmentCost(const INDEX leftNode, const INDEX rightNode, const REAL cost)
{
   unaryCostLeft_.resize(leftGraph_.size());
   assert(leftNode < leftGraph_.size());
   assert(rightNode < rightGraph_.size());
   leftGraph_[leftNode].push_back(rightNode);
   assert(HasUniqueValues(leftGraph_[leftNode]));
   unaryCostLeft_[leftNode].push_back(0.5*cost);

   unaryCostRight_.resize(rightGraph_.size());
   rightGraph_[rightNode].push_back(leftNode);
   assert(HasUniqueValues(rightGraph_[rightNode]));
   unaryCostRight_[rightNode].push_back(0.5*cost);

   for(INDEX i=0; i<leftGraph_.size(); i++) {
      assert(HasUniqueValues(leftGraph_[i]));
   }
   for(INDEX i=0; i<rightGraph_.size(); i++) {
      assert(HasUniqueValues(rightGraph_[i]));
   }

   return *this;
}

Cosegmentation& Cosegmentation::AddLeftPottsTerm(const INDEX leftNode1, const INDEX leftNode2, const REAL cost)
{
   leftPottsCost_.insert(std::make_pair(std::make_pair(leftNode1,leftNode2), cost) );
   return *this;
}
Cosegmentation& Cosegmentation::AddRightPottsTerm(const INDEX rightNode1, const INDEX rightNode2, const REAL cost)
{
   rightPottsCost_.insert(std::make_pair(std::make_pair(rightNode1,rightNode2), cost) );
   return *this;
}

std::vector<int> Cosegmentation::GetLeftSegmentation()
{
   assert(lp_ != NULL);
   std::vector<std::vector<REAL> > repam = lp_->GetReparametrizedModel();
   repam.resize(leftUnaryFactor_.size());  // do zrobienia: assumes first left unaries were added, then right ones
   std::vector<int> segmentation(repam.size(),0);
   for(INDEX i=0; i<repam.size(); ++i) {
      assert(repam[i].size() == 2);
      if(repam[i][0] < repam[i][1]) {
         segmentation[i] = 1;
      }
   }
   
   return segmentation;
}

std::vector<int> Cosegmentation::GetRightSegmentation()
{
   assert(lp_ != NULL);
   std::vector<std::vector<REAL> > repamFull = lp_->GetReparametrizedModel();
   const INDEX leftSize = leftUnaryFactor_.size() + leftPairwiseFactor_.size(); 
   std::vector<std::vector<REAL> > repam(repamFull.begin() + leftSize, repamFull.begin() + leftSize + rightUnaryFactor_.size());  // do zrobienia: assumes first left unaries were added, then right ones
   std::vector<int> segmentation(repam.size(),0);
   for(INDEX i=0; i<repam.size(); ++i) {
      assert(repam[i].size() == 2);
      if(repam[i][0] < repam[i][1]) {
         segmentation[i] = 1;
      }
   }
   return segmentation;
}

REAL Cosegmentation::evalLeft(const std::vector<int>& segmentation)
{
   return 0.0;
   //return eval(assignment, leftGraph_, unaryCostLeft_, leftPottsCost_);
}

REAL Cosegmentation::evalRight(const std::vector<int>& segmentation)
{
   return 0.0;
   //return eval(assignment, rightGraph_, unaryCostRight_, rightPairwisePotentials_);
}

REAL Cosegmentation::evalAssignment(const std::vector<int>& assignment)
{
   return 0.0;  
}

REAL Cosegmentation::eval(
      const std::vector<int>& assignment,
      const std::vector<std::vector<INDEX> >& graph, 
      const std::vector<std::vector<REAL> >& unaryCost,
      const std::map<std::pair<INDEX,INDEX>, REAL>& pottsCost)
{
   return 0.0;
   /*
   assert(assignment.size() == graph.size());
   std::vector<int> assignment_index(assignment);
   for(INDEX i=0; i<assignment_index.size(); ++i) {
      if(assignment_index[i] != -1) {
         assignment_index[i] = find_index(INDEX(assignment[i]), graph[i]);
         assert(assignment_index[i] < graph[i].size());
      }
   }
   REAL cost = 0.0;
   INDEX maxLabel = 0;
   for(INDEX i=0; i<assignment.size(); ++i) {
      if(assignment_index[i] != -1) {
         cost += unaryCost[i][assignment_index[i]];
         maxLabel = std::max(maxLabel, graph[i][assignment_index[i]]);
      }
   }
   // check if no assignment constraint is violated
   std::vector<bool> labelsTaken(maxLabel+1,false);
   for(INDEX i=0; i<assignment.size(); ++i) {
      if( assignment[i] == -1) continue;
      if( labelsTaken[ graph[i][assignment_index[i]] ] == true ) return std::numeric_limits<REAL>::max(); // do zrobienia: relax this and account for capacities > 1
      labelsTaken[ graph[i][assignment_index[i]] ] = true;
   }
   for(auto pot=pairwisePotentials.begin(); pot!=pairwisePotentials.end(); ++pot) {
      const INDEX i1 = pot->first.first;
      const INDEX i2 = pot->first.second;
      assert(i1<i2);
      const INDEX j1 = assignment[i1];
      const INDEX j2 = assignment[i2];
      auto costIt = pot->second.find(std::make_pair(j1,j2));
      if(costIt != pot->second.end()) {
         cost += costIt->second;
      }
   }
   // do zrobienia: add triplet potentials
   // do zrobienia: cost need not be true anymore, as it is not shared on left and right uniformly. consider left and right simultaneously
   return 2.0*cost;
   */
}

void Cosegmentation::ConstructGraph(
      const std::map<std::pair<INDEX,INDEX>, REAL>& PottsCost,
      std::vector<UnaryFactor*>& unaryFactor, 
      std::vector<PairwiseFactor*>& pairwiseFactor)
{
   unaryFactor.clear();
   pairwiseFactor.clear();
   // first construct all unaries
   INDEX numberUnaries = 0;
   for(auto PottsIt : PottsCost) {
      const INDEX i = PottsIt.first.first;
      const INDEX j = PottsIt.first.second;
      numberUnaries = std::max(numberUnaries, i);
      numberUnaries = std::max(numberUnaries, j);
   }
   ++numberUnaries; // previously, only indices were counted. They start at 0, hence add 1
   unaryFactor.resize(numberUnaries);
   for(INDEX i=0; i<numberUnaries; ++i) {
      unaryFactor[i] = AddSimplexFactor<UnaryFactor>(lp_, std::vector<REAL>(2,0.0));
   }
   // add pairwise Potts factors
   pairwiseFactor.reserve(PottsCost.size());
   for(auto PottsIt : PottsCost) {
      const REAL c = PottsIt.second;
      const INDEX i = PottsIt.first.first;
      const INDEX j = PottsIt.first.second;

      pairwiseFactor.push_back( AddSimplexFactor<PairwiseFactor>(lp_, std::vector<REAL>({0.0,c,c,0.0})) );
      LinkUnaryPairwiseMultiplexFactors<UnaryPairwiseMessageLeft>(lp_, unaryFactor[i], pairwiseFactor.back(), 2, 2, 0);
      LinkUnaryPairwiseMultiplexFactors<UnaryPairwiseMessageRight>(lp_, unaryFactor[j], pairwiseFactor.back(), 2, 2, 1);
   }
}

void Cosegmentation::ConstructAssignmentFactors(
      const std::vector<std::vector<REAL> >& cost,
      std::vector<AssignmentFactor*>& assignmentFactor)
{
   assignmentFactor.resize(cost.size());
   for(INDEX i=0; i<cost.size(); ++i) {
      std::vector<REAL> cur_cost = cost[i];
      cur_cost.push_back(0.0); // cost for non-assignment
      assignmentFactor[i] = AddSimplexFactor<AssignmentFactor>(lp_, cur_cost);
   }
}

void Cosegmentation::ConstructSummationConstraints(
   const std::vector<UnaryFactor*>& unaryFactor,
   const std::vector<AssignmentFactor*>& assignmentFactor)
{
   assert(unaryFactor.size() == assignmentFactor.size());
   for(INDEX i=0; i<unaryFactor.size(); ++i) {
      auto m = new MarginalSummationMessageContainer(MarginalSummationMessage(), unaryFactor[i], assignmentFactor[i], 2);
      lp_->AddMessage(m);
   }
}

void Cosegmentation::ConstructAssignmentConstraints(
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
         LinkSingleVariableEqualityMessage<AssignmentConstraintMessage>(lp_, assignmentFactor1[i_0], i_0_index, assignmentFactor2[i_1], i_1_index);
      }
   }
}


Cosegmentation& Cosegmentation::Solve(const INDEX nIter)
{
   if(lp_ != nullptr) delete lp_;
   lp_ = new LP();

   // Construct the MAP-MRF without the assignment
   ConstructGraph(leftPottsCost_,leftUnaryFactor_,leftPairwiseFactor_);
   ConstructGraph(rightPottsCost_, rightUnaryFactor_, rightPairwiseFactor_);

   ConstructAssignmentFactors(unaryCostLeft_, leftAssignmentFactor_);
   ConstructAssignmentFactors(unaryCostRight_, rightAssignmentFactor_);

   ConstructSummationConstraints(leftUnaryFactor_, leftAssignmentFactor_);
   ConstructSummationConstraints(rightUnaryFactor_, rightAssignmentFactor_);

   ConstructAssignmentConstraints(leftGraph_,leftAssignmentFactor_, rightGraph_, rightAssignmentFactor_);

   lp_->Solve(nIter);
   //delete lp_;
   return *this;
}


} // end namespace LP_MP

