#include "graph_matching.h"

namespace LP_MP {

std::vector<INDEX> GraphMatching::ConstructMultiplexCapacityVector(const std::vector<INDEX>& graph, const std::pair<INDEX,INDEX> nodeCapacity) 
{
   const INDEX dim = nodeCapacity.first < nodeCapacity.second ? graph.size()+1 : graph.size();
   std::vector<INDEX> nodeCapacityVec(dim,1);
   if(dim > graph.size()) { nodeCapacityVec.back() = nodeCapacity.second - nodeCapacity.first; }
   return nodeCapacityVec;
}

GraphMatching::GraphMatching()
   : lp_(nullptr)
{}
GraphMatching::~GraphMatching()
{
   if(lp_ != nullptr)
      delete lp_;
}

// approach for node capacities: we use MultiplexFactor with the corresponding node capacities. For min_capacitiy < max_capacity, we add an additional node handling possibly non-assignment

void GraphMatching::AddNode(const INDEX nodeNumber, const INDEX min_capacity, const INDEX max_capacity, std::vector<std::vector<INDEX> >& graph, std::vector<std::pair<INDEX,INDEX> >& nodeCapacity)
{
   //assert(min_capacity == 0);
   assert(max_capacity == 1); // do zrobienia: increase this to two, but add ternary nodes then
   assert(min_capacity <= max_capacity);
   assert(graph.size() == nodeNumber); // we only allow nodes to be added in subsequent order
   graph.push_back(std::vector<INDEX>(0));
   nodeCapacity.push_back(std::make_pair(min_capacity,max_capacity));
}

void GraphMatching::AddLeftNode(const INDEX nodeNumber, const INDEX min_capacity, const INDEX max_capacity)
{
   AddNode(nodeNumber,min_capacity,max_capacity,leftGraph_,leftNodeCapacity_);
}

void GraphMatching::AddRightNode(const INDEX nodeNumber, const INDEX min_capacity, const INDEX max_capacity)
{
   AddNode(nodeNumber,min_capacity,max_capacity,rightGraph_,rightNodeCapacity_);
}

// note: this function has no symmetric formulation, as we need the unaries for the assignments either way.
void GraphMatching::AddAssignmentCost(const INDEX leftNode, const INDEX rightNode, const REAL cost)
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

}

void GraphMatching::AddToPairwisePotential(
         INDEX leftNode1, INDEX leftNode2, 
         INDEX rightNode1, INDEX rightNode2, 
         const REAL cost, 
         std::map<std::pair<INDEX,INDEX>, std::map<std::pair<INDEX,INDEX>, REAL> >& pairwisePotentials)
{
   // do zrobienia: generally true?
   assert(leftNode1 != leftNode2);
   assert(rightNode1 != rightNode2);
   if(leftNode1 > leftNode2) { std::swap(leftNode1,leftNode2); std::swap(rightNode1,rightNode2); }
   auto potIt = pairwisePotentials.find(std::make_pair(leftNode1,leftNode2));
   if(potIt == pairwisePotentials.end()) {
      std::map<std::pair<INDEX,INDEX>, REAL> leftPot;
      pairwisePotentials.insert( std::make_pair(std::make_pair(leftNode1,leftNode2), leftPot));
      potIt = pairwisePotentials.find(std::make_pair(leftNode1,leftNode2));
      assert(potIt != pairwisePotentials.end());
   }
   assert(potIt->second.find( std::make_pair(rightNode1, rightNode2) ) == potIt->second.end());
   potIt->second.insert(std::make_pair(std::make_pair(rightNode1,rightNode2), cost));
}

void GraphMatching::AddPairwiseCost(INDEX leftNode1, INDEX rightNode1, INDEX leftNode2, INDEX rightNode2, const REAL cost)
{
   assert(leftNode1 != leftNode2);
   assert(rightNode1 != rightNode2);
   assert( 0 <= leftNode1 && leftNode1 < leftGraph_.size()); 
   assert( 0 <= leftNode2 && leftNode2 < leftGraph_.size()); 
   assert( 0 <= rightNode1 && rightNode1 < rightGraph_.size()); 
   assert( 0 <= rightNode2 && rightNode2 < rightGraph_.size()); 

   // check if (leftNode1,rightNode1) and (leftNode2,rightNode2) are present
   assert( 1 == std::count(leftGraph_[leftNode1].begin(), leftGraph_[leftNode1].end(), rightNode1));
   assert( 1 == std::count(leftGraph_[leftNode2].begin(), leftGraph_[leftNode2].end(), rightNode2));

   AddToPairwisePotential(leftNode1,leftNode2,rightNode1,rightNode2, 0.5*cost, leftPairwisePotentials_);
   AddToPairwisePotential(rightNode1,rightNode2,leftNode1,leftNode2, 0.5*cost, rightPairwisePotentials_);
}
void GraphMatching::AddLeftPairwiseCost(INDEX leftNode1, INDEX rightNode1, INDEX leftNode2, INDEX rightNode2, const REAL cost)
{
   assert(leftNode1 != leftNode2);
   assert(rightNode1 != rightNode2);
   assert( 0 <= leftNode1 && leftNode2 < leftGraph_.size()); 
   assert( 0 <= rightNode1 && rightNode1 < rightGraph_.size()); 
   assert( 0 <= rightNode2 && rightNode2 < rightGraph_.size()); 
   
   AddToPairwisePotential(leftNode1,leftNode2,rightNode1,rightNode2, cost, leftPairwisePotentials_);
}
void GraphMatching::AddRightPairwiseCost(INDEX rightNode1, INDEX leftNode1, INDEX rightNode2, INDEX leftNode2, const REAL cost)
{
   assert(leftNode1 != leftNode2);
   assert(rightNode1 != rightNode2);
   assert( 0 <= leftNode1 && leftNode2 < leftGraph_.size()); 
   assert( 0 <= rightNode1 && rightNode1 < rightGraph_.size()); 
   assert( 0 <= rightNode2 && rightNode2 < rightGraph_.size()); 

   AddToPairwisePotential(rightNode1,rightNode2,leftNode1,leftNode2, cost, rightPairwisePotentials_);
}

void GraphMatching::AddToTernaryPotential(
         INDEX leftNode1, INDEX leftNode2, INDEX leftNode3,
         INDEX rightNode1, INDEX rightNode2, INDEX rightNode3, 
         const REAL cost, 
         std::map<std::tuple<INDEX,INDEX,INDEX>, std::map<std::tuple<INDEX,INDEX,INDEX>, REAL> >& ternaryPotentials)
{
   // order such that leftNode1 < leftNode2 < leftNode3
   if(leftNode1 > leftNode2) { std::swap(leftNode1,leftNode2); std::swap(rightNode1,rightNode2); }
   if(leftNode2 > leftNode3) { std::swap(leftNode2,leftNode3); std::swap(rightNode2,rightNode3); }
   if(leftNode1 > leftNode2) { std::swap(leftNode1,leftNode2); std::swap(rightNode1,rightNode2); }
   auto potIt = ternaryPotentials.find(std::make_tuple(leftNode1,leftNode2,leftNode3));
   if(potIt == ternaryPotentials.end()) {
      std::map<std::tuple<INDEX,INDEX,INDEX>, REAL> leftPot;
      ternaryPotentials.insert( std::make_pair(std::make_tuple(leftNode1,leftNode2,leftNode3), leftPot));
      potIt = ternaryPotentials.find(std::make_tuple(leftNode1,leftNode2,leftNode3));
      assert(potIt != ternaryPotentials.end());
   }
   potIt->second.insert(std::make_pair(std::make_tuple(rightNode1,rightNode2,rightNode2), cost));

}
void GraphMatching::AddTernaryCost(const INDEX leftNode1, const INDEX rightNode1, const INDEX leftNode2, const INDEX rightNode2, const INDEX leftNode3, const INDEX rightNode3, const REAL cost)
{
   assert(leftNode1 != leftNode2 && leftNode1 != leftNode3 && leftNode2 != leftNode3);
   assert( 0 <= leftNode1 && leftNode1 < leftGraph_.size()); 
   assert( 0 <= leftNode2 && leftNode2 < leftGraph_.size()); 
   assert( 0 <= leftNode3 && leftNode3 < leftGraph_.size()); 
   assert(rightNode1 != rightNode2 && rightNode1 != rightNode3 && rightNode2 != rightNode3);
   assert( 0 <= rightNode1 && rightNode1 < rightGraph_.size()); 
   assert( 0 <= rightNode2 && rightNode2 < rightGraph_.size()); 
   assert( 0 <= rightNode3 && rightNode3 < rightGraph_.size()); 

   assert( 1 == std::count(leftGraph_[leftNode1].begin(), leftGraph_[leftNode1].end(), rightNode1));
   assert( 1 == std::count(leftGraph_[leftNode2].begin(), leftGraph_[leftNode2].end(), rightNode2));
   assert( 1 == std::count(leftGraph_[leftNode3].begin(), leftGraph_[leftNode3].end(), rightNode3));

   AddLeftTernaryCost(leftNode1,leftNode2,leftNode3,rightNode1,rightNode2,leftNode3, 0.5*cost);
   AddRightTernaryCost(rightNode1,rightNode2,rightNode3,leftNode1,leftNode2,rightNode3, 0.5*cost);
}
void GraphMatching::AddLeftTernaryCost(const INDEX leftNode1, const INDEX rightNode1, const INDEX leftNode2, const INDEX rightNode2, const INDEX leftNode3, const INDEX rightNode3, const REAL cost)
{
   AddToTernaryPotential(leftNode1,leftNode2,leftNode3,rightNode1,rightNode2,rightNode3, 0.5*cost, leftTernaryPotentials_);
}
void GraphMatching::AddRightTernaryCost(const INDEX rightNode1, const INDEX leftNode1, const INDEX rightNode2, const INDEX leftNode2, const INDEX rightNode3, const INDEX leftNode3, const REAL cost)
{
   AddToTernaryPotential(rightNode1,rightNode2,rightNode3,leftNode1,leftNode2,leftNode3, 0.5*cost, rightTernaryPotentials_);
}

// do zrobienia: think about potts terms for multiplex unaries
void GraphMatching::AddLeftPottsTerm(const INDEX leftNode1, const INDEX leftNode2, const REAL cost)
{
   leftPottsCost_.insert(std::make_pair(std::make_pair(leftNode1,leftNode2), cost) );
}
void GraphMatching::AddRightPottsTerm(const INDEX rightNode1, const INDEX rightNode2, const REAL cost)
{
   rightPottsCost_.insert(std::make_pair(std::make_pair(rightNode1,rightNode2), cost) );
}

void GraphMatching::ConstructUnaryCost(
      std::vector<UnaryFactor*>& unaryMultiplex,
      std::vector<std::vector<REAL> >& unaryCost,
      const std::vector<std::pair<INDEX,INDEX> >& nodeCapacity,
      const std::vector<std::vector<INDEX> >& graph)
{
   unaryMultiplex.resize(graph.size());
   for(INDEX i=0; i<graph.size(); i++) {
      assert(unaryCost[i].size() == graph[i].size());
      std::vector<INDEX> factor_capacities(unaryCost[i].size(), 1);
      if( nodeCapacity[i].first < nodeCapacity[i].second ) {
         unaryCost[i].push_back(0);
         factor_capacities.push_back( nodeCapacity[i].second - nodeCapacity[i].first );
      }
      unaryMultiplex[i] = AddMultiplexFactor<UnaryFactor>( lp_, unaryCost[i], factor_capacities, nodeCapacity[i].second ); // do zrobienia: incorporate node capacities
      if(i>1) { // do zrobienia: niewiem czemu > 1 zamiast > 0
         //lp_->AddFactorRelation(unaryMultiplex[i-1], unaryMultiplex[i]); 
         //lp_->AddFactorRelation(unaryMultiplex[i], unaryMultiplex[i-1]); 
      }
   }
}

// construct pairwise cost from quadratic potentials as well as from Potts-terms
void GraphMatching::ConstructPairwiseCost(
      std::map<std::pair<INDEX,INDEX>, std::map<std::pair<INDEX,INDEX>, REAL> >& pairwisePotentials,
      std::map<std::pair<INDEX,INDEX>, REAL>& PottsCost,
      const std::vector<std::vector<INDEX> >& graph,
      const std::vector<std::pair<INDEX,INDEX> >& nodeCapacity,
      const std::vector<std::vector<INDEX> >& inverseGraph,
      //const std::vector<std::pair<INDEX,INDEX> >& inverseNodeCapacity,
      const std::vector<UnaryFactor*>& unaryMultiplex,
      std::map<std::pair<INDEX,INDEX>, PairwiseFactor*>& pairwiseMultiplex)
{
   auto PottsCostCopy = PottsCost; // we copy PottsCost to remove elements already taken care of. The elements not taken care of are put in separate factors
   
   // first construct pairwise cost from pairwisePotentials_
   for(auto pIt = pairwisePotentials.begin(); pIt!=pairwisePotentials.end(); ++pIt) {
      const INDEX i1 = pIt->first.first;
      const INDEX i2 = pIt->first.second;
      assert(i1<i2);
      // account for unmatched node
      const INDEX leftDim = nodeCapacity[i1].first < nodeCapacity[i1].second ? graph[i1].size()+1 : graph[i1].size();
      const INDEX rightDim = nodeCapacity[i2].first < nodeCapacity[i2].second ? graph[i2].size()+1 : graph[i2].size();
      std::vector<REAL> cost(leftDim * rightDim, 0.0);
      std::vector<bool> cost_touched(graph[i1].size()*graph[i2].size(), false); // do zrobienia: remove this!
      REAL min_cost = std::numeric_limits<REAL>::max();
      for(auto it = pIt->second.begin(); it!= pIt->second.end(); ++it) {
         const INDEX j1_index = find_index(it->first.first, graph[i1]);
         const INDEX j2_index = find_index(it->first.second, graph[i2]);
         assert(j1_index + j2_index*leftDim < cost.size());
         cost[j1_index + j2_index*leftDim] += it->second;
         //std::cout << "pairwise cost = " << it->second << std::endl;
         min_cost = std::min(it->second, min_cost);
         cost_touched[j1_index + j2_index*graph[i1].size()] = true;
      }
      //std::cout << "min cost = " << min_cost << std::endl;
      // set diagonal to infinity
      // do zrobienia: do not necessarily set diagonal to infinity, if allowing nodeCapacity > 1
      // but then the ternary cost has to be set to infinity for the diagonal
      for(INDEX i1_index=0; i1_index<graph[i1].size(); i1_index++) {
         for(INDEX i2_index=0; i2_index<graph[i2].size(); i2_index++) {
            if(graph[i1][i1_index] == graph[i2][i2_index]) {
               assert(i1_index + i2_index*leftDim < cost.size());
               cost[i1_index + i2_index*leftDim] = std::numeric_limits<REAL>::max();
            } else {
               // here we check whether all entries for pairwise costs had been set. This is not necessary!
               //assert( cost_touched[i1_index + i2_index*graph[i1].size()] == true );
            }
         }
      }
      std::vector<INDEX> leftNodeCapacity(ConstructMultiplexCapacityVector(graph[i1], nodeCapacity[i1]));
      std::vector<INDEX> rightNodeCapacity(ConstructMultiplexCapacityVector(graph[i2], nodeCapacity[i2]));

      // see if there is a potts term with the same indices
      if(PottsCost.count(std::make_pair(i1,i2)) > 0) {
         // directly encode Potts cost INDEXo quadratic
         assert( nodeCapacity[i1].first == 0 || nodeCapacity[i2].first == 0); // otherwise the term would be constant 
         const REAL p = PottsCost[std::make_pair(i1,i2)];

         if( nodeCapacity[i1].first == 0 ) {
            for(INDEX i2_index = 0; i2_index < graph[i2].size(); i2_index++) {
               assert(graph[i1].size() + i2_index*leftDim < cost.size());
               cost[graph[i1].size() + i2_index*leftDim] += p;
            }
         }
         if( nodeCapacity[i2].first == 0 ) {
            for(INDEX i1_index = 0; i1_index < graph[i1].size(); i1_index++) {
               assert(i1_index + graph[i2].size()*leftDim < cost.size()); 
               cost[i1_index + graph[i2].size()*leftDim] += p;
            }
         }
         PottsCostCopy.erase( PottsCostCopy.find(std::make_pair(i1,i2)) );
      }

      PairwiseFactor* f = AddPairwiseMultiplexFactor<PairwiseFactor>(lp_, cost, leftNodeCapacity, rightNodeCapacity, nodeCapacity[i1].second, nodeCapacity[i2].second);
      assert( pairwiseMultiplex.count( std::make_pair(i1,i2) ) == 0);
      pairwiseMultiplex.insert(std::make_pair( std::make_pair(i1,i2), f));

      
      LinkUnaryPairwiseMultiplexFactors<UnaryPairwiseMessageLeft>(lp_, unaryMultiplex[i1], f, leftDim, rightDim, 0);
      LinkUnaryPairwiseMultiplexFactors<UnaryPairwiseMessageRight>(lp_, unaryMultiplex[i2], f, leftDim, rightDim, 1);

      //lp_->AddFactorRelation(f,unaryMultiplex[i1]); 
      //lp_->AddFactorRelation(unaryMultiplex[i2],f);
      //lp_->AddFactorRelation(unaryMultiplex[i1],f); 
      //lp_->AddFactorRelation(f,unaryMultiplex[i2]);
   }
   // now take care of the remaining potts entries to which there was no other quadratic cost
   // do zrobienia: special factor for this
   for( auto PottsIt : PottsCostCopy ) {
      throw std::logic_error("not implemented yet");
   }

   return;
}

void GraphMatching::ConstructTernaryCost(
      std::map<std::tuple<INDEX,INDEX,INDEX>, std::map<std::tuple<INDEX,INDEX,INDEX>, REAL> > ternaryPotentials,
      const std::vector<std::vector<INDEX> >& graph,
      const std::vector<std::pair<INDEX,INDEX> >& nodeCapacity,
      const std::vector<std::vector<INDEX> >& inverseGraph,
      const std::map<std::pair<INDEX,INDEX>, PairwiseFactor*>& pairwiseMultiplexFactor)
{
   for(auto pIt = ternaryPotentials.begin(); pIt!=ternaryPotentials.end(); ++pIt) {
      assert(false); // not implemented yet
      const INDEX i1 = std::get<0>(pIt->first);
      const INDEX i2 = std::get<1>(pIt->first);
      const INDEX i3 = std::get<2>(pIt->first);
      assert(i1<i2 && i2<i3);
      // account for unmatched nodes
      const INDEX dim1 = nodeCapacity[i1].first < nodeCapacity[i1].second ? graph[i1].size()+1 : graph[i1].size();
      const INDEX dim2 = nodeCapacity[i2].first < nodeCapacity[i2].second ? graph[i2].size()+1 : graph[i2].size();
      const INDEX dim3 = nodeCapacity[i3].first < nodeCapacity[i3].second ? graph[i3].size()+1 : graph[i3].size();
      std::vector<REAL> cost(dim1 * dim2 * dim3, 0.0);
      for(auto it = pIt->second.begin(); it!= pIt->second.end(); ++it) {
         const INDEX j1_index = find_index(std::get<0>(it->first), graph[i1]);
         const INDEX j2_index = find_index(std::get<1>(it->first), graph[i2]);
         const INDEX j3_index = find_index(std::get<2>(it->first), graph[i3]);

         cost[j1_index + j2_index*dim1 + j3_index*dim1*dim2] += it->second; 

         // add infinity to diagonal
         // do zrobienia: take INDEXo account nodeCapacity = 2
         for(INDEX i=0; i<graph[i1].size(); i++) {
            for(INDEX j=0; j<graph[i2].size(); j++) {
               for(INDEX k=0; k<graph[i3].size(); k++) {
                  if(graph[i1][i] == graph[i2][j] || graph[i1][i] == graph[i3][k] || graph[i2][j] == graph[i3][k]) {
                     cost[i + j*dim1 + k*dim1*dim2] = std::numeric_limits<REAL>::max();
                  }
               }
            }
         }
      }

      std::vector<INDEX> nodeCapacity1 = ConstructMultiplexCapacityVector(graph[i1], nodeCapacity[i1]);
      std::vector<INDEX> nodeCapacity2 = ConstructMultiplexCapacityVector(graph[i2], nodeCapacity[i2]);
      std::vector<INDEX> nodeCapacity3 = ConstructMultiplexCapacityVector(graph[i3], nodeCapacity[i3]);
      //---------------------------------------------------------------------------------------------------------------------------------------------------
      //MultiplexFactor* f = AddTernaryMultiplexFactor(lp_, cost, nodeCapacity1, nodeCapacity2, nodeCapacity3, nodeCapacity[i1].second, nodeCapacity[i2].second, nodeCapacity[i3].second);
      //---------------------------------------------------------------------------------------------------------------------------------------------------

      assert( pairwiseMultiplexFactor.find( std::make_pair(i1,i2) ) != pairwiseMultiplexFactor.end());
      assert( pairwiseMultiplexFactor.find( std::make_pair(i1,i3) ) != pairwiseMultiplexFactor.end());
      assert( pairwiseMultiplexFactor.find( std::make_pair(i2,i3) ) != pairwiseMultiplexFactor.end());

      PairwiseFactor* factor_12 = pairwiseMultiplexFactor.find( std::make_pair(i1,i2) )->second;
      PairwiseFactor* factor_13 = pairwiseMultiplexFactor.find( std::make_pair(i1,i3) )->second;
      PairwiseFactor* factor_23 = pairwiseMultiplexFactor.find( std::make_pair(i2,i3) )->second;
      //---------------------------------------------------------------------------------------------------------------------------------------------------
      //LinkTernaryMultiplexFactors(lp_, f, factor_12, factor_13, factor_23, dim1 , dim2, dim3);
      //---------------------------------------------------------------------------------------------------------------------------------------------------

      // do zrobienia: add ordering
   }
}

void GraphMatching::ConstructAssignmentConstraints()
{
   for(INDEX i_0=0; i_0<leftGraph_.size(); i_0++) {
      for(INDEX i_0_index=0; i_0_index<leftGraph_[i_0].size(); i_0_index++) {
         const INDEX i_1 = leftGraph_[i_0][i_0_index];
         const INDEX i_1_index = find(rightGraph_[i_1].begin(), rightGraph_[i_1].end(), i_0) - rightGraph_[i_1].begin();
         assert(i_1_index < rightGraph_[i_1].size());
         assert(rightGraph_[i_1][i_1_index] == i_0);
         // do zrobienia: account for node capacities > 1
         LinkSingleVariableEqualityMessage<AssignmentConstraintMessage>(lp_, leftGraphUnaryMultiplex_[i_0], i_0_index, rightGraphUnaryMultiplex_[i_1], i_1_index);
      }
   }
}


std::vector<int> GraphMatching::GetAssignment(
      const std::vector<std::vector<REAL> >& cost, 
      const std::vector<std::vector<INDEX> >& graph,
      const std::vector<std::pair<INDEX,INDEX> >& nodeCapacity,
      const std::vector<std::pair<INDEX,INDEX> >& inverseNodeCapacity
      )
{
   assert(nodeCapacity.size() == graph.size());

   INDEX totalLeftCapacity = 0;
   for(INDEX i=0; i<nodeCapacity.size(); i++) {
      totalLeftCapacity += nodeCapacity[i].second;
   }
   INDEX totalRightCapacity = 0;
   for(INDEX j=0; j<inverseNodeCapacity.size(); j++) {
      totalRightCapacity += inverseNodeCapacity[j].second;
   }

   INDEX nodeNum = graph.size() + 2; // + left and right surplus node
   INDEX rightNodeNum = 0;
   INDEX edgeNum = 0;
   assert(graph.size() == cost.size());
   for(INDEX i=0; i<graph.size(); ++i) {
      for(INDEX j=0; j<graph[i].size(); ++j) {
         rightNodeNum = std::max(INDEX(graph[i][j]+1), rightNodeNum);
      }
      edgeNum += graph[i].size();
   }
   nodeNum += rightNodeNum;
   edgeNum += rightNodeNum + 2*graph.size() + 3;
   MinCost<int,REAL> m( nodeNum, edgeNum ); // note: must be signed, otherwise does not terminate
   // regular edges
   for(INDEX i=0; i<graph.size(); i++) {
      for(INDEX j=0; j<graph[i].size(); j++) {
         m.AddEdge(i, graph.size() + graph[i][j], 1.0, 0.0, cost[i][j]);
      }
   }
   //edges to surplus node
   for(INDEX i=0; i<graph.size(); i++) {
      if(nodeCapacity[i].second > nodeCapacity[i].first) {
         assert(cost[i].size() == graph[i].size() + 1);
         m.AddEdge(i, nodeNum-2, nodeCapacity[i].second - nodeCapacity[i].first, 0.0, cost[i].back());
      }
   }
   for(INDEX j=0; j<rightNodeNum; j++) {
      if(inverseNodeCapacity[j].second > inverseNodeCapacity[j].first) {
         m.AddEdge(nodeNum-1, graph.size() + j, inverseNodeCapacity[j].second - inverseNodeCapacity[j].first, 0, 0.0); // do zrobienia: is cost 0.0 correct?
      }
   }
   m.AddEdge(nodeNum-1, nodeNum-2, std::max(totalLeftCapacity, totalRightCapacity),  std::max(totalLeftCapacity, totalRightCapacity), 0.0);
   m.AddNodeExcess(nodeNum-1, totalRightCapacity);
   m.AddNodeExcess(nodeNum-2, -totalLeftCapacity);

   // ensure assignment
   for(INDEX i=0; i<graph.size(); i++) {
      m.AddNodeExcess(i, nodeCapacity[i].second);
      //m.AddEdge(nodeNum-2, i, nodeCapacity[i].second, nodeCapacity[i].second, 0.0);
   }
   for(INDEX j=0; j<rightNodeNum; j++) {
      m.AddNodeExcess(graph.size() + j, -inverseNodeCapacity[j].second);
      //m.AddEdge(graph.size() + j, nodeNum-1, inverseNodeCapacity[j].second, 0, 0.0);
   }

   m.Solve();
   std::vector<int> assignment(graph.size(), -1); // we must take int, as negative values are allowed as well
   INDEX c=0;
   for(INDEX i=0; i<graph.size(); i++) {
      for(INDEX j=0; j<graph[i].size(); j++) {
         if(std::abs(m.GetRCap(c)) < 0.00001) {
            assignment[i] = graph[i][j];
         }
         c++;
      }
   }
   //assert(HasUniqueValues(assignment));
   for(INDEX i=0; i<assignment.size(); i++) { 
      assert(assignment[i] >= -1); 
   }
   return assignment;
}

std::vector<int> GraphMatching::GetLeftAssignment()
{
   assert(lp_ != NULL);
   std::vector<std::vector<REAL> > repam = lp_->GetReparametrizedModel();
   repam.resize(leftGraph_.size());  // do zrobienia: assumes first left unaries were added, then right ones
   std::vector<int> assignment = GetAssignment(repam, leftGraph_, leftNodeCapacity_, rightNodeCapacity_);
   //for(INDEX i=0; i<assignment.size(); i++) {
   //   std::cout << i << " -> " << assignment[i] << std::endl;
   //}
   return assignment;
}

std::vector<int> GraphMatching::GetRightAssignment()
{
   assert(lp_ != NULL);
   std::vector<std::vector<REAL> > repamFull = lp_->GetReparametrizedModel();
   std::vector<std::vector<REAL> > repam(repamFull.begin() + leftGraph_.size(), repamFull.begin() + leftGraph_.size() + rightGraph_.size());  // do zrobienia: assumes first left unaries were added, then right ones
   std::vector<int> assignment = GetAssignment(repam, rightGraph_, rightNodeCapacity_, leftNodeCapacity_);
   //for(INDEX i=0; i<assignment.size(); i++) {
   //   std::cout << i << " -> " << assignment[i] << std::endl;
   //}
   return assignment;
}

REAL GraphMatching::evalLeft(const std::vector<int>& assignment)
{
   return eval(assignment, leftGraph_, unaryCostLeft_, leftPairwisePotentials_);
}

REAL GraphMatching::evalRight(const std::vector<int>& assignment)
{
   return eval(assignment, rightGraph_, unaryCostRight_, rightPairwisePotentials_);
}

REAL GraphMatching::eval(
      const std::vector<int>& assignment,
      const std::vector<std::vector<INDEX> >& graph, 
      const std::vector<std::vector<REAL> >& unaryCost,
      const std::map<std::pair<INDEX,INDEX>, std::map<std::pair<INDEX,INDEX>, REAL> >& pairwisePotentials)
   // ternary potentials
{
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
}

void GraphMatching::Solve(const INDEX nIter)
{
   if(lp_ != nullptr) delete lp_;
   lp_ = new LP();


   ConstructUnaryCost(leftGraphUnaryMultiplex_, unaryCostLeft_, leftNodeCapacity_, leftGraph_);
   ConstructUnaryCost(rightGraphUnaryMultiplex_, unaryCostRight_, rightNodeCapacity_, rightGraph_);
   std::cout << leftGraphUnaryMultiplex_.size() + rightGraphUnaryMultiplex_.size() << " unary potentials = " << leftGraphUnaryMultiplex_.size() << " left unary potentials + " << rightGraphUnaryMultiplex_.size() << " right unary potentials." << std::endl;

   ConstructAssignmentConstraints();

   ConstructPairwiseCost(leftPairwisePotentials_,  leftPottsCost_,  leftGraph_,  leftNodeCapacity_,  rightGraph_, leftGraphUnaryMultiplex_,  leftGraphPairwiseMultiplex_);
   ConstructPairwiseCost(rightPairwisePotentials_, rightPottsCost_, rightGraph_, rightNodeCapacity_, leftGraph_,  rightGraphUnaryMultiplex_, rightGraphPairwiseMultiplex_);
   std::cout << leftGraphPairwiseMultiplex_.size() + rightGraphPairwiseMultiplex_.size() << " pairwise potentials = " << leftGraphPairwiseMultiplex_.size() << " left pairwise potentials + " << rightGraphPairwiseMultiplex_.size() << " right pairwise potentials." << std::endl;

   //ConstructTernaryCost(leftTernaryPotentials_,leftGraph_,leftNodeCapacity_,rightGraph_,leftGraphPairwiseMultiplex_);
   //ConstructTernaryCost(rightTernaryPotentials_,rightGraph_,rightNodeCapacity_,leftGraph_,rightGraphPairwiseMultiplex_);
   //std::cout << leftGraphTernaryMultiplex_.size() + rightGraphTernaryMultiplex_.size() << " ternary potentials = " << leftGraphTernaryMultiplex_.size() << " left ternary potentials + " << rightGraphTernaryMultiplex_.size() << " right ternary potentials." << std::endl;

   lp_->Solve(nIter);
   //delete lp_;
}


} // end namespace LP_MP
