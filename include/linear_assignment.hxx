#ifndef LP_MP_LINEAR_ASSIGNMENT_HXX
#define LP_MP_LINEAR_ASSIGNMENT_HXX

#include "instances.inc"
#include "LP_MP.h"
#include "MinCost/MinCost.h"

namespace LP_MP {

// do zrobienia: remove this file, functionality is in assignment_problem_constructor.hxx
// AddEdge has been changed in MinCost.h

// template argument only needed because of hxx file. do zrobienia: templatize for cost and graph vector to allow other representations as well
template<bool unneeded = true>
std::vector<int> GetAssignment(
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
   MinCost<int,REAL> m( nodeNum, edgeNum ); // note: first template argument must be signed, otherwise algorithm does not terminate

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

// for linear assignment problem with capacities one everywhere
template<bool unneeded = true>
std::vector<int> GetAssignment(
      const std::vector<std::vector<REAL> >& cost,
      const std::vector<std::vector<INDEX> >& graph)
{
   std::vector<std::pair<INDEX,INDEX> > nodeCapacity(graph.size(), std::pair<INDEX,INDEX>(0,1));
   INDEX inverseGraphSize = 0;
   for(auto it : graph) {
      inverseGraphSize = std::max(inverseGraphSize, *std::max_element(it.begin(), it.end()));
   }
   std::vector<std::pair<INDEX,INDEX> > inverseNodeCapacity(inverseGraphSize, std::pair<INDEX,INDEX>(0,1));
   return GetAssignment(cost,graph,nodeCapacity,inverseNodeCapacity);
}

} // end namespace LP_MP

#endif // LP_MP_LINEAR_ASSIGNMENT_HXX
