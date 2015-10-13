// Compute topological sorting of a DAG
#include <iostream>
#include <list>
#include <stack>
#include <queue> 
#include <assert.h>
#include "help_functions.hxx"
#include "instances.inc"

namespace LP_MP {
namespace Topological_Sort {

// do zrobienia: detect, when graph is not DAG
class Graph
{
    INDEX V;    // number of vertices
    std::vector<std::list<INDEX> > adj;
public:
    inline Graph(INDEX V);   
    inline void addEdge(INDEX v, INDEX w);
    inline std::vector<INDEX> topologicalSort();
};
 
Graph::Graph(INDEX V)
   : adj(V)
{
   this->V = V;
}

void Graph::addEdge(INDEX v, INDEX w)
{
   assert(v<V && w<V);
   adj[v].push_back(w); 
}

std::vector<INDEX> Graph::topologicalSort()
{
   /*
   std::stack<INDEX> nodeList; 
   std::vector<INDEX> pred(V,-1);
   std::vector<bool> visited(V,false);
   std::vector<INDEX> finishingTime(V,0);

   INDEX f = 1;
   for(INDEX i=0; i<V; i++) {
      if(visited[i] == false) {
         nodeList.push(i);
         while(!nodeList.empty()) {
            const INDEX node = nodeList.top();
            nodeList.pop();
            std::list<INDEX> newVec = adj[node]; //vector of neighbors
            for(auto it=newVec.begin(); it!=newVec.end(); ++it) {
               if(visited[*it] == false) {
                  //path[depth+1] = node;
                  pred[*it] = node;
                  visited[*it] = true;
                  finishingTime[*it] = f++;
                  nodeList.push(*it);
               }
            }
         }
      }
   }

   // now sort with respect to finishing times. This will determine the order
   std::vector<INDEX> nodeOrder(V);
   for(INDEX i=0; i<nodeOrder.size(); i++) nodeOrder[i] = i;
   std::sort(nodeOrder.begin(), nodeOrder.end(), [&](const INDEX i, const INDEX j)->bool {return finishingTime[i] < finishingTime[j];});
   return nodeOrder;
   */



   std::vector<bool> visited(V,false);
   std::stack<std::pair<bool,INDEX> > dfs;
   std::stack<INDEX> postOrder;
   std::list<INDEX> newVec;
   std::list<INDEX>::iterator it;
   for(INDEX i=0;i<V;i++){
      if(visited[i]==false){
         dfs.push(std::make_pair(false,i));
      }   
      while(!dfs.empty()){
         std::pair<bool,INDEX> node=dfs.top();
         dfs.pop();
         if (node.first) {
            postOrder.push(node.second);
         } else {
            visited[node.second]=true;
            dfs.push(std::make_pair(true, node.second));
            newVec=adj[node.second]; //vector of neighbors
            for(it=newVec.begin();it!=newVec.end();it++){
               INDEX son=*it;
               if(visited[son]==false){
                  dfs.push(std::make_pair(false, son));
               }
            }
         }
      }
   }
   std::vector<INDEX> sorted_vertices(0);
   sorted_vertices.reserve(V);
   while (postOrder.empty() == false)
   {
      //std::cout << postOrder.top() << " ";
      sorted_vertices.push_back(postOrder.top());
      postOrder.pop();
   }
   //std::cout << std::endl;
   std::cout << sorted_vertices.size() << ", " << V << std::endl;
   assert(LP_MP::HasUniqueValues(sorted_vertices));
   assert(sorted_vertices.size() == INDEX(V));
   return sorted_vertices;
}

} // end namespace Topological_Sort
} // end namespace LP_MP
