// Compute topological sorting of a DAG
#include <iostream>
#include <list>
#include <stack>
#include <queue> 
#include <assert.h>
#include "help_functions.hxx"
#include "config.hxx"

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

// add directed edge from v to w
void Graph::addEdge(INDEX v, INDEX w)
{
   assert(v<V && w<V);
   adj[v].push_back(w); 
}

std::vector<INDEX> Graph::topologicalSort()
{
   constexpr unsigned char notMarked = 0;
   constexpr unsigned char tempMarked = 1;
   constexpr unsigned char permMarked = 2;
   
   std::vector<unsigned char> mark(V,notMarked); 
   std::stack<std::pair<INDEX,decltype(adj[0].begin())> > dfs;
   std::vector<INDEX> postOrder;

   for(INDEX i=0;i<V;i++){
      if(mark[i] == notMarked) {
         dfs.push(std::make_pair(i,adj[i].begin()));
         while(!dfs.empty()){
            INDEX node = dfs.top().first;
            auto& it = dfs.top().second;
            //if(mark[node] == tempMarked) {
            //   throw std::runtime_error("graph not a dag");
            //}
            if(it == adj[node].begin()) { // first visit
               assert(mark[node] == notMarked);
               mark[node] = permMarked;
            }
            while(it != adj[node].end() && mark[*it] != notMarked) {
               ++it;
            }
            if(it != adj[node].end()) {
               INDEX next = *it;
               ++it;
               dfs.top().second = it;
               dfs.push(std::make_pair(next,adj[next].begin()));
            } else {
               dfs.pop();
               mark[node] = permMarked;
               postOrder.push_back(node);
            }
         }
      }
   }

   assert(postOrder.size() == INDEX(V));
   assert(LP_MP::HasUniqueValues(postOrder));
   std::reverse(postOrder.begin(),postOrder.end());
   return std::move(postOrder);
}

} // end namespace Topological_Sort
} // end namespace LP_MP
