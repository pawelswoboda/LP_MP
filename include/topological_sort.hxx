#ifndef LP_MP_TOPOLOGICAL_SORT
#define LP_MP_TOPOLOGICAL_SORT

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
    bool sorting_valid(const std::vector<INDEX>& ordering) const;

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

bool Graph::sorting_valid(const std::vector<INDEX>& ordering) const
{
  std::vector<INDEX> inverse_ordering(ordering.size());
  for(INDEX i=0; i<ordering.size(); ++i) {
    inverse_ordering[ordering[i]] = i; 
  }

   // check validity of sorting
   for(INDEX i=0; i<adj.size(); ++i) {
     for(const INDEX j : adj[i]) {
       assert(inverse_ordering[i] != inverse_ordering[j]);
       if(inverse_ordering[i] > inverse_ordering[j]) {
         return false; 
       }
     }
   }
   return true; 
}


std::vector<INDEX> Graph::topologicalSort()
{
  if(debug()) {
    std::cout << "sort " << adj.size() << " elements subject to " << std::accumulate(adj.begin(), adj.end(), INDEX(0), [](INDEX s, auto& a) { return s + a.size(); }) << " ordering constraints\n";
  }
   /*
   std::vector<bool> visited(V,false);
   std::vector<INDEX> order;
   std::deque<SIGNED_INDEX> stack; // change to stack

   for(SIGNED_INDEX i=0; i<V; ++i) {
      if(visited[i] == false) {
         visited[i] = true;
         stack.push_back(i);
         while(!stack.empty()) {
            auto v = stack.back();
            stack.pop_back();
            if(v>=0) {
               stack.push_back(-v-1);
               for(SIGNED_INDEX w : adj[v]) {
                  // if graph is acyclic, then w may not already be
                  if(!visited[w]) {
                     visited[w] = true;
                     stack.push_back(w); 
                  } 
               }
            } else {
               order.push_back(-v-1); 
            } 
         }
      }
   }
   assert(order.size() == INDEX(V));
   assert(LP_MP::HasUniqueValues(order));
   std::reverse(order.begin(),order.end());
   return std::move(order);
   */

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
            if(mark[node] == tempMarked) {
               throw std::runtime_error("graph not a dag");
            }
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

   assert(sorting_valid(postOrder));

   return std::move(postOrder);
}

} // end namespace Topological_Sort
} // end namespace LP_MP

#endif // LP_MP_TOPOLOGICAL_SORT
