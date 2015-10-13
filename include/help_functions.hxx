#ifndef LP_MP_HELP_FUNCTIONS
#define LP_MP_HELP_FUNCTIONS

#include <map>
#include <set>
#include <vector>
#include <limits>
#include <numeric>
#include <algorithm>
#include <assert.h>

namespace LP_MP {

inline int binary_compl(const int i)
{
   assert(i == 0 && i == 1);
   return 1-i;
}

template<class T>
int find_index(const T i, const std::vector<T>& vec)
{
   const int index = find(vec.begin(), vec.end(), i) - vec.begin();
   assert(index <= vec.size());
   return index;
}

template<class F, class I>
void BuildIndexMaps(const std::vector<F>& f, std::map<F,I>& elemToIndex, std::map<I,F>& indexToElem)
{
   // do zrobienia: - reserve space, 
   //               - possibly use hash_map for speed
   for(size_t i=0; i<f.size(); i++) {
      elemToIndex.insert(std::make_pair(f[i],i));
      indexToElem.insert(std::make_pair(i,f[i]));
   }
}

template<class T>
std::vector<T> GetSubVector(const std::vector<T>& v, const std::vector<size_t>& inds)
{
   std::vector<T> subVec(inds.size());
   for(size_t i=0; i<inds.size(); i++) {
      subVec[i] = v[inds[i]];
   }
   return subVec;
}

template<class T>
bool HasUniqueValues(const std::vector<T>& v)
{
   std::set<T> values;
   for(size_t i=0; i<v.size(); i++) {
      if(values.find( v[i] ) != values.end()) return false;
      else values.insert( v[i] );
   }
   return true;
}

// return indices belonging to the three smallest entries
template<class T>
std::tuple<size_t,size_t,size_t> MinThreeIndices(const std::vector<T>& v)
{
   assert(v.size() > 2);
   std::vector<size_t> ind(v.size());
   for(size_t i=0; i<ind.size(); i++)  ind[i] = i;
   // do zrobienia: use partial_sort here
   std::sort(ind.begin(), ind.end(), [&](const size_t i, const size_t j)->bool {return v[i] < v[j];} );
   return std::tuple<size_t,size_t,size_t>(ind[0], ind[1], ind[2]);
}

template<class T, class A>
std::pair<size_t,size_t> MinIndices(const A& v)
{
   assert(v.size() > 1);
   T min_val = std::numeric_limits<T>::max();
   size_t min_index, second_min_index;
   T second_min_val = std::numeric_limits<T>::max();
   for(size_t i=0; i<v.size(); i++) {
      if(v[i] <= min_val) {
         min_val = v[i];
         min_index = i;
      }
   }
   for(size_t i=0; i<v.size(); i++) {
      if(v[i] < second_min_val && i != min_index) {
         second_min_val = v[i];
         second_min_index = i;
      }
   }
   assert(min_index != second_min_index);
   return std::pair<T,T>(min_index, second_min_index);
}

template<class T, class A>
std::pair<T,T> SmallestValues(const A& v)
{
   std::pair<size_t,size_t> min_indices = MinIndices<T>(v);
   return std::pair<T,T>(v[min_indices.first], v[min_indices.second]);
}


template<class T>
void NormalizeVector(std::vector<T>& v)
{
   //const T mean = std::accumulate(v.begin(), v.end(),0)/v.size();
   return; // do zrobienia: why is it not applicable in marg_message?
   const T mean = *std::min_element(v.begin(), v.end());
   for(size_t i=0; i<v.size(); i++) {
      v[i] -= mean;
   }
}

} // end namespace LP_MP

#endif // LP_MP_HELP_FUNCTIONS
