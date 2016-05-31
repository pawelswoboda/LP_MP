#ifndef LP_MP_UNION_FIND_HXX
#define LP_MP_UNION_FIND_HXX

namespace LP_MP {
class UnionFind {
   INDEX *id, cnt, *sz;
public:
   // Create an empty union find data structure with N isolated sets.
   UnionFind(INDEX N) {
      cnt = N;
      id = new INDEX[N];
      sz = new INDEX[N];
      for(INDEX i=0; i<N; i++) {
         id[i] = i;
         sz[i] = 1;
      }
   }
   ~UnionFind() {
      delete [] id;
      delete [] sz;
   }
   // Return the id of component corresponding to object p.
   INDEX find(INDEX p) {
      INDEX root = p;
      while (root != id[root])
         root = id[root];
      while (p != root) {
         INDEX newp = id[p];
         id[p] = root;
         p = newp;
      }
      return root;
   }
   // Replace sets containing x and y with their union.
   void merge(INDEX x, INDEX y) {
      INDEX i = find(x);
      INDEX j = find(y);
      if(i == j) return;

      // make smaller root point to larger one
      if(sz[i] < sz[j])	{ 
         id[i] = j; 
         sz[j] += sz[i]; 
      } else	{ 
         id[j] = i; 
         sz[i] += sz[j]; 
      }
      cnt--;
   }
   // Are objects x and y in the same set?
   bool connected(INDEX x, INDEX y) {
      return find(x) == find(y);
   }
   // Return the number of disjoint sets.
   INDEX count() {
      return cnt;
   }
};

}; // end namespace LP_MP
#endif // LP_MP_UNION_FIND_HXX
