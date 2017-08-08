#ifndef LP_MP_UNION_FIND_HXX
#define LP_MP_UNION_FIND_HXX

namespace LP_MP {
class UnionFind {
   INDEX *id, cnt, *sz, N; // it is not necessary to hold sz!
public:
   // Create an empty union find data structure with N isolated sets.
   UnionFind(INDEX _N) : N(_N) {
      id = new INDEX[2*N];
      sz = id + N;
      reset();
   }
   ~UnionFind() {
      delete [] id;
   }
   void reset() {
      cnt = N;
      for(INDEX i=0; i<N; ++i) { id[i] = i; }
      for(INDEX i=0; i<N; ++i) { sz[i] = 1; }
   }
   // Return the id of component corresponding to object p.
   INDEX find(INDEX p) {
      assert(p < N);
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

   INDEX thread_safe_find(INDEX p) const {
      INDEX root = p;
      while (root != id[root])
         root = id[root];
      return root;
   }
   bool thread_safe_connected(INDEX x, INDEX y) const {
      return thread_safe_find(x) == thread_safe_find(y);
   }
   // Return the number of disjoint sets.
   INDEX count() {
      return cnt;
   }

   std::vector<INDEX> get_contiguous_ids()
   {
      std::vector<INDEX> contiguous_ids(N);
      std::vector<INDEX> id_mapping(N, std::numeric_limits<REAL>::max());
      //INDEX* id_mapping = new INDEX[N];
      //std::fill(contiguous_ids.begin(), contiguous_idx.end(), std::numeric_limits<INDEX>::max());
      for(INDEX i=0; i<N; ++i) {
         INDEX d = find(i);
         id_mapping[d] = 1; 
      }
      INDEX next_id = 0;
      for(INDEX d=0; d<N; ++d) {
         if(id_mapping[d] == 1) {
            id_mapping[d] = next_id;
            ++next_id;
         }
      }

      for(INDEX i=0; i<N; ++i) {
         INDEX d = find(i);
         assert(id_mapping[d] != std::numeric_limits<INDEX>::max());
         contiguous_ids[i] = id_mapping[d];
      }
      return std::move(id_mapping);
   }
};

}; // end namespace LP_MP

#endif // LP_MP_UNION_FIND_HXX
