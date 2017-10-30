#ifndef MINCONV_HXX
#define MINCONV_HXX

#include <limits>
#include <cmath>
#include <functional>
#include <queue>
#include <stdexcept>
#include <algorithm>
#include <tuple>
#include <cassert>
#include "vector.hxx"

// do zrobienia: possibly make own project out of min convolution, but then used vector class and associated memory allocator also needs to be a separate project

namespace LP_MP {
   namespace discrete_tomo{

      // functions naively computing min convolutions
      // possibly encapsulate in class as for efficient min conv

      template<typename ITERATOR_1, typename ITERATOR_2>
      void min_conv_naive(ITERATOR_1 a_begin, ITERATOR_1 a_end, ITERATOR_2 b_begin, ITERATOR_2 b_end, vector<REAL>& result)
      {
         const INDEX a_size = std::distance(a_begin, a_end);
         assert(a_size > 1);
         const INDEX b_size = std::distance(b_begin, b_end);
         assert(b_size > 1);
         assert(result.size() <= a_size + b_size - 1);
         std::fill(result.begin(), result.end(), std::numeric_limits<REAL>::infinity());
         for(INDEX i=0; i<std::min(result.size(), a_size); ++i) {
            //assert(false); // bounds below correct?
            for(INDEX j=0; j<std::min(b_size, result.size() - i); ++j) {
               result[i+j] = std::min(result[i+j], a_begin[i] + b_begin[j]); 
            }
         }
      }

      template<typename ITERATOR_1, typename ITERATOR_2>
      vector<REAL> min_conv_naive(ITERATOR_1 a_begin, ITERATOR_1 a_end, ITERATOR_2 b_begin, ITERATOR_2 b_end, const INDEX sum)
      {
         vector<REAL> result(sum);
         min_conv_naive(a_begin, a_end, b_begin, b_end, result);
         return result;
      }

      template<typename ITERATOR_1, typename ITERATOR_2>
      vector<REAL> min_conv_naive(ITERATOR_1 a_begin, ITERATOR_1 a_end, ITERATOR_2 b_begin, ITERATOR_2 b_end)
      {
         const INDEX a_size = std::distance(a_begin, a_end);
         const INDEX b_size = std::distance(b_begin, b_end);
         return min_conv_naive(a_begin, a_end, b_begin, b_end, a_size + b_size - 1); 
      }

      using arg_min_conv_vector = std::vector<std::array<unsigned char,2>>;

      using arg_min_conv_type = vector<std::array<INDEX,2>>;

      // returns the min convolution and extends the given vector with the indices that result in the given min convolution
      template<typename ITERATOR_1, typename ITERATOR_2>
      std::tuple<vector<REAL>, arg_min_conv_type> arg_min_conv_naive( ITERATOR_1 a_begin, ITERATOR_1 a_end, ITERATOR_2 b_begin, ITERATOR_2 b_end, const INDEX min_conv_size)
      {
         const INDEX a_size = std::distance(a_begin, a_end);
         assert(a_size > 1);
         const INDEX b_size = std::distance(b_begin, b_end);
         assert(b_size > 1);
         vector<REAL> result(min_conv_size, std::numeric_limits<REAL>::infinity());
         vector<std::array<INDEX,2>> arg_result(min_conv_size);
         assert(result.size() <= a_size + b_size - 1);
         for(INDEX i=0; i<std::min(result.size(), a_size); ++i) {
            for(INDEX j=0; j<std::min(b_size, result.size() - i); ++j) {
               const REAL cur_val = a_begin[i] + b_begin[j];
               if(cur_val <= result[i+j]) {
                  result[i+j] = a_begin[i] + b_begin[j];
                  arg_result[i+j] = std::array<INDEX,2>({i,j});
               }
            }
         } 
         return std::make_tuple(std::move(result), std::move(arg_result));
      }





      template<typename ITERATOR_1, typename ITERATOR_2>
      REAL min_sum(ITERATOR_1 a_begin, ITERATOR_1 a_end, ITERATOR_2 b_begin, ITERATOR_2 b_end, const INDEX sum)
      {
         REAL val = std::numeric_limits<REAL>::infinity();

         const INDEX a_size = std::distance(a_begin, a_end);
         assert(a_size > 1);
         const INDEX b_size = std::distance(b_begin, b_end);
         assert(b_size > 1);

         assert(sum <= a_size-1 + b_size-1);

         for(INDEX i=std::max(0, int(b_size)-int(sum)-1); i<std::min(sum+1, a_size); ++i) {
            val = std::min(val, a_begin[i] + b_begin[sum-i]);
         }
         return val;
      }

      // find out the two indices (left,right) whose sum is minimal
      template<typename ITERATOR_1, typename ITERATOR_2>
      std::tuple<INDEX,INDEX> arg_min_sum(ITERATOR_1 a_begin, ITERATOR_1 a_end, ITERATOR_2 b_begin, ITERATOR_2 b_end, const INDEX sum)
      {
         REAL val = std::numeric_limits<REAL>::infinity();
         INDEX a = std::numeric_limits<INDEX>::max();
         INDEX b = std::numeric_limits<INDEX>::max();

         const INDEX a_size = std::distance(a_begin, a_end);
         assert(a_size > 1);
         const INDEX b_size = std::distance(b_begin, b_end);
         assert(b_size > 1);

         assert(sum <= a_size-1 + b_size-1);

         for(INDEX i=0; i<std::min(a_size, sum+1); ++i) {
            const INDEX j = sum-i;
            if(j < b_size) { // more efficient is possible
               const REAL cur_val = a_begin[i] + b_begin[sum-i];
               if(cur_val <= val) {
                  a = i;
                  b = sum - i;
                  val = cur_val;
               }
            }
         }
         assert(a < std::numeric_limits<INDEX>::max() && b < std::numeric_limits<INDEX>::max());
         return std::make_tuple(a,b);;
      }

      /*
      // matrix has dimensions no_variables/cardinality
      // we want to compute the minimum cost configuration so that \sum_{i} x_i = sum

      // compute cost of each possible value of sub sum
      // iteratively equipartition cost matrix into left and right part and compute via min convolution cost of sub sums
      vector<REAL> min_sub_sum_naive(const matrix<REAL>& m, const INDEX start_col, const INDEX end_col)
      {
         assert(end_col < m.dim2() && start_col < end_col);
         const INDEX n = end_col - start_col;

         // allocate result already here, so intermediate vectors can be released immediately (remember: we use a stack allocator)
         vector<REAL> result((m.dim1()-1) * n + 1, std::numeric_limits<REAL>::infinity());

         if(n == 2) {
            //return discrete_tomo::min_conv_naive(&m(start_col,0), &m(start_col, m.dim2()-1), &m(start_col+1, 0), &m(start_col+1, m.dim2()-1));
            auto slice_1 = m.slice_left(start_col);
            auto slice_2 = m.slice_left(start_col+1);
            discrete_tomo::min_conv_naive(slice_1.begin(), slice_1.end(), slice_2.begin(), slice_2.end(), result);
         } else if(n == 3) {
            const auto left_sum = sub_sum(m, start_col+1, end_col);
            auto slice = m.slice_left(start_col);
            discrete_tomo::min_conv_naive(left_sum.begin(), left_sum.end(), slice.begin(), slice.end(), result); 
         } else {
            assert(n > 3);
            const auto left_sum = sub_sum(m, start_col, start_col + n/2);
            const auto right_sum = sub_sum(m, start_col + n/2, end_col);
            discrete_tomo::min_conv_naive(left_sum.begin(), left_sum.end(), right_sum.begin(), right_sum.end(), result);
         }

         return std::move(result);
      }

      REAL min_sum_naive(const matrix<REAL>& m, const INDEX sum)
      {
         auto left = min_sub_sum_naive(m, 0, m.dim2()/2);
         auto right = min_sub_sum_naive(m, m.dim2()/2, m.dim2());
         return min_sum_naive(left.begin(), left.end(), right.begin(), right.end(), sum); 
      }
      */

      // class for more efficient computation of min convolution. Algorithm by Bauschke et al.  
      template<typename ITERATOR1, typename ITERATOR2>
         class MinConv
         {
            public:

               template<class T1, class T2>
                  MinConv(T1 a, T2 b, Index n, Index m,Index t);

               void init(Index,Index,Index);
               Value getConv(Index k) const { return c_[k]; };
               Value getMin() const { return minimum_; };
               Index getIdxA(Index k) const { assert(0 <= k);assert(k < outA_.size()); return outA_[k]; };
               Index getIdxB(Index k) const { assert(0 <= k);assert(k < outB_.size()); return outB_[k]; };

               template<class T1,class T2,class T3>
                  void CalcConv(T1 op,T2 a,T3 b,bool onlyMin = false);

            private:
               vector<Index> idxa_;
               vector<Index> idxb_;
               vector<Index> outA_;
               vector<Index> outB_;
               vector<REAL> c_;
               vector<Index> cp_;
               Value minimum_ = std::numeric_limits<Value>::infinity();

               template<class T>
                  void sort(const T&,Index,std::vector<Index>&);
         };

      // idx stores the sorted indices according to T
      template<class Value,class Index>
         template<class T>
         void MinConv<Value,Index>::sort(const T &V,Index n,std::vector<Index> &idx){
            idx.resize(n);
            std::iota(idx.begin(),idx.end(),0);
            auto compare = [&](Index i,Index j){ return (V(i) < V(j)); };
            std::sort(idx.begin(),idx.end(),compare);
         }

      template<class Value,class Index>
         template<class VEC1,class VEC2>
         MinConv<Value,Index>::MinConv(VEC1 a, VEC2 b, Index n, Index m,Index t)
         : c_(t+1)
         {
            if( n < 1 || m < 1 || t < 1){ 
               throw std::runtime_error("m,n and t should > 0!");  }
            //assert(n-1+m-1 >= t-1);
            sort(a,n,idxa_);
            sort(b,m,idxb_);
            //c_.resize(t+1); 
            cp_.resize(t+1);
            outA_.resize(t+1); outB_.resize(t+1);
            std::fill(c_.begin(),c_.end(),std::numeric_limits<Value>::infinity());
            std::fill(cp_.begin(),cp_.end(),0);
         };

      /* 
         The input functions defines the operation for the indices
         k = i o j (e.g. i + j or i - j). 
         if not 0 <= i o j < t, then return t! 
         */
      template<class Value,class Index>
         template<class T1,class T2,class T3>
         void
         MinConv<Value,Index>::CalcConv(T1 op,T2 a,T3 b,bool onlyMin)
         {
            // access function to the sorted matrix      
            auto M = [&](Index i,Index j){
               //assert( a(idxa_[i]) < std::numeric_limits<Value>::max() || b(idxb_[j]) > -std::numeric_limits<Value>::max() );
               //assert( a(idxa_[i]) > -std::numeric_limits<Value>::max() || b(idxb_[j]) < std::numeric_limits<Value>::max() );
               //assert( b(idxb_[j]) < std::numeric_limits<Value>::max() || a(idxa_[i]) > -std::numeric_limits<Value>::max() );
               //assert( b(idxb_[j]) > -std::numeric_limits<Value>::max() || a(idxa_[i]) < std::numeric_limits<Value>::max() );
               assert( !std::isnan(a(idxa_[i])) );
               assert( !std::isnan(b(idxb_[j])) );
               // do zrobienia: below two check are not needed anymore
               Value ai = ( a(idxa_[i]) > -std::numeric_limits<Value>::max() ) ?  a(idxa_[i]) : std::numeric_limits<Value>::max();
               Value bj = ( b(idxb_[j]) > -std::numeric_limits<Value>::max() ) ?  b(idxb_[j]) : std::numeric_limits<Value>::max();
               return ai + bj;
            };

            Index n = idxa_.size();
            Index m = idxb_.size();
            Index open = c_.size();

            struct indices{
               indices(Index ii, Index jj) : i(ii),j(jj) {}; 
               Index i = 0; Index j = 0;
            };

            auto compare = [&](indices x,indices y){
               return M(x.i,x.j) > M(y.i,y.j);
            };
            // do zrobienia: use static memory pool allocator and better datastructore for indices
            std::priority_queue<indices,std::vector<indices>,decltype(compare) > queue(compare);
            std::vector<Index> cover(m+n,0);

            queue.push(indices(0,0));
            cover[0] = 1; cover[n] = 1;

            /* START define cover operations */
            auto remCover = [&](Index i,Index j){
               //std::cout << "RemCover: (" << i;
               //std::cout << "," << j << ") = " << op(idxa_[i],idxb_[j]);
               //std::cout << std::endl;
               if( cover[i] > 0 )
               { cover[i]--; }
               if( cover[n+j] > 0 )
               { cover[n+j]--; }
            };
            auto addCover = [&](Index i,Index j){
               if( cover[i] == 0 && cover[n+j] == 0 ){
                  cover[i]++; 
                  cover[n+j]++;
                  queue.push(indices(i,j));
                  //std::cout << "AddCover: (" << i;
                  //std::cout << "," << j << ") = " << op(idxa_[i],idxb_[j]);
                  //std::cout << std::endl;
               }
            };

            auto addCoverI = [&](Index i,Index j){
               if( cover[i+1] == 0 && cover[n+j] == 0 ){
                  Index inc = 1;
                  while( i+inc < n ){
                     Index idx = op(idxa_[i+inc],idxb_[j]);//idxa_[i+inc]+idxb_[j];
                     idx = std::min(idx, Index(c_.size()-1));
                     assert(idx<cp_.size());
                     if( cp_[idx] == 0 ){
                        addCover(i+inc,j);
                        break;
                     }
                     if( cover[i+inc] > 0 ){ break; }
                     inc++;
                  }
               }
            };

            auto addCoverJ = [&](Index i,Index j){
               if( cover[i] == 0 && cover[n+j+1] == 0 ){
                  Index inc = 1;
                  while( j+inc < m ){
                     Index idx = op(idxa_[i],idxb_[j+inc]);//idxa_[i]+idxb_[j+inc];
                     idx = std::min(idx, Index(c_.size()-1));
                     assert(idx<cp_.size());
                     if( cp_[idx] == 0 ){
                        addCover(i,j+inc);
                        break;
                     }
                     if( cover[n+j+inc] > 0 ){ break; }
                     inc++;
                  }
               }
            };
            /* END define cover operations */

            //std::cout << "START" << std::endl;
            /*
               for(Index i=0;i<n;i++){
               for(Index j=0;j<m;j++){
               printf("%03d .. ",op(idxa_[i],idxb_[j]));  
               }
               printf("\n");
               }
               */      
            while( open > 0 && !queue.empty() ){
               indices it = queue.top();
               Index i=it.i, j=it.j;
               assert( i < n ); assert( j < m );

               Value minV = M(i,j);
               assert(!std::isnan(minV));
               Index mink = op(idxa_[i],idxb_[j]);//idxa_[i]+idxb_[j]
               mink = std::min(mink, Index(c_.size()-1)); // all out of bounds go to entry last entry of c_, i.e. c_.size()-1
               queue.pop();

               assert(0 <= mink); assert( mink < cp_.size() ); // operation seems to be wrong

               if( c_[mink] > minV || cp_[mink] == 0 ){
                  assert(cp_[mink] == 0);

                  c_[mink] = minV; open--;

                  if( mink != c_.size()-1 && minV < minimum_ ){
                     minimum_ = minV;
                     if( onlyMin ){ break; }
                  }

                  outA_[mink] = idxa_[i];
                  outB_[mink] = idxb_[j];
                  cp_[mink] = 1;
               }

               remCover(i,j);

               if( i+1 < n && j+1 < m )
               {
                  if( cover[i+1] == 0 && cover[n+j] == 0 
                        && cover[i] == 0 && cover[n+j+1] == 0 ){
                     addCover(i,j+1);
                     addCoverI(i,j);
                  }
                  else if( cover[i+1] == 0 && cover[n+j] == 0 ){
                     addCoverI(i,j);
                  }
                  else if( cover[i] == 0 && cover[n+j+1] == 0 ){
                     addCoverJ(i,j);
                  }
               }
               else if( i+1 < n )
               {
                  addCoverI(i,j);
               }
               else if( j+1 < m )
               {
                  addCoverJ(i,j);
               }
            }
            //assert(open == 1 || open == 0);
         }

      static const INDEX min_conv_threshold = 500; // when more elements than indicated by threshold need to be computed, use heuristic, otherwise use naive implementation.

      // automatically choose between naive and efficient version of min convolution
      template<typename ITERATOR_1, typename ITERATOR_2>
      void min_conv(ITERATOR_1 a_begin, ITERATOR_1 a_end, ITERATOR_2 b_begin, ITERATOR_2 b_end, vector<REAL>& result)
      {
         if(result.size() < min_conv_threshold) {
            return min_conv_naive(a_begin, a_end, b_begin, b_end, result);
         } else {
            assert(false);
         }
      }

      template<typename ITERATOR_1, typename ITERATOR_2>
      vector<REAL> min_conv(ITERATOR_1 a_begin, ITERATOR_1 a_end, ITERATOR_2 b_begin, ITERATOR_2 b_end, const INDEX sum)
      {
         if(sum < min_conv_threshold) {
            return min_conv_naive(a_begin, a_end, b_begin, b_end, sum);
         } else {
            assert(false);
         }
      }

      template<typename ITERATOR_1, typename ITERATOR_2>
      vector<REAL> min_conv(ITERATOR_1 a_begin, ITERATOR_1 a_end, ITERATOR_2 b_begin, ITERATOR_2 b_end)
      {
         if(std::distance(a_begin, a_end) + std::distance(b_begin, b_end) < min_conv_threshold) {
            return min_conv_naive(a_begin, a_end, b_begin, b_end);
         } else {
            assert(false);
         }
      }

      template<typename ITERATOR_1, typename ITERATOR_2>
      std::tuple<vector<REAL>, arg_min_conv_type> arg_min_conv( ITERATOR_1 a_begin, ITERATOR_1 a_end, ITERATOR_2 b_begin, ITERATOR_2 b_end, const INDEX min_conv_size)
      {
         if(min_conv_size < min_conv_threshold) {
            return min_conv_naive(a_begin, a_end, b_begin, b_end, min_conv_size);
         } else {
            assert(false);
         }
      }



   } // end namespace discrete_tomo
} // end namespace LP_MP

#endif
