#ifndef MINCONV_HXX
#define MINCONV_HXX

#include <limits>
#include <cmath>
#include <functional>
#include <queue>
#include <stdexcept>
#include <algorithm>

namespace LP_MP {
  namespace discrete_tomo{

    template<class Value,class Index = int>
    class MinConv
    {
    public:

      template<class T1, class T2>
      MinConv(T1 a, T2 b, Index n, Index m,Index t);

      void init(Index,Index,Index);
      Value getConv(Index k){ return c_[k]; };
      Value getMin(){ return minimum_; };
      Index getIdxA(Index k){ assert(0 <= k);assert(k < outA_.size()); return outA_[k]; };
      Index getIdxB(Index k){ assert(0 <= k);assert(k < outB_.size()); return outB_[k]; };

      template<class T1,class T2,class T3>
      void CalcConv(T1 op,T2 a,T3 b,bool onlyMin = false);
    
    private:
      std::vector<Index> idxa_;
      std::vector<Index> idxb_;
      std::vector<Index> outA_;
      std::vector<Index> outB_;
      std::vector<Value> c_;
      std::vector<Index> cp_;
      Value minimum_ = std::numeric_limits<Value>::infinity();
      Value inf = std::numeric_limits<Value>::infinity();

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
    template<class T1,class T2>
    MinConv<Value,Index>::MinConv(T1 a, T2 b, Index n, Index m,Index t){
      if( n < 1 || m < 1 || t < 1){ 
	throw std::runtime_error("m,n and t should > 0!");  }
      sort(a,n,idxa_);
      sort(b,m,idxb_);
      c_.resize(t+1); cp_.resize(t+1);
      outA_.resize(t+1); outB_.resize(t+1);
      std::fill(c_.begin(),c_.end(),inf);
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
	return a(idxa_[i]) + b(idxb_[j]);
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
      std::priority_queue<indices,std::vector<indices>,decltype(compare) > queue(compare);
      std::vector<Index> cover(m+n,0);

      queue.push(indices(0,0));
      cover[0] = 1; cover[n] = 1;

      /* START define cover operations */
      auto remCover = [&](Index i,Index j){
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
	}
      };
          
      auto addCoverI = [&](Index i,Index j){
	if( cover[i+1] == 0 && cover[n+j] == 0 ){
	  Index inc = 1;
	  while( i+inc < n ){
	    Index idx = op(idxa_[i+inc],idxb_[j]);//idxa_[i+inc]+idxb_[j];
	    if( cp_[idx] == 0 ){
	      addCover(i+inc,j);
	      break;
	    }
	    inc++;
	  }
	}
      };

      auto addCoverJ = [&](Index i,Index j){
	if( cover[i] == 0 && cover[n+j+1] == 0 ){
	  Index inc = 1;
	  while( j+inc < m ){
	    Index idx = op(idxa_[i],idxb_[j+inc]);//idxa_[i]+idxb_[j+inc];
	    if( cp_[idx] == 0 ){
	      addCover(i,j+inc);
	      break;
	    }
	    inc++;
	  }
	}
      };
      /* END define cover operations */
    
      while( open > 0 && !queue.empty() ){
	indices it = queue.top();
	Index i=it.i, j=it.j;
	Value minV = M(i,j);
	Index mink = op(idxa_[i],idxb_[j]);//idxa_[i]+idxb_[j]
       	queue.pop();

	assert(0 <= mink); assert( mink < cp_.size() ); // operation seems to be wrong
	
	if( c_[mink] > minV || cp_[mink] == 0 ){
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
      
    }
  }
}

#endif
