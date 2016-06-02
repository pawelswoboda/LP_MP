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

    template<class T,class Value,class Index = int>
    class MinConv
    {
    public:

      MinConv(T &a,T &b, Index n, Index m,Index t) : a_(a),b_(b)
      { init(n,m,t); };
      MinConv(T &a,T &b) : a_(a),b_(b){};

      void init(Index,Index,Index);
      Value getConv(Index k){ return c_[k]; };
      Value getMin(){ return minimum_; };
      Index getIdxA(Index k){ return outA_[k]; };
      Index getIdxB(Index k){ return outB_[k]; };

      void CalcConv(std::function<Index(Index,Index)>,bool onlyMin = false);
    
    private:
      T &a_;
      T &b_;
      std::vector<Index> idxa_;
      std::vector<Index> idxb_;
      std::vector<Index> outA_;
      std::vector<Index> outB_;
      std::vector<Value> c_;
      std::vector<Index> cp_;
      Value minimum_ = std::numeric_limits<Value>::infinity();
      Value inf = std::numeric_limits<Value>::infinity();
    
      void sort(const T&,Index,std::vector<Index>&);
    };

    // idx stores the sorted indices according to T
    template<class T,class Value,class Index>
    void MinConv<T,Value,Index>::sort(const T &V,Index n,std::vector<Index> &idx){
      idx.resize(n);
      std::iota(idx.begin(),idx.end(),0);
      auto compare = [&](Index i,Index j){ return (V(i) < V(j)); };
std::sort(idx.begin(),idx.end(),compare);
  }

  // initialize the size of a,b and c
  template<class T,class Value,class Index>
  void MinConv<T,Value,Index>::init(Index n,Index m,Index t){
    if( n < 1 || m < 1 || t < 1){ 
      throw std::runtime_error("m,n and t should > 0!");  }
    sort(a_,n,idxa_);
    sort(b_,m,idxb_);
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
  template<class T,class Value,class Index>
  void
  MinConv<T,Value,Index>::CalcConv(std::function<Index(Index,Index)> op,bool onlyMin)
  {
    // access function to the sorted matrix      
    auto M = [this](Index i,Index j){
      return a_(idxa_[i]) + b_(idxb_[j]);
    };

    Index n = idxa_.size();
    Index m = idxb_.size();
    Index open = c_.size();

    struct indices{
      indices(Index ii, Index jj) : i(ii),j(jj) {}; 
      Index i = 0; Index j = 0;
    };
      
    std::function<bool(indices,indices)> compare = [&](indices a,indices b){
      return M(a.i,a.j) > M(b.i,b.j);
    };
    std::priority_queue<indices,std::vector<indices>,
			std::function<bool(indices,indices)> > queue(compare);
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
	
	if( c_[mink] > minV || cp_[mink] == 0 ){ 
	  c_[mink] = minV; open--;
	  if( mink != c_.size() && minV < minimum_ ){
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
