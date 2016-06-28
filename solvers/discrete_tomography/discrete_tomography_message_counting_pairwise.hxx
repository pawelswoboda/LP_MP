#ifndef LP_MP_DT_COUNTING_PAIRWISE_MESSAGE_HXX
#define LP_MP_DT_COUNTING_PAIRWISE_MESSAGE_HXX

#include "LP_MP.h"
#include "minConv.hxx"
#include <math.h> 

namespace LP_MP {
    
  //using MinConv = MinConv<REAL,INDEX>;

  class DiscreteTomographyMessageCountingPairwise{

  public:

    DiscreteTomographyMessageCountingPairwise(INDEX numberOfLabels,INDEX numberOfVarsLeft,INDEX numberOfVarsRight,INDEX SumBound);

    // RIGHT -> TOP Ternary
    // LEFT  -> Pairwise

    /* repam.GetFactor()->&f  */
      
    template<typename RIGHT_FACTOR, typename G1, typename G2>
    void ReceiveMessageFromRight(RIGHT_FACTOR* const f_right, const G1& repam_right, G2& msg);

    template<typename LEFT_FACTOR, typename G1, typename G3>
    void SendMessageToRight(LEFT_FACTOR* const f_left, const G1& repam_left, G3& msg, const REAL omega);

    /*------*/
      
    //template<typename LEFT_FACTOR, typename G1, typename G2>
    //void ReceiveMessageFromLeft(LEFT_FACTOR* const f_left, const G1& repam_left, G2& msg);

    //template<typename LEFT_FACTOR, typename RIGHT_FACTOR, typename G2, typename G3>
    //void SendMessageToLeft(RIGHT_FACTOR* const f_right, const G2& repam_right, G3& msg, const REAL omega);

    /*------*/
    // Update repam with message for each factor ?
      
    template<typename G>
    void RepamLeft(G& repam, const REAL msg, const INDEX msg_dim);
    
    template<typename G>
    void RepamRight(G& repam, const REAL msg, const INDEX msg_dim);

    void ComputeRightFromLeftPrimal(const typename PrimalSolutionStorage::Element left, typename PrimalSolutionStorage::Element right);

    
  private:
    const INDEX numberOfLabels_,numberOfVarsLeft_,numberOfVarsRight_,SumBound_;
    INDEX upSize_,leftSize_,rightSize_,regSize_;
  };

  DiscreteTomographyMessageCountingPairwise::DiscreteTomographyMessageCountingPairwise(INDEX numberOfLabels,INDEX numberOfVarsLeft,INDEX numberOfVarsRight,INDEX SumBound)
    : numberOfLabels_(numberOfLabels),numberOfVarsLeft_(numberOfVarsLeft),numberOfVarsRight_(numberOfVarsRight),SumBound_(SumBound) {
    assert(numberOfLabels_ > 1);
    assert(numberOfVarsLeft_ > 0);
    assert(numberOfVarsRight_ > 0);
    assert(SumBound_ > 0);
    
    upSize_ = pow(numberOfLabels_,2)*std::min(((numberOfVarsLeft_+numberOfVarsRight_)*(numberOfLabels_-1)+1),SumBound);
    leftSize_ = pow(numberOfLabels_,2)*std::min((numberOfVarsLeft_*(numberOfLabels_-1)+1),SumBound);
    rightSize_ = pow(numberOfLabels_,2)*std::min((numberOfVarsRight_*(numberOfLabels_-1)+1),SumBound);
    
    regSize_ = pow(numberOfLabels_,2);
  }

  template<typename LEFT_FACTOR, typename G1, typename G3>
  void DiscreteTomographyMessageCountingPairwise::SendMessageToRight(LEFT_FACTOR* const f_left, const G1& repam_left, G3& msg, const REAL omega){
    //printf("START\n");
    assert(msg.size() == pow(numberOfLabels_,2));
    assert(repam_left.size() == pow(numberOfLabels_,2));

    for(INDEX i=0;i<pow(numberOfLabels_,2);i++){
      //printf("send reg(%d) --> count: %.3e * %.3e\n",i,omega,repam_left[i]);
      msg[i] -= omega*repam_left[i];
    }     
  }
    
  template<typename RIGHT_FACTOR, typename G1, typename G2>
  void DiscreteTomographyMessageCountingPairwise::ReceiveMessageFromRight(RIGHT_FACTOR* const f_right, const G1& repam_right, G2& msg){
    assert(msg.size() == pow(numberOfLabels_,2));
    assert(repam_right.size() == ((*f_right).getSize(DiscreteTomographyFactorCounting::NODE::up) +
				  (*f_right).getSize(DiscreteTomographyFactorCounting::NODE::left) +
				  (*f_right).getSize(DiscreteTomographyFactorCounting::NODE::right) +
				  (*f_right).getSize(DiscreteTomographyFactorCounting::NODE::reg)));

    INDEX left_size = (*f_right).getSize(DiscreteTomographyFactorCounting::NODE::left)/pow(numberOfLabels_,2);
    INDEX right_size =(*f_right).getSize(DiscreteTomographyFactorCounting::NODE::right)/pow(numberOfLabels_,2);
    INDEX up_size = (*f_right).getSize(DiscreteTomographyFactorCounting::NODE::up)/pow(numberOfLabels_,2);
	
    auto op = [&](INDEX i,INDEX j){ return (i+j < up_size) ? i+j : up_size; };
	
    for(INDEX i=0;i<pow(numberOfLabels_,2);i++){
      REAL m = std::numeric_limits<REAL>::infinity();
      INDEX b = i % numberOfLabels_;
      INDEX c = ((i-b)/numberOfLabels_) % numberOfLabels_;
      
      REAL reg = repam_right[up_size*pow(numberOfLabels_,2) + left_size*pow(numberOfLabels_,2) + right_size*pow(numberOfLabels_,2) + b + c*numberOfLabels_];
      
      for(INDEX j=0;j<pow(numberOfLabels_,2);j++){
	INDEX a = j % numberOfLabels_;
	INDEX d = ((j-a)/numberOfLabels_) % numberOfLabels_;

	auto z_up = [&](INDEX k){ return repam_right[a + d*numberOfLabels_ + k*pow(numberOfLabels_,2)];  };
	auto z_left = [&](INDEX k){ return repam_right[up_size*pow(numberOfLabels_,2) + a + b*numberOfLabels_ + k*pow(numberOfLabels_,2)];  };
	auto z_right = [&](INDEX k){ return repam_right[up_size*pow(numberOfLabels_,2) + left_size*pow(numberOfLabels_,2) + c + d*numberOfLabels_ + k*pow(numberOfLabels_,2)];  };
	  	  
	MinConv mc(z_left,z_right,left_size,right_size,up_size);
	mc.CalcConv(op,z_left,z_right);

	for(INDEX k=0;k<up_size;k++){
	  assert( k == (mc.getIdxA(k) + mc.getIdxB(k)));
	  
	  REAL val = mc.getConv(k) + z_up(k);
	  m = std::min(m,val);
	}
      }
      assert( (m + reg) > -1.0e-02);
      msg[i] -= (m + reg);
    }    
  }

  template<typename G>
  void DiscreteTomographyMessageCountingPairwise::RepamLeft(G& repam, const REAL msg, const INDEX msg_dim){

    auto f = repam.GetFactor();
    assert( repam.size() == pow(numberOfLabels_,2));

    assert(msg_dim < pow(numberOfLabels_,2));
    assert(repam[msg_dim] > -std::numeric_limits<REAL>::max());
    if( std::isfinite(msg) ){ repam[msg_dim] += msg; } 
    
    else{ repam[msg_dim] = std::numeric_limits<REAL>::infinity(); }
    assert(repam[msg_dim] > -1.0e-02);
  }

  template<typename G>
  void DiscreteTomographyMessageCountingPairwise::RepamRight(G& repam, const REAL msg, const INDEX msg_dim){

    auto f = repam.GetFactor();
    assert( repam.size() == (f->getSize(DiscreteTomographyFactorCounting::NODE::up) +
			     f->getSize(DiscreteTomographyFactorCounting::NODE::left) +
			     f->getSize(DiscreteTomographyFactorCounting::NODE::right) +
			     f->getSize(DiscreteTomographyFactorCounting::NODE::reg))
	    );
    
    assert(msg_dim < pow(numberOfLabels_,2));
    assert(repam[f->getSize(DiscreteTomographyFactorCounting::NODE::up) +
		 f->getSize(DiscreteTomographyFactorCounting::NODE::left) +
		 f->getSize(DiscreteTomographyFactorCounting::NODE::right) +
		 msg_dim] > -std::numeric_limits<REAL>::max());
    
    if( std::isfinite(msg) ){
	repam[f->getSize(DiscreteTomographyFactorCounting::NODE::up) +
	      f->getSize(DiscreteTomographyFactorCounting::NODE::left) +
	      f->getSize(DiscreteTomographyFactorCounting::NODE::right) +
	      msg_dim] +=  msg;
	assert(repam[f->getSize(DiscreteTomographyFactorCounting::NODE::up) +
		     f->getSize(DiscreteTomographyFactorCounting::NODE::left) +
		     f->getSize(DiscreteTomographyFactorCounting::NODE::right) +
		     msg_dim] > -std::numeric_limits<REAL>::max());
      }
    else{
      repam[f->getSize(DiscreteTomographyFactorCounting::NODE::up) +
	    f->getSize(DiscreteTomographyFactorCounting::NODE::left) +
	    f->getSize(DiscreteTomographyFactorCounting::NODE::right) +
	    msg_dim] = std::numeric_limits<REAL>::infinity();
    }
  }

  void DiscreteTomographyMessageCountingPairwise::ComputeRightFromLeftPrimal(const typename PrimalSolutionStorage::Element left, typename PrimalSolutionStorage::Element right){

    INDEX opt = 0;
    INDEX count = 0;

    for(INDEX i=0;i<pow(numberOfLabels_,2);i++){
      if( left[i] == 1 ){
	opt = i; count++;
	right[upSize_ + leftSize_ + rightSize_ + i] = 1;
      } //else{
	//right[upSize_ + leftSize_ + rightSize_ + i] = 0;
      //}
    }
    assert(count <= 1);
    
    INDEX a = 0;
    INDEX b = opt % numberOfLabels_;
    INDEX c = ((opt - b)/numberOfLabels_) % numberOfLabels_;
    INDEX d = 0;
    
    INDEX lIdx = 0;
    INDEX rIdx = 0;

    if( count == 1 ){
      if(numberOfVarsLeft_ == 1){
	lIdx = b + b*numberOfLabels_ + b*pow(numberOfLabels_,2);
	//for(INDEX i=0;i<leftSize_;i++){
	//right[upSize_ + i]=0;
	//}
	right[upSize_ + lIdx]=1;
	count++;
	a = b;
	lIdx = b;
      }
      else{
	INDEX tcount = 0;
	opt = 0;
	for(INDEX i=0;i<leftSize_;i++){
	  if( right[upSize_ + i] == 1 ){
	    opt = i; tcount++; count++;
	  }
	  assert(tcount <= 1);
	}
	a = opt % numberOfLabels_;
	opt = (opt-a)/numberOfLabels_;
	INDEX tb = opt % numberOfLabels_;
	assert(tb == b);
      
	opt = (opt-tb)/numberOfLabels_;
	lIdx = opt % numberOfLabels_;
      }
    
      if(numberOfVarsRight_ == 1){
	rIdx = c + c*numberOfLabels_ + c*pow(numberOfLabels_,2);
	//for(INDEX i=0;i<rightSize_;i++){
	// right[upSize_ + leftSize_ + i]=0;
	//}
	right[upSize_ + leftSize_ + rIdx]=1;
	count++;
	d = c;
	rIdx = c;
      }
      else{
	INDEX tcount = 0;
	opt = 0;
	for(INDEX i=0;i<rightSize_;i++){
	  if( right[upSize_ + leftSize_ + i] == 1 ){
	    opt = i; tcount++; count++;
	  }
	  assert(tcount <= 1);
	}
	INDEX tc = opt % numberOfLabels_;
	opt = (opt-tc)/numberOfLabels_;
	d = opt % numberOfLabels_;
	assert(tc == c);

	opt = (opt-d)/numberOfLabels_;
	rIdx = opt % numberOfLabels_;
      }
    }
    
    if( count == 3 ){
      //for(INDEX i=0;i<upSize_;i++){
      // right[i] = 0;
      //}
      INDEX z = lIdx + rIdx;
      assert(z < (upSize_/pow(numberOfLabels_,2)));
      INDEX idx = a + d*numberOfLabels_ + z*pow(numberOfLabels_,2);
      assert(idx < upSize_);
      right[idx] = 1;
    }
  }
  
}


#endif // LP_MP_DT_COUNTING_MESSAGE_HXX
