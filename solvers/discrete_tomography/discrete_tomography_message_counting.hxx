#ifndef LP_MP_DT_COUNTING_MESSAGE_HXX
#define LP_MP_DT_COUNTING_MESSAGE_HXX

#include "LP_MP.h"
#include "minConv.hxx"
#include <math.h> 

namespace LP_MP {
    
  //using MinConv = MinConv<std::function<REAL(INDEX)>,REAL,INDEX>;

  enum class DIRECTION {left,right};

  template<DIRECTION DR>
  class DiscreteTomographyMessageCounting{

  public:

    DiscreteTomographyMessageCounting(INDEX numberOfLabels,INDEX numberOfVarsLeft,INDEX numberOfVarsRight,INDEX SumBound);

    // RIGHT -> TOP Ternary
    // LEFT  -> BOTTOM Ternary

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
      
    template<typename G>
    void RepamLeft(G& repam, const REAL msg, const INDEX msg_dim);
    
    template<typename G>
    void RepamRight(G& repam, const REAL msg, const INDEX msg_dim);

    /*------*/

    void ComputeRightFromLeftPrimal(const typename PrimalSolutionStorage::Element left, typename PrimalSolutionStorage::Element right);
    
  private:
    const INDEX numberOfLabels_,numberOfVarsLeft_,numberOfVarsRight_,SumBound_;
    INDEX upSize_,leftSize_,rightSize_,regSize_;
  };

  template<DIRECTION DR>
  DiscreteTomographyMessageCounting<DR>::DiscreteTomographyMessageCounting(INDEX numberOfLabels,INDEX numberOfVarsLeft,INDEX numberOfVarsRight,INDEX SumBound)
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

  template<DIRECTION DR>
  template<typename LEFT_FACTOR, typename G1, typename G3>
  void DiscreteTomographyMessageCounting<DR>::SendMessageToRight(LEFT_FACTOR* const f_left, const G1& repam_left, G3& msg, const REAL omega){    
    assert(msg.size() == (*f_left).getSize(DiscreteTomographyFactorCounting::NODE::up));
    assert(repam_left.size() == ((*f_left).getSize(DiscreteTomographyFactorCounting::NODE::up) +
				 (*f_left).getSize(DiscreteTomographyFactorCounting::NODE::left) +
				 (*f_left).getSize(DiscreteTomographyFactorCounting::NODE::right) +
				 (*f_left).getSize(DiscreteTomographyFactorCounting::NODE::reg)));

    INDEX up_size = (*f_left).getSize(DiscreteTomographyFactorCounting::NODE::up)/pow(numberOfLabels_,2);
    INDEX right_size = (*f_left).getSize(DiscreteTomographyFactorCounting::NODE::right)/pow(numberOfLabels_,2);
    INDEX left_size = (*f_left).getSize(DiscreteTomographyFactorCounting::NODE::left)/pow(numberOfLabels_,2);

    auto op = [&](INDEX i,INDEX j){ return (i+j < up_size) ? i+j : up_size;  }; // 0 <= i+j < up_size

    std::vector<REAL> msg_v(up_size*pow(numberOfLabels_,2),std::numeric_limits<REAL>::infinity());
    for(INDEX i=0;i<pow(numberOfLabels_,4);i++){
      INDEX idx = i;
      INDEX a = idx % numberOfLabels_;
      idx = ( idx - a )/numberOfLabels_;
      INDEX b = idx % numberOfLabels_;
      idx = ( idx - b )/numberOfLabels_;
      INDEX c = idx % numberOfLabels_;
      idx = ( idx - c )/numberOfLabels_;
      INDEX d = idx % numberOfLabels_;

      auto z_up = [&](INDEX k){
        return repam_left[a + d*numberOfLabels_ + k*pow(numberOfLabels_,2)];  };
      auto z_left = [&](INDEX k){
	return repam_left[up_size*pow(numberOfLabels_,2) + a + b*numberOfLabels_ + k*pow(numberOfLabels_,2)];  };
      auto z_right = [&](INDEX k){
	return repam_left[up_size*pow(numberOfLabels_,2) + left_size*pow(numberOfLabels_,2) + c + d*numberOfLabels_ + k*pow(numberOfLabels_,2)];  };
      
      REAL reg = repam_left[up_size*pow(numberOfLabels_,2) + left_size*pow(numberOfLabels_,2) + right_size*pow(numberOfLabels_,2) + b + c*numberOfLabels_];
      assert(reg > -std::numeric_limits<REAL>::max());
      
      MinConv mc(z_left,z_right,left_size,right_size,up_size);
      mc.CalcConv(op,z_left,z_right);

      for(INDEX k=0;k<up_size;k++){
	assert(k == (mc.getIdxA(k) + mc.getIdxB(k)));
	assert(!std::isnan(reg));
	
	REAL val = mc.getConv(k) + z_up(k)  + reg;
	INDEX kidx = a + numberOfLabels_*d + k*pow(numberOfLabels_,2);
	
	assert(kidx < (up_size*pow(numberOfLabels_,2)));
	assert(!std::isnan(val));
	msg_v[kidx] = std::min(msg_v[kidx],val);
      }
    }
    for(INDEX i=0;i<msg_v.size();i++){
      assert(msg_v[i] > -1.0e-02);
      msg[i] -= omega*msg_v[i];
    }
  }
    
  template<DIRECTION DR>
  template<typename RIGHT_FACTOR, typename G1, typename G2>
  void DiscreteTomographyMessageCounting<DR>::ReceiveMessageFromRight(RIGHT_FACTOR* const f_right, const G1& repam_right, G2& msg){
    
    assert(repam_right.size() == ((*f_right).getSize(DiscreteTomographyFactorCounting::NODE::up) +
				  (*f_right).getSize(DiscreteTomographyFactorCounting::NODE::left) +
				  (*f_right).getSize(DiscreteTomographyFactorCounting::NODE::right) +
				  (*f_right).getSize(DiscreteTomographyFactorCounting::NODE::reg)));

    INDEX left_size = (*f_right).getSize(DiscreteTomographyFactorCounting::NODE::left)/pow(numberOfLabels_,2);
    INDEX right_size = (*f_right).getSize(DiscreteTomographyFactorCounting::NODE::right)/pow(numberOfLabels_,2);
    INDEX up_size = (*f_right).getSize(DiscreteTomographyFactorCounting::NODE::up)/pow(numberOfLabels_,2);
    
    if( DR == DIRECTION::left ){
      assert(msg.size() == (*f_right).getSize(DiscreteTomographyFactorCounting::NODE::left));
	
      auto op = [&](INDEX i,INDEX j){ // 0 <= i-j < left_size
	if( i < j ){ // i-j < 0
	  return left_size;
	}
	else{
	  return (i-j < left_size) ? i-j : left_size;
	}
      };

      std::vector<REAL> msg_v(left_size*pow(numberOfLabels_,2),std::numeric_limits<REAL>::infinity());
      for(INDEX i=0;i<pow(numberOfLabels_,4);i++){
	INDEX idx = i;
	INDEX a = idx % numberOfLabels_;
	idx = ( idx - a )/numberOfLabels_;
	INDEX b = idx % numberOfLabels_;
	idx = ( idx - b )/numberOfLabels_;
	INDEX c = idx % numberOfLabels_;
	idx = ( idx - c )/numberOfLabels_;
	INDEX d = idx % numberOfLabels_;

	auto z_up = [&](INDEX k){
	  return repam_right[a + d*numberOfLabels_ + k*pow(numberOfLabels_,2)];  };
	auto z_left = [&](INDEX k){
	  return repam_right[up_size*pow(numberOfLabels_,2) + a + b*numberOfLabels_ + k*pow(numberOfLabels_,2)];  };
	auto z_right = [&](INDEX k){
	  return repam_right[up_size*pow(numberOfLabels_,2) + left_size*pow(numberOfLabels_,2) + c + d*numberOfLabels_ + k*pow(numberOfLabels_,2)];  };
	
	REAL reg = repam_right[up_size*pow(numberOfLabels_,2) + left_size*pow(numberOfLabels_,2) + right_size*pow(numberOfLabels_,2) + b + c*numberOfLabels_];
	assert(reg > -std::numeric_limits<REAL>::max());
	
	MinConv mc(z_up,z_right,up_size,right_size,left_size);
	mc.CalcConv(op,z_up,z_right);

	for(INDEX k=0;k<left_size;k++){
	  assert(k == (mc.getIdxA(k) - mc.getIdxB(k)));
	  assert(!std::isnan(reg));

	  REAL val = mc.getConv(k) + z_left(k) + reg;
	  INDEX kidx = a + numberOfLabels_*b + k*pow(numberOfLabels_,2);

	  assert(kidx < (left_size*pow(numberOfLabels_,2)));
	  assert(!std::isnan(val));

	  msg_v[kidx] = std::min(msg_v[kidx],val);
	}
      }
      for(INDEX i=0;i<msg_v.size();i++){
	assert(msg_v[i] > -1.0e-02);
	msg[i] -= msg_v[i];
      }
    }
    else{
      assert(msg.size() == (*f_right).getSize(DiscreteTomographyFactorCounting::NODE::right));

      auto op = [&](INDEX i,INDEX j){ // 0 <= i-j <= right_size
	if( i < j ){ 
	  return right_size;
	}
	else{
	  return (i-j < right_size) ? i-j : right_size;
	}
      };
	
      std::vector<REAL> msg_v(right_size*pow(numberOfLabels_,2),std::numeric_limits<REAL>::infinity());
      for(INDEX i=0;i<pow(numberOfLabels_,4);i++){
	INDEX idx = i;
	INDEX a = idx % numberOfLabels_;
	idx = ( idx - a )/numberOfLabels_;
	INDEX b = idx % numberOfLabels_;
	idx = ( idx - b )/numberOfLabels_;
	INDEX c = idx % numberOfLabels_;
	idx = ( idx - c )/numberOfLabels_;
	INDEX d = idx % numberOfLabels_;

	auto z_up = [&](INDEX k){
	  return repam_right[a + d*numberOfLabels_ + k*pow(numberOfLabels_,2)];  };
	auto z_left = [&](INDEX k){
	  return repam_right[up_size*pow(numberOfLabels_,2) + a + b*numberOfLabels_ + k*pow(numberOfLabels_,2)];  };
	auto z_right = [&](INDEX k){
	  return repam_right[up_size*pow(numberOfLabels_,2) + left_size*pow(numberOfLabels_,2) + c + d*numberOfLabels_ + k*pow(numberOfLabels_,2)];  };
	
	REAL reg = repam_right[up_size*pow(numberOfLabels_,2) + left_size*pow(numberOfLabels_,2) + right_size*pow(numberOfLabels_,2) + b + c*numberOfLabels_];
	assert(reg > -std::numeric_limits<REAL>::max());
	
	MinConv mc(z_up,z_left,up_size,left_size,right_size);
	mc.CalcConv(op,z_up,z_left);

	for(INDEX k=0;k<right_size;k++){
	  assert(k == op(mc.getIdxA(k),mc.getIdxB(k)));//(mc.getIdxA(k) - mc.getIdxB(k)));
	  assert(!std::isnan(reg));

	  REAL val = mc.getConv(k) + z_right(k) + reg;
	  INDEX kidx = c + numberOfLabels_*d + k*pow(numberOfLabels_,2);

	  assert(kidx < (right_size*pow(numberOfLabels_,2)));
	  assert(!std::isnan(val));
	  
	  msg_v[kidx] = std::min(msg_v[kidx],val);
	}
      }
      for(INDEX i=0;i<msg_v.size();i++){
	assert(msg_v[i] > -1.0e-02);
	msg[i] -= msg_v[i];
      }
    }      
  }

  template<DIRECTION DR>
  template<typename G>
  void DiscreteTomographyMessageCounting<DR>::RepamLeft(G& repam, const REAL msg, const INDEX msg_dim){
    
    assert(!std::isnan(msg));
    auto f = repam.GetFactor();
    assert( repam.size() == (f->getSize(DiscreteTomographyFactorCounting::NODE::left) +
			     f->getSize(DiscreteTomographyFactorCounting::NODE::right) +
			     f->getSize(DiscreteTomographyFactorCounting::NODE::up) +
			     f->getSize(DiscreteTomographyFactorCounting::NODE::reg)));

    assert(msg_dim < f->getSize(DiscreteTomographyFactorCounting::NODE::up));
    assert(repam[msg_dim] > -std::numeric_limits<REAL>::max());
    if( std::isfinite(msg) ){ 
      repam[msg_dim] += msg; 
      assert(repam[msg_dim] > -std::numeric_limits<REAL>::max() );
    }
    else{ repam[msg_dim] = std::numeric_limits<REAL>::infinity(); }
      
  }

  template<DIRECTION DR>    
  template<typename G>
  void DiscreteTomographyMessageCounting<DR>::RepamRight(G& repam, const REAL msg, const INDEX msg_dim){
    
    assert(!std::isnan(msg));
    auto f = repam.GetFactor();
    assert( repam.size() == (f->getSize(DiscreteTomographyFactorCounting::NODE::left) +
			     f->getSize(DiscreteTomographyFactorCounting::NODE::right) +
			     f->getSize(DiscreteTomographyFactorCounting::NODE::up) +
			     f->getSize(DiscreteTomographyFactorCounting::NODE::reg)));

    if( DR == DIRECTION::left ){
      assert(msg_dim < f->getSize(DiscreteTomographyFactorCounting::NODE::left));
      assert(repam[f->getSize(DiscreteTomographyFactorCounting::NODE::up) + msg_dim] > -std::numeric_limits<REAL>::max());
      if( std::isfinite(msg) ){
	repam[f->getSize(DiscreteTomographyFactorCounting::NODE::up) + msg_dim] +=  msg; 
	assert(repam[f->getSize(DiscreteTomographyFactorCounting::NODE::up) + msg_dim] > -std::numeric_limits<REAL>::max());
      }
      else{
	repam[f->getSize(DiscreteTomographyFactorCounting::NODE::up) + msg_dim] = std::numeric_limits<REAL>::infinity();
      }
    }
    else{
      assert(msg_dim < f->getSize(DiscreteTomographyFactorCounting::NODE::right));
      assert(repam[f->getSize(DiscreteTomographyFactorCounting::NODE::up) + f->getSize(DiscreteTomographyFactorCounting::NODE::left) + msg_dim] > -std::numeric_limits<REAL>::max());
      if( std::isfinite(msg) ){      
	repam[f->getSize(DiscreteTomographyFactorCounting::NODE::up) + f->getSize(DiscreteTomographyFactorCounting::NODE::left) + msg_dim] += msg; 
	assert(repam[f->getSize(DiscreteTomographyFactorCounting::NODE::up) + f->getSize(DiscreteTomographyFactorCounting::NODE::left) + msg_dim] > -std::numeric_limits<REAL>::max());
      }
      else{
	repam[f->getSize(DiscreteTomographyFactorCounting::NODE::up) + f->getSize(DiscreteTomographyFactorCounting::NODE::left) + msg_dim] = std::numeric_limits<REAL>::infinity();
      }
    }
  }

      
  template<DIRECTION DR>  
  void DiscreteTomographyMessageCounting<DR>::ComputeRightFromLeftPrimal(const typename PrimalSolutionStorage::Element left, typename PrimalSolutionStorage::Element right){
    INDEX count = 0;

    INDEX opt = 0;
    INDEX a = 0;
    INDEX b = 0;
    INDEX c = 0;
    INDEX d = 0;
    INDEX kl = 0;
    INDEX kr = 0;
    
    if( DR == DIRECTION::left ){
      for(INDEX i=0;i<leftSize_;i++){
	if( left[i] == true ){
	  opt = i;
	  count++;
	}
	right[upSize_ + i] = left[i];
      }
            
      a = opt % numberOfLabels_;
      opt = (opt - a)/numberOfLabels_;
      b = opt % numberOfLabels_;
      opt = (opt - b)/numberOfLabels_;
      kl = opt % numberOfLabels_;

      for(INDEX i=0;i<rightSize_;i++){
	if(right[upSize_ + leftSize_ + i] == true){
	  opt = i; count++;
	}
      }
      
      c = opt % numberOfLabels_;
      opt = (opt - c)/numberOfLabels_;
      d = opt % numberOfLabels_;
      opt = (opt - d)/numberOfLabels_;
      kr = opt % numberOfLabels_;
    }
    else{
      for(INDEX i=0;i<rightSize_;i++){
	if( left[i] == true ){
	  opt = i;
	  count++;
	}
	right[upSize_ + rightSize_ + i] = left[i];
      }
      
      c = opt % numberOfLabels_;
      opt = (opt - c)/numberOfLabels_;
      d = opt % numberOfLabels_;
      opt = (opt - d)/numberOfLabels_;
      kr = opt % numberOfLabels_;

      for(INDEX i=0;i<leftSize_;i++){
	if(right[upSize_ + i] == true){
	  opt = i; count++;
	}
      }

      a = opt % numberOfLabels_;
      opt = (opt - a)/numberOfLabels_;
      b = opt % numberOfLabels_;
      opt = (opt - d)/numberOfLabels_;
      kl = opt % numberOfLabels_;					  
    }

    INDEX z = kl+kr;
    if(count == 2 && z<(upSize_/pow(numberOfLabels_,2)) ){
      //assert(z<(upSize_/pow(numberOfLabels_,2)));
      INDEX idx = a + d*numberOfLabels_ + z*pow(numberOfLabels_,2);
      assert(idx < upSize_);
      for(INDEX i=0;i<upSize_;i++){ right[i]=false; }
      right[idx]=true;
    }
    
  }
    
}

#endif // LP_MP_DT_COUNTING_MESSAGE_HXX
