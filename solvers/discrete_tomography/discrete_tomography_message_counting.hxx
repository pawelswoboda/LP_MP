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

    DiscreteTomographyMessageCounting(INDEX numberOfLabels,INDEX numberOfVars);

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

    //void ComputeRightFromLeftPrimal(const bool leftPrimal, MulticutTripletFactor::LabelingType& rightPrimal);

  private:
    const INDEX numberOfLabels_,numberOfVars_;
  };

  template<DIRECTION DR>
  DiscreteTomographyMessageCounting<DR>::DiscreteTomographyMessageCounting(INDEX numberOfLabels,INDEX numberOfVars)
    : numberOfLabels_(numberOfLabels),numberOfVars_(numberOfVars) {
    assert(numberOfLabels_ > 1);
    assert(numberOfVars_ > 0);
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
      
    for(INDEX i=0;i<pow(numberOfLabels_,4);i++){
      INDEX idx = i;
      INDEX a = idx % numberOfLabels_;
      idx = ( idx - a )/numberOfLabels_;
      INDEX b = idx % numberOfLabels_;
      idx = ( idx - b )/numberOfLabels_;
      INDEX c = idx % numberOfLabels_;
      idx = ( idx - c )/numberOfLabels_;
      INDEX d = idx % numberOfLabels_;

      //auto z_up = [&](INDEX k){
      //  return repam_left[a + d*numberOfLabels_ + k*pow(numberOfLabels_,2)];  };
      auto z_left = [&](INDEX k){
	return repam_left[up_size*pow(numberOfLabels_,2) + a + b*numberOfLabels_ + k*pow(numberOfLabels_,2)];  };
      auto z_right = [&](INDEX k){
	return repam_left[up_size*pow(numberOfLabels_,2) + left_size*pow(numberOfLabels_,2) + c + d*numberOfLabels_ + k*pow(numberOfLabels_,2)];  };
      
      REAL reg = repam_left[up_size*pow(numberOfLabels_,2) + left_size*pow(numberOfLabels_,2) + right_size*pow(numberOfLabels_,2) + b + c*numberOfLabels_];
      if( reg <= -std::numeric_limits<REAL>::max() ){ reg = std::numeric_limits<REAL>::max(); }
      
      MinConv mc(z_left,z_right,left_size,right_size,up_size);
      mc.CalcConv(op,z_left,z_right);

      for(INDEX k=0;k<up_size;k++){
	assert(k == (mc.getIdxA(k) + mc.getIdxB(k)));
	assert(!std::isnan(reg));
	REAL val = mc.getConv(k) + reg;
	INDEX kidx = a + numberOfLabels_*d + k*pow(numberOfLabels_,2);
	assert(kidx < (up_size*pow(numberOfLabels_,2)));
	assert(!std::isnan(msg[kidx]));
	assert(!std::isnan(val));
	//printf(" --> msg %.3e ... val %.3e ... omega %.3e \n",msg[kidx],val,omega);
	if( msg[kidx] < std::numeric_limits<REAL>::max() &&
	    msg[kidx] > -std::numeric_limits<REAL>::max()){ msg[kidx] -= omega*val; assert(!std::isnan(msg[kidx])); }
      }
	
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
	auto z_right = [&](INDEX k){
	  return repam_right[up_size*pow(numberOfLabels_,2) + left_size*pow(numberOfLabels_,2) + c + d*numberOfLabels_ + k*pow(numberOfLabels_,2)];  };
	
	REAL reg = repam_right[up_size*pow(numberOfLabels_,2) + left_size*pow(numberOfLabels_,2) + right_size*pow(numberOfLabels_,2) + b + c*numberOfLabels_];
	if( reg <= -std::numeric_limits<REAL>::max() ){ reg = std::numeric_limits<REAL>::max(); }
	
	MinConv mc(z_up,z_right,up_size,right_size,left_size);
	mc.CalcConv(op,z_up,z_right);

	for(INDEX k=0;k<left_size;k++){
	  assert(k == (mc.getIdxA(k) - mc.getIdxB(k)));
	  //printf(" --> c[k]=%.3e ... reg %.3e\n",mc.getConv(k),reg);
	  assert(!std::isnan(reg));
	  REAL val = mc.getConv(k) + reg;
	  INDEX kidx = a + numberOfLabels_*b + k*pow(numberOfLabels_,2);
	  assert(kidx < (left_size*pow(numberOfLabels_,2)));
	  assert(!std::isnan(msg[kidx]));
	  assert(!std::isnan(val));
	  if( msg[kidx] < std::numeric_limits<REAL>::max() &&
	      msg[kidx] > -std::numeric_limits<REAL>::max()){ msg[kidx] -= val; assert(!std::isnan(msg[kidx]));}
	}
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
	  return repam_right[up_size*pow(numberOfLabels_,2) + c + d*numberOfLabels_ + k*pow(numberOfLabels_,2)];  };
	
	REAL reg = repam_right[up_size*pow(numberOfLabels_,2) + left_size*pow(numberOfLabels_,2) + right_size*pow(numberOfLabels_,2) + b + c*numberOfLabels_];
	if( reg <= -std::numeric_limits<REAL>::max() ){ reg = std::numeric_limits<REAL>::max(); }
	
	MinConv mc(z_up,z_left,up_size,left_size,right_size);
	mc.CalcConv(op,z_up,z_left);

	for(INDEX k=0;k<right_size;k++){
	  assert(k == op(mc.getIdxA(k),mc.getIdxB(k)));//(mc.getIdxA(k) - mc.getIdxB(k)));
	  assert(!std::isnan(reg));
	  REAL val = mc.getConv(k) + reg;
	  INDEX kidx = c + numberOfLabels_*d + k*pow(numberOfLabels_,2);
	  assert(kidx < (right_size*pow(numberOfLabels_,2)));
	  assert(!std::isnan(msg[kidx]));
	  assert(!std::isnan(val));
	  if( msg[kidx] < std::numeric_limits<REAL>::max() &&
	      msg[kidx] > -std::numeric_limits<REAL>::max() ){ msg[kidx] -= val; assert(!std::isnan(msg[kidx]));}
	}
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
    if( msg < std::numeric_limits<REAL>::max() &&
	msg > -std::numeric_limits<REAL>::max()){ repam[msg_dim] += msg; }
    else{ repam[msg_dim] = std::numeric_limits<REAL>::max(); }
      
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
      if( msg < std::numeric_limits<REAL>::max() &&
	  msg > -std::numeric_limits<REAL>::max() ){
	repam[f->getSize(DiscreteTomographyFactorCounting::NODE::up) + msg_dim] +=  msg;
      }
      else{
	repam[f->getSize(DiscreteTomographyFactorCounting::NODE::up) + msg_dim] = std::numeric_limits<REAL>::max();
      }
    }
    else{
      assert(msg_dim < f->getSize(DiscreteTomographyFactorCounting::NODE::right));
      assert(repam[f->getSize(DiscreteTomographyFactorCounting::NODE::up) + f->getSize(DiscreteTomographyFactorCounting::NODE::left) + msg_dim] > -std::numeric_limits<REAL>::max());
      if( msg < std::numeric_limits<REAL>::max() &&
	  msg > -std::numeric_limits<REAL>::max() ){
	repam[f->getSize(DiscreteTomographyFactorCounting::NODE::up) + f->getSize(DiscreteTomographyFactorCounting::NODE::left) + msg_dim] += msg;
      }
      else{
	repam[f->getSize(DiscreteTomographyFactorCounting::NODE::up) + f->getSize(DiscreteTomographyFactorCounting::NODE::left) + msg_dim] = std::numeric_limits<REAL>::max();
      }
    }
  }
    
  
}

#endif // LP_MP_DT_COUNTING_MESSAGE_HXX
