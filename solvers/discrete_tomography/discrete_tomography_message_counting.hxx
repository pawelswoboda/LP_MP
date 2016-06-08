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
  DiscreteTomographyMessageCounting<DR>::DiscreteTomographyMessageCounting(INDEX numberOfLabels,INDEX numberOfVars);
  : numberOfLabels_(numberOfLabels),numberOfVars_(numberOfVars) {
    assert(numberOfLabels_ > 1);
    assert(numberOfVars_ > 0);
  }

  template<DIRECTION DR>
  template<typename LEFT_FACTOR, typename G1, typename G3>
  void DiscreteTomographyMessageCounting::SendMessageToRight(LEFT_FACTOR* const f_left, const G1& repam_left, G3& msg, const REAL omega){
    assert(msg.size() == f_left.getSize(f_left::up));
    assert(repam_right.size() == (f_left.getSize(f_left::up) + f_left.getSize(f_left::left) +
				  f_left.getSize(f_left::right) + f_left.getSize(f_left::reg)));

     
    INDEX up_size = f_left.getSize(f_left::up)/pow(numberOfLabels_,2);
    INDEX right_size = f_left.getSize(f_left::right)/pow(numberOfLabels_,2);
    INDEX left_size = f_left.getSize(f_left::left)/pow(numberOfLabels_,2);

    op = [&](INDEX i,INDEX j){ return (i+j < up_size) ? i+j : up_size;  }; // 0 <= i+j < up_size
      
    for(INDEX i=0;i<pow(numberOfLabels_,4);i++){
      INDEX idx = i;
      INDEX a = idx % numberOfLabels_;
      idx = ( idx - a )/numberOfLabels_;
      INDEX b = idx % numberOfLabels_;
      idx = ( idx - b )/numberOfLabels_;
      INDEX c = idx % numberOfLabels_;
      idx = ( idx - c )/numberOfLabels_;
      INDEX d = idx % numberOfLabels_;

      //auto z_up = [&](INDEX k){ return repam_right[a + d*numberOfLabels_ + k*pow(numberOfLabels_,2)];  };
      auto z_left = [&](INDEX k){ return repam_right[f_right.getSize(f_right::up) + a + b*numberOfLabels_ + k*pow(numberOfLabels_,2)];  };
      auto z_right = [&](INDEX k){ return repam_right[f_right.getSize(f_right::up) + f_right.getSize(f_right::left) + c + d*numberOfLabels_ + k*pow(numberOfLabels_,2)];  };
      REAL reg = repam_right[f_right.getSize(f_right::up) + f_right.getSize(f_right::left) + f_right.getSize(f_right::right) + b + c*numberOfLabels_];
		  
      MinConv mc(z_left,z_right,left_size,right_size,up_size);
      mc.CalcConv(op);

      for(INDEX k=0;k<up_size;k++){
	assert(k == (mc.getIdxA(k) + mc.getIdxB(k));
	       REAL val = mc.getConv(k) + reg;
	       INDEX kidx = a + numberOfLabels_*d + k*pow(numberOfLabels_,2);
	       assert(kidx < (up_size*pow(numberOfLabels_,2)));
	       msg[kidx] -= omega*val;
	       }
	
      }
    }
    
    template<DIRECTION DR>
      template<typename RIGHT_FACTOR, typename G1, typename G2>
      void DiscreteTomographyMessageCounting::ReceiveMessageFromRight(RIGHT_FACTOR* const f_right, const G1& repam_right, G2& msg){

      if( DR == DIRECTION::left ){
	assert(msg.size() == f_right.getSize(f_right::left));
	assert(repam_right.size() == (f_right.getSize(f_right::up) + f_right.getSize(f_right::left) +
				      f_right.getSize(f_right::right) + f_right.getSize(f_right::reg)));

	INDEX left_size = f_right.getSize(f_right::left)/pow(numberOfLabels_,2);
	INDEX right_size = f_right.getSize(f_right::right)/pow(numberOfLabels_,2);
	INDEX up_size = f_right.getSize(f_right::up)/pow(numberOfLabels_,2);
	
	op = [&](INDEX i,INDEX j){ // 0 <= i-j < left_size
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

	  auto z_up = [&](INDEX k){ return repam_right[a + d*numberOfLabels_ + k*pow(numberOfLabels_,2)];  };
	  //auto z_left = [&](INDEX k){ return repam_right[f_right.getSize(f_right::up) + a + b*numberOfLabels_ + k*pow(numberOfLabels_,2)];  };
	  auto z_right = [&](INDEX k){ return repam_right[f_right.getSize(f_right::up) + f_right.getSize(f_right::left) + c + d*numberOfLabels_ + k*pow(numberOfLabels_,2)];  };
	  REAL reg = repam_right[f_right.getSize(f_right::up) + f_right.getSize(f_right::left) + f_right.getSize(f_right::right) + b + c*numberOfLabels_];
	  
	  MinConv mc(z_up,z_right,up_size,right_size,left_size);
	  mc.CalcConv(op);

	  for(INDEX k=0;k<left_size;k++){
	    assert(k == (mc.getIdxA(k) - mc.getIdxB(k));
		   REAL val = mc.getConv(k) + reg;
		   INDEX kidx = a + numberOfLabels_*b + k*pow(numberOfLabels_,2);
		   assert(kidx < (left_size*pow(numberOfLabels_,2)));
		   msg[kidx] -= val;
		   }
	  }
	}
	else{
	  assert(msg.size() == f_right.getSize(f_right::right));
	  assert(repam_right.size() == (f_right.getSize(f_right::up) + f_right.getSize(f_right::left) +
					f_right.getSize(f_right::right) + f_right.getSize(f_right::reg)));

	  INDEX left_size = f_right.getSize(f_right::left)/pow(numberOfLabels_,2);
	  INDEX right_size = f_right.getSize(f_right::right)/pow(numberOfLabels_,2);
	  INDEX up_size = f_right.getSize(f_right::up)/pow(numberOfLabels_,2);

	  op = [&](INDEX i,INDEX j){ // 0 <= i-j < right_size
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

	    auto z_up = [&](INDEX k){ return repam_right[a + d*numberOfLabels_ + k*pow(numberOfLabels_,2)];  };
	    auto z_left = [&](INDEX k){ return repam_right[f_right.getSize(f_right::up) + a + b*numberOfLabels_ + k*pow(numberOfLabels_,2)];  };
	    //auto z_right = [&](INDEX k){ return repam_right[f_right.getSize(f_right::up) + f_right.getSize(f_right::left) + c + d*numberOfLabels_ + k*pow(numberOfLabels_,2)];  };
	    REAL reg = repam_right[f_right.getSize(f_right::up) + f_right.getSize(f_right::left) + f_right.getSize(f_right::right) + b + c*numberOfLabels_];

	    MinConv mc(z_up,z_left,up_size,left_size,right_size);
	    mc.CalcConv(op);

	    for(INDEX k=0;k<right_size;k++){
	      assert(f_right.eval(mc.getIdxA(k),mc.getIdxB(k),k) < std::numeric_limits<REAL>::max() );
	      REAL val = mc.getConv(k) + reg;
	      INDEX kidx = c + numberOfLabels_*d + k*pow(numberOfLabels_,2);
	      assert(kidx < (right_size*pow(numberOfLabels_,2)));
	      msg[kidx] -= val; 
	    }
	  }
	}
            
      }

      template<DIRECTION DR>
	template<typename G>
	void DiscreteTomographyMessageCounting::RepamLeft(G& repam, const REAL msg, const INDEX msg_dim){

	auto f = repam.GetFactor();
	assert( repam.size() == (f.getSize(f::left) + f.getSize(f::right) + f.getSize(f::up) + f.getSize(f::reg)));

	assert(msg_dim < f.getSize(f::up));
	repam[msg_dim] += msg;
      
      }

      template<DIRECTION DR>    
	template<typename G>
	void DiscreteTomographyMessageCounting::RepamRight(G& repam, const REAL msg, const INDEX msg_dim){

	auto f = repam.GetFactor();
	assert( repam.size() == (f.getSize(f::left) + f.getSize(f::right) + f.getSize(f::up) + f.getSize(f::reg)));

	if( DR == DIRECTION::left ){
	  assert(msg_dim < f.getSize(f::left));
	  repam[f.getSize(f::up) + msg_dim] +=  msg;
	}
	else{
	  assert(msg_dim < f.getSize(f::right));
	  repam[f.getSize(f::up) + f.getSize(f::left) + msg_dim] += msg;
	}
      }
    
  
    }

#endif // LP_MP_DT_COUNTING_MESSAGE_HXX
