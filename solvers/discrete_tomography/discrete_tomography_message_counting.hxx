#ifndef LP_MP_DT_COUNTING_MESSAGE_HXX
#define LP_MP_DT_COUNTING_MESSAGE_HXX

#include "LP_MP.h"
#include "minConv.hxx"
#include <math.h> 

namespace LP_MP {
  namespace discrete_tomo{
    
    using MinConv = MinConv<std::function<REAL(INDEX)>,REAL,INDEX>;

    enum class DIRECTION {left,right};

    template<DIRECTION D>
    class DiscreteTomographyMessageCounting{

    public:

      DiscreteTomographyMessageCounting(INDEX numberOfLabels,INDEX,DIRECTION);

      // RIGHT -> TOP Ternary
      // LEFT  -> BOTTOM Ternary

      /* repam.GetFactor()->&f  */
      
      template<typename RIGHT_FACTOR, typename G1, typename G2>
      void ReceiveMessageFromRight(RIGHT_FACTOR* const f_right, const G1& repam_right, G2& msg);

      template<typename LEFT_FACTOR, typename RIGHT_FACTOR, typename G1, typename G3>
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

      //void ComputeRightFromLeftPrimal(const bool leftPrimal, MulticutTripletFactor::LabelingType& rightPrimal);
    
    private:
      const INDEX numberOfLabels_;
      const DIRECTION dr_;
      INDEX msgSize_;
      INDEX leftSize_,rightSize_,upSize_,regSize_;
    
    };

    DiscreteTomographyMessageCounting::DiscreteTomographyMessageCounting(INDEX numberOfLabels,std::vector<INDEX> top,std::vector<INDEX> bottom,DIRECTION dr);
      : numberOfLabels_(numberOfLabels),dr_(dr) {
	assert(numberOfLabels_ > 1);

	
	
	upSize_ = pow(numberOfLabels_,2)*((d_-a_+1)(numberOfLabels_-1)+numberOfLabels_);
	leftSize_ = pow(numberOfLabels_,2)*((b_-a_+1)(numberOfLabels_-1)+numberOfLabels_);
	rightSize_ = pow(numberOfLabels_,2)*((c_-b_+1)(numberOfLabels_-1)+numberOfLabels_);
	regSize_ = pow(numberOfLabels_,2);
	
	if( dr_ == DIRECTION::up ){
	  msgSize_ = upSize_;
	}
	else if( dr_ == DIRECTION::left ){
	  msgSize_ = leftSize_;
	}
	else if( dr_ == DIRECTION::right ){
	  msgSize_ = rightSize_;
	}	
	
      }

    template<typename RIGHT_FACTOR, typename G1, typename G2>
    void DiscreteTomographyMessageCounting::ReceiveMessageFromRight(RIGHT_FACTOR* const f_right, const G1& repam_right, G2& msg){

      op = [&](INDEX i,INDEX j){ return i-j < 0 ? f_right.getSize(f_right::left) : i-j;  };
      if( dr_ == DIRECTION::left ){
	assert(msg.size() == f_right.getSize(f_right::left));
	assert(repam_right.size() == (f_right.getSize(f_right::up) + f_right.getSize(f_right::left) +
				      f_right.getSize(f_right::right) + f_right.getSize(f_right::reg)));

	op = [&](INDEX i,INDEX j){ return i-j < 0 ? f_right.getSize(f_right::left) : i-j;  };
	INDEX z_size = f_right.getSize(f_right::left)/pow(numberOfLabels_,2);
	INDEX up_size = f_right.getSize(f_right::up)/pow(numberOfLabels_,2);
	INDEX right_size = f_right.getSize(f_right::right)/pow(numberOfLabels_,2);
	INDEX left_size = f_right.getSize(f_right::left)/pow(numberOfLabels_,2);
	
	for(INDEX i=0;i<pow(numberOfLabels_,4);i++){
	  INDEX idx = i;
	  INDEX a = idx % numberOfLabels_;
	  idx = ( idx - a )/numberOfLabels_;
	  INDEX b = idx % numberOfLabels_;
	  idx = ( idx - b )/numberOfLabels_;
	  INDEX c = idx % numberOfLabels_;
	  idx = ( idx - c )/numberOfLabels_;
	  INDEX d = idx % numberOfLabels_;

	  auto z_up = [&](INDEX k){ return repam_right[c + d*numberOfLabels_ + k*pow(numberOfLabels_,2)];  };
	  auto z_left = [&](INDEX k){ return repam_right[f_right.getSize(f_right::up) + a + b*numberOfLabels_ + k*pow(numberOfLabels_,2)];  };
	  auto z_right = [&](INDEX k){ return repam_right[f_right.getSize(f_right::up) + f_right.getSize(f_right::left) + c + d*numberOfLabels_ + k*pow(numberOfLabels_,2)];  };
	  auto z_reg = repam_right[f_right.getSize(f_right::up) + f_right.getSize(f_right::left) + f_right.getSize(f_right::right) + b + c*numberOfLabels_];

	  
	  
	  MinConv mc(z_up,z_right,up_size,right_size,left_size);
	  mc.CalcConv(op);

	  for(INDEX k=0;k<z_size;k++){
	    REAL val = mc.getConv(k) + f_right.eval(mc.getIdxA(k),k,mc.getIdxB(k));
	    assert(val < std::numeric_limits<REAL>::max() );
	    INDEX kidx = a + numberOfLabels_*b + k*pow(numberOfLabels_,2);
	    msg[kidx] -= val;
	  }
	  
	}
	
      }
      
      
      
      
    }

    template<typename G>
    void DiscreteTomographyMessageCounting::RepamLeft(G& repam, const REAL msg, const INDEX msg_dim){

      auto f = repam.GetFactor();
      assert( repam.size() == (f.getSize(f::left) + f.getSize(f::right) + f.getSize(f::up) + f.getSize(f::reg))); // <-- wie kann man das machen?

      assert(msg_dim < f.getSize(f::up));
      repam[msg_dim] += msg;
      
    }

    template<typename G>
    void DiscreteTomographyMessageCounting::RepamRight(G& repam, const REAL msg, const INDEX msg_dim){

      auto f = repam.GetFactor();
      assert( repam.size() == (f.getSize(f::left) + f.getSize(f::right) + f.getSize(f::up) + f.getSize(f::reg))); // <-- wie kann man das machen?

      if( dr_ == DIRECTION::left ){
	assert(msg_dim < f.getSize(f::left));
	repam[f.getSize(f::up) + msg_dim] +=  msg;
      }
      else{
	assert(msg_dim < f.getSize(f::right));
	repam[f.getSize(f::up) + f.getSize(f::left) + msg_dim] += msg;
      }
    }
    
  }
}

#endif // LP_MP_DT_COUNTING_MESSAGE_HXX
