#ifndef LP_MP_DT_COUNTING_PAIRWISE_MESSAGE_HXX
#define LP_MP_DT_COUNTING_PAIRWISE_MESSAGE_HXX

#include "LP_MP.h"
#include "minConv.hxx"
#include <math.h> 

namespace LP_MP {
    
  //using MinConv = MinConv<REAL,INDEX>;

  class DiscreteTomographyMessageCountingPairwise{

  public:

    DiscreteTomographyMessageCountingPairwise(INDEX numberOfLabels);

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

    //void ComputeRightFromLeftPrimal(const bool leftPrimal, MulticutTripletFactor::LabelingType& rightPrimal);

  private:
    const INDEX numberOfLabels_;    
  };

  DiscreteTomographyMessageCountingPairwise::DiscreteTomographyMessageCountingPairwise(INDEX numberOfLabels)
    : numberOfLabels_(numberOfLabels) {
    assert(numberOfLabels_ > 1);
  }

  template<typename LEFT_FACTOR, typename G1, typename G3>
  void DiscreteTomographyMessageCountingPairwise::SendMessageToRight(LEFT_FACTOR* const f_left, const G1& repam_left, G3& msg, const REAL omega){
    assert(msg.size() == pow(numberOfLabels_,2));
    assert(repam_left.size() == pow(numberOfLabels_,2));

    for(INDEX i=0;i<pow(numberOfLabels_,2);i++){
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
	  if( m > val ){ m = val; }
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
    if( msg > -std::numeric_limits<REAL>::max() &&
	msg < std::numeric_limits<REAL>::max()){ repam[msg_dim] += msg; } 
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
    if( msg > -std::numeric_limits<REAL>::max() &&
	msg < std::numeric_limits<REAL>::max())
      {
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
    
}


#endif // LP_MP_DT_COUNTING_MESSAGE_HXX
