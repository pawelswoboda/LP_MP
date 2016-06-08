#ifndef LP_MP_DT_COUNTING_MESSAGE_HXX
#define LP_MP_DT_COUNTING_MESSAGE_HXX

#include "LP_MP.h"
#include "minConv.hxx"
#include <math.h> 

namespace LP_MP {
  namespace discrete_tomo{
    
    using MinConv = MinConv<std::function<REAL(INDEX)>,REAL,INDEX>;

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

    DiscreteTomographyMessageCountingPairwise::DiscreteTomographyMessageCountingPairwise(INDEX numberOfLabels);
    : numberOfLabels_(numberOfLabels) {
      assert(numberOfLabels_ > 1);
    }

    template<typename LEFT_FACTOR, typename G1, typename G3>
    void DiscreteTomographyMessageCountingPairwise::SendMessageToRight(LEFT_FACTOR* const f_left, const G1& repam_left, G3& msg, const REAL omega){
      assert(msg.size() == pow(numberOfLabels_,2));
      assert(repam_right.size() == pow(numberOfLabels_,2));

      for(INDEX i=0;i<pow(numberOfLabels_,2);i++){
	msg[i] -= omega*repam_left[i];
      }     
    }
    
    template<typename RIGHT_FACTOR, typename G1, typename G2>
    void DiscreteTomographyMessageCountingPairwise::ReceiveMessageFromRight(RIGHT_FACTOR* const f_right, const G1& repam_right, G2& msg){

      assert(msg.size() == pow(numberOfLabels_,2));
      assert(repam_right.size() == (f_right.getSize(f_right::up) + f_right.getSize(f_right::left) +
				    f_right.getSize(f_right::right) + f_right.getSize(f_right::reg)));

      INDEX left_size = f_right.getSize(f_right::left)/pow(numberOfLabels_,2);
      INDEX right_size = f_right.getSize(f_right::right)/pow(numberOfLabels_,2);
      INDEX up_size = f_right.getSize(f_right::up)/pow(numberOfLabels_,2);
	
      op = [&](INDEX i,INDEX j){ return i+j; };
	
      for(INDEX i=0;i<pow(numberOfLabels_,2);i++){
	REAL m = std::numeric_limits<REAL>::max();
	INDEX b = i % numberOfLabels_;
	INDEX c = ((i-b)/numberOfLabels_) % numberOfLabels_;
	REAL reg = repam_right[f_right.getSize(f_right::up) + f_right.getSize(f_right::left) + f_right.getSize(f_right::right) + b + c*numberOfLabels_];
	for(INDEX j=0;j<pow(numberOfLabels_,2);j++){
	  INDEX a = j % numberOfLabels_;
	  INDEX d = ((j-a)/numberOfLabels_) % numberOfLabels_;

	  auto z_up = [&](INDEX k){ return repam_right[a + d*numberOfLabels_ + k*pow(numberOfLabels_,2)];  };
	  auto z_left = [&](INDEX k){ return repam_right[f_right.getSize(f_right::up) + a + b*numberOfLabels_ + k*pow(numberOfLabels_,2)];  };
	  auto z_right = [&](INDEX k){ return repam_right[f_right.getSize(f_right::up) + f_right.getSize(f_right::left) + c + d*numberOfLabels_ + k*pow(numberOfLabels_,2)];  };
	  	  
	  MinConv mc(z_left,z_right,left_size,right_size,up_size);
	  mc.CalcConv(op);

	  for(INDEX k=0;k<up_size;k++){
	    assert( k == (mc.getIdxA(k) + mc.getIdxB(k)));
	    REAL val = mc.getConv(k) + reg;
	    if( m > val ){ m = val; }
	  }
	  msg[i] -= m;
	}
      }    
    }

    template<typename G>
    void DiscreteTomographyMessageCountingPairwise::RepamLeft(G& repam, const REAL msg, const INDEX msg_dim){

      auto f = repam.GetFactor();
      assert( repam.size() == pow(numberOfLabels_,2));

      assert(msg_dim < pow(numberOfLabels_,2));
      repam[msg_dim] += msg;
      
    }

    template<typename G>
    void DiscreteTomographyMessageCountingPairwise::RepamRight(G& repam, const REAL msg, const INDEX msg_dim){

      auto f = repam.GetFactor();
      assert( repam.size() == (f.getSize(f::up) + f.getSize(f::left) + f.getSize(f::right) + f.getSize(f::reg))); // <-- wie kann man das machen?

      assert(msg_dim < pow(numberOfLabels_,2));
      repam[f.getSize(f::up) + f.getSize(f::left) + f.getSize(f::right) + msg_dim] +=  msg;
    }
    
  }
}

#endif // LP_MP_DT_COUNTING_MESSAGE_HXX
