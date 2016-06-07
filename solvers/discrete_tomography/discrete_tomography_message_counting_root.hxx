#ifndef LP_MP_DT_COUNTING_MESSAGE_HXX
#define LP_MP_DT_COUNTING_MESSAGE_HXX

#include "LP_MP.h"
#include "minConv.hxx"
#include <math.h> 

namespace LP_MP {
  namespace discrete_tomo{
    
    using MinConv = MinConv<std::function<REAL(INDEX)>,REAL,INDEX>;

    class DiscreteTomographyMessageCountingRoot{

    public:

      DiscreteTomographyMessageCountingRoot(INDEX numberOfLabels,INDEX numberOfVars);

      // RIGHT -> Unary
      // LEFT  -> Ternary

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
    
      void setRHS(INDEX rhs){ isRhs_ = true; rhs_ = rhs; }
      void setHCS(bool hcs){ hardCS_ = hcs; };
      
      
    private:
      const INDEX numberOfLabels_;
      INDEX msgSize_;
      INDEX leftSize_,rightSize_,upSize_,regSize_;
      INDEX rhs_;
      bool isRhs_ = false, hardCS_ = false;
    
    };

    DiscreteTomographyMessageCountingRoot::DiscreteTomographyMessageCountingRoot(INDEX numberOfLabels,INDEX numberOfVars);
    : numberOfLabels_(numberOfLabels), numberOfVars_(numberOfVars) {
      assert(numberOfLabels_ > 1);
      assert(numberOfVars_ > 1)
    }

    template<typename LEFT_FACTOR, typename G1, typename G3>
    void DiscreteTomographyMessageCountingRoot::SendMessageToRight(LEFT_FACTOR* const f_left, const G1& repam_left, G3& msg, const REAL omega){
      assert(msg.size() == f_left.getSize(f_left::up));
      assert(repam_left.size() == (f_left.getSize(f_left::up) + f_left.getSize(f_left::left) +
				    f_left.getSize(f_left::right) + f_left.getSize(f_left::reg)));

      INDEX up_size = f_left.getSize(f_left::up)/pow(numberOfLabels_,2);
      INDEX right_size = f_left.getSize(f_left::right)/pow(numberOfLabels_,2);
      INDEX left_size = f_left.getSize(f_left::left)/pow(numberOfLabels_,2);

      if( isRhs_ && up_size > rhs_ + 1){ up_size = rhs_ + 1; }
      op = [&](INDEX i,INDEX j){ return (i+j < up_size ) ? i+j : up_size;  };
      
      for(INDEX i=0;i<pow(numberOfLabels_,4);i++){
        INDEX idx = i;
        INDEX a = idx % numberOfLabels_;
        idx = ( idx - a )/numberOfLabels_;
        INDEX b = idx % numberOfLabels_;
        idx = ( idx - b )/numberOfLabels_;
        INDEX c = idx % numberOfLabels_;
        idx = ( idx - c )/numberOfLabels_;
        INDEX d = idx % numberOfLabels_;

        //auto z_up = [&](INDEX k){ return repam_left[a + d*numberOfLabels_ + k*pow(numberOfLabels_,2)];  };
        auto z_left = [&](INDEX k){ return repam_left[f_left.getSize(f_left::up) + a + b*numberOfLabels_ + k*pow(numberOfLabels_,2)];  };
        auto z_right = [&](INDEX k){ return repam_left[f_left.getSize(f_left::up) + f_left.getSize(f_left::left) + c + d*numberOfLabels_ + k*pow(numberOfLabels_,2)];  };
        auto z_reg = repam_left[f_left.getSize(f_left::up) + f_left.getSize(f_left::left) + f_left.getSize(f_left::right) + b + c*numberOfLabels_];
            
        MinConv mc(z_left,z_right,left_size,right_size,up_size);
        mc.CalcConv(op);

        for(INDEX k=0;k<up_size;k++){
          REAL val = mc.getConv(k);
          INDEX kidx = a + numberOfLabels_*d + k*pow(numberOfLabels_,2);
          assert(kidx < (up_size*pow(numberOfLabels_,2)));
          msg[kidx] -= omega*val;
        }
      }
    }
    
    template<typename RIGHT_FACTOR, typename G1, typename G2>
    void DiscreteTomographyMessageCountingRoot::ReceiveMessageFromRight(RIGHT_FACTOR* const f_right, const G1& repam_right, G2& msg){
      assert(f_right.getSize() == repam_right.size());
      assert(f_right.getSize() == msg.size());
    
      if( hardCS_ ){
        assert(isRhs_); assert(rhs_ >= 0); assert(rhs_ < (f_right.getSize()/pow(numberOfLabels_,2)));
        for(INDEX i=0;i<pow(numberOfLabels_,2);i++){
          INDEX a = i % numberOfLabels_;
          INDEX b = ((i - a)/numberOfLabels_) % numberOfLabels_;
          msg[a + numberOfLabels_*b + i*pow(numberOfLabels_,2)] -= repam_right[a + numberOfLabels_*b + i*pow(numberOfLabels_,2)];
        }
      }
      else{
        
      }
      
    }

    template<DIRECTION DR>
    template<typename G>
    void DiscreteTomographyMessageCountingRoot::RepamLeft(G& repam, const REAL msg, const INDEX msg_dim){

      auto f = repam.GetFactor();
      assert( repam.size() == (f.getSize(f::left) + f.getSize(f::right) + f.getSize(f::up) + f.getSize(f::reg))); // <-- wie kann man das machen?

      assert(msg_dim < f.getSize(f::up));
      repam[msg_dim] += msg;
      
    }

    template<DIRECTION DR>    
    template<typename G>
    void DiscreteTomographyMessageCountingRoot::RepamRight(G& repam, const REAL msg, const INDEX msg_dim){

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
