#ifndef LP_MP_DT_COUNTING_UP_MESSAGE_HXX
#define LP_MP_DT_COUNTING_UP_MESSAGE_HXX

#include "LP_MP.h"

namespace LP_MP {

  class DiscreteTomographyMessageCountingUP{

  public:

    DiscreteTomographyMessageCountingUP(){};

    // UP
    template<typename RIGHT_FACTOR, typename G1, typename G2>
    void ReceiveMessageFromRight(RIGHT_FACTOR* const r, const G1& rightPot, G2& msg);

    template<typename LEFT_FACTOR, typename RIGHT_FACTOR, typename G1, typename G3>
    void SendMessageToRight(LEFT_FACTOR* const l, const G1& leftPot, G3& msg, const REAL omega);

    // DOWN
    template<typename LEFT_FACTOR, typename G1, typename G2>
    void ReceiveMessageFromLeft(LEFT_FACTOR* const r, const G1& rightPot, G2& msg);

    template<typename LEFT_FACTOR, typename RIGHT_FACTOR, typename G2, typename G3>
    void SendMessageToLeft(RIGHT_FACTOR* const r, const G2& rightPot, G3& msg, const REAL omega);

    /*------*/
    template<typename G>
    void RepamLeft(G& repamPot, const REAL msg, const INDEX msg_dim);
    
    template<typename G>
    void RepamRight(G& repamPot, const REAL msg, const INDEX msg_dim);

    //void ComputeRightFromLeftPrimal(const bool leftPrimal, MulticutTripletFactor::LabelingType& rightPrimal);
    
  private:
    
    
  };

}

#endif // LP_MP_DT_COUNTING_UP_MESSAGE_HXX
