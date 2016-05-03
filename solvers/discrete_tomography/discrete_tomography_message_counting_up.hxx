#ifndef LP_MP_DT_COUNTING_UP_MESSAGE_HXX
#define LP_MP_DT_COUNTING_UP_MESSAGE_HXX

#include "LP_MP.h"

namespace LP_MP {

  enum class MessageSending { SRMP, MPLP };
  template<MessageSending MST = MessageSending::SRMP>
  class DiscreteTomographyMessageCountingUP{

  public:

    DiscreteTomographyMessageCountingUP(){};

    // UP
    template<typename RIGHT_FACTOR, typename G1, typename G2, MessageSending MST_TMP = MST>
    typename std::enable_if<MST_TMP == MessageSending::SRMP,void>::type
    ReceiveMessageFromRight(RIGHT_FACTOR* const r, const G1& rightPot, G2& msg);

    template<typename LEFT_FACTOR, typename RIGHT_FACTOR, typename G1, typename G2, typename G3, MessageSending MST_TMP = MST>
    typename std::enable_if<MST_TMP == MessageSending::SRMP,void>::type
    SendMessageToRight(LEFT_FACTOR* const l, RIGHT_FACTOR* const r, const G1& leftPot, const G2& rightPot, G3& msg, const REAL omega);

    // DOWN
    template<typename LEFT_FACTOR, typename G1, typename G2, MessageSending MST_TMP = MST>
    typename std::enable_if<MST_TMP == MessageSending::SRMP,void>::type
    ReceiveMessageFromLeft(LEFT_FACTOR* const r, const G1& rightPot, G2& msg);

    template<typename LEFT_FACTOR, typename RIGHT_FACTOR, typename G1, typename G2, typename G3, MessageSending MST_TMP = MST>
    typename std::enable_if<MST_TMP == MessageSending::SRMP,void>::type
    SendMessageToLeft(LEFT_FACTOR* const l, RIGHT_FACTOR* const r, const G1& leftPot, const G2& rightPot, G3& msg, const REAL omega);

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
