#ifndef LP_MP_MIN_COST_FLOW_MESSAGE_HXX
#define LP_MP_MIN_COST_FLOW_MESSAGE_HXX

// message between variable 0 <= x <= 1 and edge of min cost flow graph with capacities [0,1]
// left factor is simplex with variable index, right factor is min cost flow factor and also variable index

// do zrobienia: apply inverse LP theory to "SendMessagesToLeft"

namespace LP_MP {

class MinCostFlowMessage {
public:
   MinCostFlowMessage(const INDEX simplexIdx, const INDEX minCostFlowIdx)
      :
      simplexIdx_(simplexIdx),
      minCostFlowIdx_(minCostFlowIdx)
   {}

   template<typename RIGHT_FACTOR, typename G1, typename G2>
   void ReceiveMessageFromRight(RIGHT_FACTOR* const r, const G1& rightPot, G2& msg)
   {
      // make cost in edge such that having flow 0 or 1 will give the same cost
      // breakpoint is the value by which we can change the current cost in the flow network at most while still letting the current flow be optimal
      const REAL breakpoint = r->GetBreakpoint(rightPot, minCostFlowIdx_); // here we assume that rightPot is the current reparametrization
      msg[0] += breakpoint;
   }

   template<typename LEFT_FACTOR, typename G1, typename G2>
   void ReceiveMessageFromLeft(LEFT_FACTOR* l, const G1& leftPot, G2& msg)
   {
      //make left factor uniform
      REAL minVal = std::numeric_limits<REAL>::max();
      for(INDEX i=0; i<leftPot.size(); ++i) {
         if(i!= simplexIdx_) {
            minVal = std::min(leftPot[i], minVal);
         }
      }
      msg[0] += (leftPot[simplexIdx_] - minVal);
   }

   template<typename LEFT_FACTOR, typename RIGHT_FACTOR, typename G1, typename G2, typename G3>
   void SendMessageToRight(LEFT_FACTOR* const l, RIGHT_FACTOR* const r, const G1& leftPot, const G2& rightPot, G3& msg, const REAL omega)
   {
      assert(simplexIdx_ == 0); // do zrobienia: for now
      //make left factor uniform
      REAL minVal = std::numeric_limits<REAL>::max();
      for(INDEX i=0; i<leftPot.size(); ++i) {
         if(i!= simplexIdx_) {
            minVal = std::min(leftPot[i], minVal);
         }
      }
      msg[0] += omega*(leftPot[simplexIdx_] - minVal);
   }

   // note: computing this message would not give us much, as very many messages are computed on the graph which would need sharing, thus making each change extremely small
   // additionally, this would entail copying of the whole cost of the network, which would be slow. Use Inverse LP theory
   /*
   template<typename LEFT_FACTOR, typename RIGHT_FACTOR, typename G1, typename G2, typename G3>
   void SendMessageToLeft(LEFT_FACTOR* const l, RIGHT_FACTOR* const r, const G1& leftPot, const G2& rightPot, G3& msg, const REAL omega)
   {}
   */

   template<typename REPAM_ARRAY>
   void RepamLeft(REPAM_ARRAY& leftRepamPot, const REAL msg, const INDEX dim)
   {  
      assert(dim == 0);
      leftRepamPot[simplexIdx_] -= msg;
   }
   template<typename REPAM_ARRAY>
   void RepamRight(REPAM_ARRAY& rightRepamPot, const REAL msg, const INDEX dim)
   {
      assert(dim == 0);
      rightRepamPot[minCostFlowIdx_] += msg;
   }

private:
   const INDEX simplexIdx_;
   const INDEX minCostFlowIdx_;
};

} // end namespace LP_MP

#endif // LP_MP_MIN_COST_FLOW_MESSAGE_HXX
