#ifndef LP_MP_MARGINAL_SUMMATION_MESSAGE_HXX
#define LP_MP_MARGINAL_SUMMATION_MESSAGE_HXX

#include "instances.inc"

namespace LP_MP {

// implement constraint mu(0) = nu(0) + ... + nu(n-1)
//                      mu(1) = nu(n)
// left factor is simplex mu with 2 entries, right factor is simplex nu with n entries

// do zrobienia: generalize this to other summmation cases
class MarginalSummationMessage {
public:
   template<typename LEFT_POT, typename MSG_ARRAY>
   void MakeLeftFactorUniform(const LEFT_POT& leftPot, MSG_ARRAY& msg, const double omega = 1.0) {
      msg[0] += omega*leftPot[0];
      msg[1] += omega*leftPot[1];
   }
   template<typename RIGHT_POT, typename MSG_ARRAY>
   void MakeRightFactorUniform(const RIGHT_POT& rightPot, MSG_ARRAY& msg, const double omega = 1.0) {
      REAL rightPotMin = std::numeric_limits<double>::max();
      for(INDEX i=0; i<rightPot.size()-1; ++i) {
         rightPotMin = std::min(rightPotMin, rightPot[i]);
      }
      msg[0] -= omega*rightPotMin;
      msg[1] -= omega*rightPot[rightPot.size()-1];
   }

   template<typename RIGHT_FACTOR, typename G1, typename G2>
   void ReceiveMessageFromRight(RIGHT_FACTOR* const r, const G1& rightPot, G2& msg)
   {
      MakeRightFactorUniform(rightPot,msg);
   }

   template<typename LEFT_FACTOR, typename G1, typename G2>
   void ReceiveMessageFromLeft(LEFT_FACTOR* l, const G1& leftPot, G2& msg)
   {
      MakeLeftFactorUniform(leftPot,msg);
   }

   template<typename LEFT_FACTOR, typename RIGHT_FACTOR, typename G1, typename G2, typename G3>
   void SendMessageToRight(LEFT_FACTOR* const l, RIGHT_FACTOR* const r, const G1& leftPot, const G2& rightPot, G3& msg, const REAL omega)
   { 
      MakeLeftFactorUniform(leftPot,msg,omega);
   }

   template<typename LEFT_FACTOR, typename RIGHT_FACTOR, typename G1, typename G2, typename G3>
   void SendMessageToLeft(LEFT_FACTOR* const l, RIGHT_FACTOR* const r, const G1& leftPot, const G2& rightPot, G3& msg, const REAL omega)
   {
      MakeRightFactorUniform(rightPot,msg,omega);
   }

   template<typename G>
   void RepamLeft(G& leftRepamPot, const REAL msg, const INDEX dim) {
      assert(dim == 0 || dim == 1);
      leftRepamPot[dim] = leftRepamPot[dim] - msg;
   }
   template<typename G>
   void RepamRight(G& rightRepamPot, const REAL msg, const INDEX dim) {
      assert(dim == 0 || dim == 1);
      if(dim == 0) {
         for(INDEX i=0; i<rightRepamPot.size()-1; ++i) {
            rightRepamPot[i] = rightRepamPot[i] + msg;
         }
      } else {
         rightRepamPot[rightRepamPot.size() - 1] = rightRepamPot[rightRepamPot.size() - 1] + msg;
      }
   }

private:
   
};
}

#endif // LP_MP_MARGINAL_SUMMATION_MESSAGE_HXX



