#ifndef LP_MP_MULTICUT_ODD_WHEEL_HXX
#define LP_MP_MULTICUT_ODD_WHEEL_HXX

#include "LP_MP.h"
#include "multicut_triplet_factor.hxx"
#include "permutation.hxx"

namespace LP_MP {

// solve odd wheel constraint by noting that any odd-wheel has tree-width 3, hence we can use separators of size four and directly connect them to each other. Other possiblity, use intermediate separators of size three.(pursued here)
// The separators can be made triplets as well

class MulticutTripletPlusSpokeFactor {
public:
   MulticutTripletPlusSpokeFactor() {}
   // labelings: first are the triplet edges, than the spoke.
   // the triplet edges are ordered as (n1,n2), (n1,centerNode), (n2,centerNode), (centerNode,spokeNode), where n1,n2 are the nodes opposite the spokeNode and n1<n2
   // 0110, 1010, 1100, 1110, 0111, 1011, 1101, 1111, 0001

   using TripletPlusSpokeEdges = std::array<std::array<INDEX,2>,4>;
   static TripletPlusSpokeEdges SortEdges(INDEX n1, INDEX n2, const INDEX centerNode, const INDEX spokeNode)
   {
      if(n2>n1) { 
         std::swap(n1,n2);
      }
      assert(n1 != n2 && n1 != centerNode && n1 != spokeNode);
      assert(            n2 != centerNode && n2 != spokeNode);
      assert(centerNode != spokeNode);
      TripletPlusSpokeEdges te{{{std::min(n1,centerNode), std::max(n1,centerNode)}, {std::min(n2,centerNode),std::max(n2,centerNode)}, {n1,n2}, {std::min(centerNode,spokeNode), std::max(centerNode,spokeNode)}}};
      std::sort(te.begin(), te.begin()+3); // the spoke edge is always the last one
      assert(te.back()[0] == std::min(centerNode,spokeNode) && te.back()[1] == std::max(centerNode,spokeNode));
      return te;
   }

   constexpr static INDEX size() { return 9; }
   template<typename REPAM_ARRAY>
   void MaximizePotential(const REPAM_ARRAY& repam) {};
   template<typename REPAM_ARRAY>
   REAL LowerBound(const REPAM_ARRAY& repam) const
   {
      REAL x = 0.0;
      for(INDEX i=0; i<size(); ++i) {
         x = std::min(x,repam[i]);
      }
      return x;
   }
   template<typename REPAM_ARRAY, typename PRIMAL>
   REAL EvaluatePrimal(const REPAM_ARRAY& repam, const PRIMAL primal) const
   {
      assert(repam.size() == 9);
      return 3e13;
   }

private:
};

// left factor is Triplet, right one is TripletPlusSpoke
// TripletPlusSpoke: the message is applied as follows: two edges are affected, one edge not being the spoke, and the spoke
// Triplet:  Let the corresponding labelings be 01,10,11.
// format: (triplet to triplet edge, triplet to spoke)
// The labeling 00 is not needed, as it only refers to the valid labeling 000 of the triplet, which has cost zero
class MulticutTripletPlusSpokeMessage {
public:
   // i is  the shared triplet edge index in the triplet factor
   // s is the shared spoke edge index in the triplet factor
   // j is the shared triplet edge index in the TripletPlusSpokeFactor
   MulticutTripletPlusSpokeMessage(const INDEX i, const INDEX s, const INDEX j) : i_(i), s_(s), j_(j)
   {
      assert(i_ < 3);
      assert(s_ < 3);
      assert(i_ != s_);
      std::cout << "index not taken = " << 3 - i - s << "\n"; // 3 - i - s == 2 always
      assert(j_ < 3);
   } 
   //MulticutTripletPlusSpokeMessage(Permutation<3> p, const INDEX j) : p_(p), j_(j) {}

   constexpr static INDEX size() { return 3; }

   // send message from TripletPlusSpokeFactor
   template<typename LEFT_FACTOR, typename G1, typename G3>
   void SendMessageToRight(LEFT_FACTOR* const l, const G1& leftPot, G3& msg, const REAL omega)
   {
      assert(leftPot.size() == 4);
      assert(msg.size() == 3);
      msg[0] -= omega*leftPot[i_]; 
      msg[1] -= omega*leftPot[s_];
      msg[2] -= omega*std::min(leftPot[3-i_-s_], leftPot[3]);
      //msg[0] += omega*leftPot[p_[0]];
      //msg[1] += omega*leftPot[p_[1]];
      //msg[2] += omega*std::min(leftPot[p_[2]], leftPot[3]);
   }

   template<typename RIGHT_FACTOR, typename G1, typename G2>
   void ReceiveMessageFromRight(RIGHT_FACTOR* const r, const G1& rightPot, G2& msg) 
   {
      assert(rightPot.size() == 9);
      assert(msg.size() == 3);
      const REAL x = std::min(rightPot[j_],0.0); // this entry and label 0000 is not covered by the message, hence it has to be substracted
      msg[0] -= std::min(rightPot[j_+4],rightPot[8]) - x;
      msg[1] -= std::min(rightPot[3],std::min(rightPot[(j_+1)%3],rightPot[(j_+2)%3])) - x;
      msg[2] -= std::min(std::min(rightPot[(j_+1)%3 + 4],rightPot[(j_+2)%3 + 4]),rightPot[7]) - x;
   }


   // for such factors and messages, it might be advantageous to templatize reparametrization by msg_dim. Do zrobienia.
   template<typename G>
   void RepamLeft(G& repam, const REAL msg, const INDEX msg_dim)
   {
      assert(msg_dim < 3);
      if(msg_dim == 0) { // 01, i.e. complete to 101, 011 or 110, depending on i_ = 0,1,2.
         repam[i_] += msg;
         //repam[p_[0]] -= msg;
      } else if(msg_dim == 1) { // 10, 110, 110, 011, depending on i_ = 0,1,2.
         repam[s_] += msg;
         //repam[p_[1]] -= msg;
      } else { // 11, complete to 011, 101 od 110, depending on i_+s_, and to 111
         repam[3-i_-s_] += msg;
         //repam[p_[2]] -= msg;
         repam[3] += msg; // 111 // (msg_dim==2)*msg
      }
   }

   // 0110, 1010, 1100, 1110, 0111, 1011, 1101, 1111, 0001
   template<typename G>
   void RepamRight(G& repam, const REAL msg, const INDEX msg_dim)
   {
      assert(msg_dim < 3);
      if(msg_dim == 0) { // label 01, complete to 0001 + {0111 | 1011 | 1101} depending on j_ = 0,1,2
         repam[j_+4] += msg;
         repam[8] += msg;
      } else if(msg_dim == 1) { // label 10, complete to 1110 + {1100,1010 | 1100,0110 | 0110,1010} depending on j_ = 0,1,2
         repam[(j_+1)%3] += msg; 
         repam[(j_+2)%3] += msg;
         repam[3] += msg;
      } else { // label 11, complete to 1111 + {1011,1101 | 0111,1101 | 0111,1011 } depending on j_ = 0,1,2
         repam[(j_+1)%3 + 4] += msg; 
         repam[(j_+2)%3 + 4] += msg;
         repam[7] += msg;
      }
   }
private:
   // a message is parametrized as follows: in the triplet, two edges are affected. One that is connected to a triplet edge in the TripletSpoke factor, and one that is connected to the spoke
   // Let i_ be the edge that is not affected in the triplet
   // for the TripletFactor
   // do zrobienia: let i_ and s_ come from a permutation
   //Permutation<3> p_; // permutation of triplet edges of tripletPlusSpoke to triplet edges.
   const unsigned char i_; // shared edge in triplet factor which links to the triplet edge in TripletPlusSpokeFactor
   const unsigned char s_; // shared edge in triplet factor which links to the spoke edge in TripletPlusSpokeFactor
   const unsigned char j_; // edge in triplet of spoke that !is! shared
   //const INDEX i_ : 2; // shared edge in triplet factor which links to the triplet edge in TripletPlusSpokeFactor
   //const INDEX s_ : 2; // shared edge in triplet factor which links to the spoke edge in TripletPlusSpokeFactor
   //const INDEX j_ : 2; // edge in triplet of spoke that !is! shared
   // do zrobienia: check if not using bitfields is faster due to easier execution.
};

// message between triplet and TripletPlusSpoke, when TripletPlusSpoke completely covers the triplet
// left factor is triplet, right one TripletPlusSpoke
class MulticutTripletPlusSpokeCoverMessage {
public:
   MulticutTripletPlusSpokeCoverMessage(Permutation<3> p) : p_(p) {}
   //MulticutTripletPlusSpokeCoverMessage(std::array<INDEX,3> tripletIndices, std::array<INDEX,3> tripletPlusSpokeIndices) :  p_(tripletIndices, tripletPlusSpokeIndices) {}
   MulticutTripletPlusSpokeCoverMessage(const INDEX i1, const INDEX i2, const INDEX i3) : p_({i1,i2,i3})
   {
      //assert that {i1,i2,i3} = {0,1,2}
      assert(i1 == 0 && i2 == 1 && i3 == 2);
      assert(i1 < 3);
      assert(i2 < 3);
      assert(i3 < 3);
      assert(i1 != i2 && i2 != i3 && i1 != i3);
      assert(false);
      //p_.Invert();
   }
   ~MulticutTripletPlusSpokeCoverMessage()
   {
      static_assert(MulticutTripletFactor::size() == 4,"");
   }

   constexpr static INDEX size() { return 4; } 

   // send message from TripletPlusSpokeFactor
   template<typename LEFT_FACTOR, typename G1, typename G3>
   void SendMessageToRight(LEFT_FACTOR* const l, const G1& leftPot, G3& msg, const REAL omega)
   {
      assert(leftPot.size() == 4);
      assert(msg.size() == 4);

      msg[0] -= omega*leftPot[p_[0]];
      msg[1] -= omega*leftPot[p_[1]];
      msg[2] -= omega*leftPot[p_[2]];
      msg[3] -= omega*leftPot[3];
   }

   template<typename RIGHT_FACTOR, typename G1, typename G2>
   void ReceiveMessageFromRight(RIGHT_FACTOR* const r, const G1& rightPot, G2& msg) 
   {
      assert(rightPot.size() == 9);
      assert(msg.size() == 4);
      for(INDEX i=0; i<4; ++i) {
         msg[i] -= std::min(rightPot[i], rightPot[i+4]) - std::min(rightPot[8],0.0);
      }
   }

   // for factors and messages like this one, it might be advantageous to templatize reparametrization by msg_dim. Do zrobienia.
   template<typename G>
   void RepamLeft(G& repam, const REAL msg, const INDEX msg_dim)
   {
      assert(msg_dim < 4);
      if(msg_dim < 3) {
         repam[p_[msg_dim]] += msg;
      } else {
         repam[3] += msg;
      }
   }

   // 0110, 1010, 1100, 1110, 0111, 1011, 1101, 1111, 0001
   template<typename G>
   void RepamRight(G& repam, const REAL msg, const INDEX msg_dim)
   {
      assert(msg_dim < 4);
      repam[msg_dim] += msg;
      repam[msg_dim + 4] += msg;
   }
private:
   // the triplet in the TripletPlusSpokeFactor can have edges arranged in arbitrary order. p_ permutes those
   Permutation<3> p_;
};

} // end namespace LP_MP

#endif // LP_MP_MULTICUT_ODD_WHEEL_HXX
