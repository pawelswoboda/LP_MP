#ifndef LP_MP_MULTICUT_ODD_WHEEL_HXX
#define LP_MP_MULTICUT_ODD_WHEEL_HXX

#include "LP_MP.h"
#include "multicut_triplet_factor.hxx"

namespace LP_MP {

// solve odd wheel constraint by noting that any odd-wheel has tree-width 3, hence we can use separators of size four and directly connect them to each other. Other possiblity, use intermediate separators of size three, which we pursue here.
// The separators can be made triplets as well

class MulticutTripletPlusSpokeFactor : public std::array<REAL,9> {
public:
   MulticutTripletPlusSpokeFactor() {
      std::fill(this->begin(), this->end(), 0.0);
   }
   // labelings: first are the triplet edges, than the spoke.
   // the triplet edges are ordered as (n1<n2), (n1,centerNode), (n2,centerNode), (centerNode,spokeNode), where n1,n2 are the nodes opposite the spokeNode and n1<n2
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
      TripletPlusSpokeEdges te{{{n1,n2}, {std::min(n1,centerNode), std::max(n1,centerNode)}, {std::min(n2,centerNode),std::max(n2,centerNode)}, {std::min(centerNode,spokeNode), std::max(centerNode,spokeNode)}}};
      std::sort(te.begin(), te.begin()+3); // the spoke edge is always the last one
      assert(te.back()[0] == std::min(centerNode,spokeNode) && te.back()[1] == std::max(centerNode,spokeNode));
      return te;
   }

   constexpr static INDEX size() { return 9; }
   REAL LowerBound() const
   {
      return std::min(0.0, *std::min_element(this->begin(), this->end()));
   }
   REAL EvaluatePrimal() const
   {
      const auto sum = std::count(primal_.begin(), primal_.begin()+3, true);
      if(sum == 1) { 
         return std::numeric_limits<REAL>::infinity();
      }

      if(!primal_[3]) {
         if(!primal_[0] && primal_[1] && primal_[2]) {
            return (*this)[0];
         } else if(primal_[0] && !primal_[1] && primal_[2]) {
            return (*this)[1];
         } else if(primal_[0] && primal_[1] && !primal_[2]) {
            return (*this)[2];
         } else if(primal_[0] && primal_[1] && primal_[2]) {
            return (*this)[3];
         } else if(!primal_[0] && !primal_[1] && !primal_[2]) {
            return 0.0;
         } else {
            assert(false);
            return std::numeric_limits<REAL>::infinity();
         }
      } else {
         if(!primal_[0] && primal_[1] && primal_[2]) {
            return (*this)[4];
         } else if(primal_[0] && !primal_[1] && primal_[2]) {
            return (*this)[5];
         } else if(primal_[0] && primal_[1] && !primal_[2]) {
            return (*this)[6];
         } else if(primal_[0] && primal_[1] && primal_[2]) {
            return (*this)[7];
         } else if(!primal_[0] && !primal_[1] && !primal_[2]) {
            return (*this)[8];
         } else {
            assert(false);
            return std::numeric_limits<REAL>::infinity();
         }
      }
   }

   void init_primal() {}
   template<typename ARCHIVE> void serialize_dual(ARCHIVE& ar) { ar( *static_cast<std::array<REAL,9>*>(this) ); }
   template<typename ARCHIVE> void serialize_primal(ARCHIVE& ar) { ar( primal_ ); }

   auto get_primal() const { return primal_; }
   auto& set_primal() { return primal_; }

private:
   std::array<bool,4> primal_;
};

// left factor is Triplet, right one is TripletPlusSpoke
// TripletPlusSpoke: the message is applied as follows: two edges are affected, one edge not being the spoke, and the spoke
// Triplet:  Let the corresponding labelings be 01,10,11.
// format: (triplet to triplet edge, triplet to spoke)
// The labeling 00 is not needed, as it only refers to the valid labeling 000 of the triplet, which has cost zero
//
//
// msg[0] joins first triplet edge with the non-spoke edge, value = 0, other edge value 1
// msg[1] joins the second triplet edge with the spoke edge, value = 0, other edge value 1
// msg[2] corresponds to first and second shared triplet edge taking value 1 both
class MulticutTripletPlusSpokeMessage {
public:
   // first four entries are coordinates of TripletPlusSpoke, last three ones are coordinates of Triplet
   MulticutTripletPlusSpokeMessage(const INDEX n1, const INDEX n2, const INDEX centerNode, const INDEX spokeNode, const INDEX i1, const INDEX i2, const INDEX i3) 
   {
      assert(i1 < i2 && i2 < i3);
      if(std::min(n1,centerNode) == i1 && std::max(n1,centerNode) == i2) {
         tripletPlusSpokeEdge_ = 1;
         tripletEdge_[0] = 0;
      } else if(std::min(n1,centerNode) == i1 && std::max(n1,centerNode) == i3) {
         tripletPlusSpokeEdge_ = 1;
         tripletEdge_[0] = 1;
      } else if(std::min(n1,centerNode) == i2 && std::max(n1,centerNode) == i3) {
         tripletPlusSpokeEdge_ = 1;
         tripletEdge_[0] = 2;
      } else if(std::min(n2,centerNode) == i1 && std::max(n2,centerNode) == i2) {
         tripletPlusSpokeEdge_ = 2;
         tripletEdge_[0] = 0;
      } else if(std::min(n2,centerNode) == i1 && std::max(n2,centerNode) == i3) {
         tripletPlusSpokeEdge_ = 2;
         tripletEdge_[0] = 1;
      } else if(std::min(n2,centerNode) == i2 && std::max(n2,centerNode) == i3) {
         tripletPlusSpokeEdge_ = 2;
         tripletEdge_[0] = 2;
      } else {
         assert(false);
      }
      
      if(std::min(centerNode,spokeNode) == i1 && std::max(centerNode,spokeNode) == i2) {
         tripletEdge_[1] = 0;
      } else if(std::min(centerNode,spokeNode) == i1 && std::max(centerNode,spokeNode) == i3) {
         tripletEdge_[1] = 1;
      } else if(std::min(centerNode,spokeNode) == i2 && std::max(centerNode,spokeNode) == i3) {
         tripletEdge_[1] = 2;
      } else {
         assert(false);
      }

      //i_ = 0;
      //s_ = 1;
      //j_ = 2;
   }
   // i is  the shared triplet edge index in the triplet factor
   // s is the shared spoke edge index in the triplet factor
   // j is the shared triplet edge index in the TripletPlusSpokeFactor
   /*
   MulticutTripletPlusSpokeMessage(const INDEX i, const INDEX s, const INDEX j) : i_(i), s_(s), j_(j)
   {
      assert(false); // deprecated
      assert(i_ < 3);
      assert(s_ < 3);
      assert(i_ != s_);
      std::cout << "index not taken = " << 3 - i - s << "\n"; // 3 - i - s == 2 always
      assert(j_ < 3);
   } 
   */
   //MulticutTripletPlusSpokeMessage(Permutation<3> p, const INDEX j) : p_(p), j_(j) {}

   constexpr static INDEX size() { return 3; }

   // send message from TripletPlusSpokeFactor
   template<typename LEFT_FACTOR, typename G3>
   void SendMessageToRight(const LEFT_FACTOR& l, G3& msg, const REAL omega)
   {
      assert(l.size() == 4);
      msg[0] -= omega*l[tripletEdge_[0]]; 
      msg[1] -= omega*l[tripletEdge_[1]];
      msg[2] -= omega*std::min(l[3-tripletEdge_[0] - tripletEdge_[1]], l[3]);
      //msg[0] -= omega*leftPot[i_]; 
      //msg[1] -= omega*leftPot[s_];
      //msg[2] -= omega*std::min(leftPot[3-i_-s_], leftPot[3]);
   }

   template<typename RIGHT_FACTOR, typename G2>
   void ReceiveMessageFromRight(const RIGHT_FACTOR& r, G2& msg) 
   {
      assert(r.size() == 9);
      const REAL x = std::min(r[tripletPlusSpokeEdge_],0.0); // this entry and label 0000 is not covered by the message, hence it has to be substracted
      msg[0] -= std::min(r[tripletPlusSpokeEdge_+4],r[8]) - x;
      msg[1] -= std::min(r[3],std::min(r[(tripletPlusSpokeEdge_+1)%3],r[(tripletPlusSpokeEdge_+2)%3])) - x;
      msg[2] -= std::min(std::min(r[((tripletPlusSpokeEdge_+1)%3) + 4],r[((tripletPlusSpokeEdge_+2)%3) + 4]),r[7]) - x;
      //const REAL x = std::min(rightPot[j_],0.0); // this entry and label 0000 is not covered by the message, hence it has to be substracted
      //msg[0] -= std::min(rightPot[j_+4],rightPot[8]) - x;
      //msg[1] -= std::min(rightPot[3],std::min(rightPot[(j_+1)%3],rightPot[(j_+2)%3])) - x;
      //msg[2] -= std::min(std::min(rightPot[(j_+1)%3 + 4],rightPot[(j_+2)%3 + 4]),rightPot[7]) - x;
   }


   // for such factors and messages, it might be advantageous to templatize reparametrization by msg_dim. Do zrobienia.
   template<typename G>
   void RepamLeft(G& repam, const REAL msg, const INDEX msg_dim)
   {
      assert(msg_dim < 3);
      assert(repam.size() == 4);
      if(msg_dim < 2) {
         repam[tripletEdge_[msg_dim]] += msg;
      } else {
         repam[3-tripletEdge_[0] - tripletEdge_[1]] += msg;
         repam[3] += msg;
      }
      return;


      /*
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
      */
   }

   template<typename G>
   void RepamRight(G& repam, const REAL msg, const INDEX msg_dim)
   {
      // 0110, 1010, 1100, 1110, 0111, 1011, 1101, 1111, 0001
      assert(msg_dim < 3);
      assert(repam.size() == 9);
      if(msg_dim == 0) { // 01
         repam[tripletPlusSpokeEdge_ + 4] += msg;
         repam[8] += msg;
      } else if(msg_dim == 1) { // 10
         repam[(tripletPlusSpokeEdge_+1)%3] += msg;
         repam[(tripletPlusSpokeEdge_+2)%3] += msg;
         repam[3] += msg;
      } else { // 11
         repam[((tripletPlusSpokeEdge_+1)%3) + 4] += msg;
         repam[((tripletPlusSpokeEdge_+2)%3) + 4] += msg;
         repam[7] += msg;
      }
      return;


      // 0110, 1010, 1100, 1110, 0111, 1011, 1101, 1111, 0001
      /*
      if(msg_dim == 0) { // label 01, complete to 0001 + {0111 | 1011 | 1101} depending on j_ = 0,1,2
         //repam[j_] += msg; // ???
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
      */
   }

   template<typename LEFT_FACTOR, typename RIGHT_FACTOR>
   void ComputeRightFromLeftPrimal(const LEFT_FACTOR& l, RIGHT_FACTOR& r) const
   {
      r.set_primal()[tripletPlusSpokeEdge_] = l.get_primal()[tripletEdge_[0]];
      r.set_primal()[3] = l.get_primal()[tripletEdge_[1]];
   }

   template<typename LEFT_FACTOR, typename RIGHT_FACTOR>
   bool CheckPrimalConsistency(const LEFT_FACTOR& l, const RIGHT_FACTOR& r) const
   {
      return true;
   }

private:
   INDEX tripletPlusSpokeEdge_; // the triplet edge that tripletPlusSpoke has in common, the other being the spoke edge
   std::array<INDEX,2> tripletEdge_; // the two edges triplet has in common
};

// message between triplet and TripletPlusSpoke, when TripletPlusSpoke completely covers the triplet
// left factor is triplet, right one TripletPlusSpoke
class MulticutTripletPlusSpokeCoverMessage {
public:
   MulticutTripletPlusSpokeCoverMessage(const INDEX n1, const INDEX n2, const INDEX centerNode, const INDEX spokeNode)  
   {
      assert(n1<n2);
      // find permutation such that TripletPlusSpoke and Triplet are mapped correctly onto each other
      // in triplet, the edges are sorted according to values n1,n2,centerNode
      std::array<INDEX,3> tripletNodes{n1,n2,centerNode};
      std::sort(tripletNodes.begin(),tripletNodes.end());
      p_[0] = 0;
      p_[1] = 1;
      p_[2] = 2;
      if(n1 == tripletNodes[0] && n2 == tripletNodes[1]) {
         p_[0] = 0;
      } else if(n1 == tripletNodes[0] && n2 == tripletNodes[2]) {
         p_[0] = 1;
      } else {
         assert(n1 == tripletNodes[1] && n2 == tripletNodes[2]);
         p_[0] = 2;
      }
      if(std::min(n1,centerNode) == tripletNodes[0] && std::max(n1,centerNode) == tripletNodes[1]) {
         p_[1] = 0;
      } else if(std::min(n1,centerNode) == tripletNodes[0] && std::max(n1,centerNode) == tripletNodes[2]) {
         p_[1] = 1;
      } else {
         assert(std::min(n1,centerNode) == tripletNodes[1] && std::max(n1,centerNode) == tripletNodes[2]);
         p_[1] = 2;
      }
      p_[2] = 3 - p_[0] - p_[1];
      assert(p_[0] != p_[1] && p_[0] != p_[2] && p_[1] != p_[2]);
   }
   /*
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
   */
   ~MulticutTripletPlusSpokeCoverMessage()
   {
      static_assert(MulticutTripletFactor::size() == 4,"");
   }

   constexpr static INDEX size() { return 4; } 

   // send message from TripletPlusSpokeFactor
   template<typename LEFT_FACTOR, typename G3>
   void SendMessageToRight(const LEFT_FACTOR& l, G3& msg, const REAL omega)
   {
      assert(l.size() == 4);

      msg[0] -= omega*l[p_[0]]; // 011
      msg[1] -= omega*l[p_[1]]; // 101
      msg[2] -= omega*l[p_[2]]; // 110
      msg[3] -= omega*l[3];     // 111
   }

   template<typename RIGHT_FACTOR, typename G2>
   void ReceiveMessageFromRight(const RIGHT_FACTOR& r, G2& msg) 
   {
      assert(r.size() == 9);
      for(INDEX i=0; i<4; ++i) {
         msg[i] -= std::min(r[i], r[i+4]) - std::min(r[8],0.0);
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

   template<typename LEFT_FACTOR, typename RIGHT_FACTOR>
   void ComputeRightFromLeftPrimal(const LEFT_FACTOR& l, RIGHT_FACTOR& r) const
   {
      r.set_primal()[0] = l.get_primal()[p_[0]];
      r.set_primal()[1] = l.get_primal()[p_[1]];
      r.set_primal()[2] = l.get_primal()[p_[2]];
   }

   template<typename LEFT_FACTOR, typename RIGHT_FACTOR>
   bool CheckPrimalConsistency(const LEFT_FACTOR& l, const RIGHT_FACTOR& r) const
   {
      return true;
   }

private:
   // the triplet in the TripletPlusSpokeFactor can have edges arranged in arbitrary order. p_ permutes those
   //Permutation<3> p_;
   std::array<INDEX,3> p_;
};

} // end namespace LP_MP

#endif // LP_MP_MULTICUT_ODD_WHEEL_HXX
