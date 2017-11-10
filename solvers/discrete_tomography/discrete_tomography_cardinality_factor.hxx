#ifndef LP_MP_DISCRETE_TOMOGRAPHY_CARDINALITY_HXX
#define LP_MP_DISCRETE_TOMOGRAPHY_CARDINALITY_HXX

namespace LP_MP {

#include <algorithm>
#include "vector.hxx"
#include "tropical_convolution.hxx"

// build factors recursively. Share reparametrization between factors
class discrete_tomography_cardinality_factor {
public:
   discrete_tomography_cardinality_factor(const INDEX left_size, const INDEX right_size, const INDEX min_conv_size)
      : left(left_size, 0.0), 
      right(right_size, 0.0),
      min_conv(min_conv_size, 0.0)
   {}

   discrete_tomography_cardinality_factor(const INDEX left_size, const INDEX right_size)
      : left(left_size, 0.0), 
      right(right_size, 0.0),
      min_conv(left_size + right_size - 1, 0.0)
   {}

   void set_reference_to_min_conv(const vector<REAL>& other)
   {
      min_conv.share(other); 
   }

   REAL LowerBound() const
   {
      // compute min convolution of left and right potential;
      vector<REAL> mc(min_conv.size());
      tropical_convolution::min_conv(left.begin(), left.end() ,right.begin(), right.end(), mc.begin(), mc.end());
      REAL x = std::numeric_limits<REAL>::infinity();
      for(INDEX i=0; i<min_conv.size(); ++i) {
         x = std::min(x, normalize(mc[i] - min_conv[i]));
      }
      return x; 
   }

   REAL EvaluatePrimal() const 
   {
      if(!primal_valid()) { return std::numeric_limits<REAL>::infinity(); }
      assert(left[primal_left] + right[primal_right] - min_conv[primal_sum] < std::numeric_limits<REAL>::infinity());
      return normalize(left[primal_left] + right[primal_right] - min_conv[primal_sum]);
   }

   void MaximizePotentialAndComputePrimal()
   {
      if(primal_sum == std::numeric_limits<INDEX>::max() && primal_left == std::numeric_limits<INDEX>::max() && primal_right == std::numeric_limits<INDEX>::max()) {
         vector<REAL> mc(min_conv.size());
         vector<INDEX> mc_a(min_conv.size());;
         tropical_convolution::min_conv(left.begin(), left.end() ,right.begin(), right.end(), mc.begin(), mc.end(), mc_a.begin());

         REAL val = std::numeric_limits<REAL>::infinity();
         for(INDEX i=0; i<min_conv.size(); ++i) {
            const REAL cur_val = normalize(mc[i] - min_conv[i]);
            if(cur_val <= val) {
               val = cur_val;
               primal_left = mc_a[i];
               assert(mc_a[i] <= i);
               primal_right = i - primal_left;
            }
         }
         primal_sum = primal_left + primal_right;
         assert(primal_valid());
         assert(EvaluatePrimal() < std::numeric_limits<REAL>::infinity());
      } else if(primal_sum != std::numeric_limits<INDEX>::max() && primal_left == std::numeric_limits<INDEX>::max() && primal_right == std::numeric_limits<INDEX>::max()) {
         REAL val;
         std::tie(val, primal_left, primal_right) = tropical_convolution::arg_min_sum(left.begin(), left.end(), right.begin(), right.end(), primal_sum); 
         assert(primal_valid());
      } else if(primal_sum == std::numeric_limits<INDEX>::max() && primal_left != std::numeric_limits<INDEX>::max() && primal_right != std::numeric_limits<INDEX>::max()) {
         primal_sum = primal_left + primal_right;
         assert(primal_valid());
         assert(EvaluatePrimal() < std::numeric_limits<REAL>::infinity());
      } else {
         assert(false);
         std::cout << "maximize potential with given partial information not implemented yet\n";
      } 
      assert(primal_valid());
   }

   INDEX size() const { return 1; }

   void init_primal() { primal_left = std::numeric_limits<INDEX>::max(); primal_right = std::numeric_limits<INDEX>::max(); primal_sum = std::numeric_limits<INDEX>::max(); }
   template<class ARCHIVE> void serialize_primal(ARCHIVE& ar) { ar( primal_left, primal_right, primal_sum ); }
   template<class ARCHIVE> void serialize_dual(ARCHIVE& ar) { ar( left, right ); }

   vector<REAL> left, right, min_conv; // those vectors are shared with factors left below and right below and above respectively.
   // min conv here is taken with negative sing, and becomes left or right potential in the factor above
   INDEX primal_left, primal_right, primal_sum;

   bool primal_valid() const 
   {
      if(primal_left >= left.size()) { return false; }
      if(primal_right >= right.size()) { return false; }
      if(primal_sum >= min_conv.size()) { return false; }
      if(primal_left + primal_right != primal_sum) { return false; }
      return true;
   } 

private:
};

// left is lower, right is upper cardinality factor
// possibly hard code chirality as template parameter
class discrete_tomography_cardinality_message {
public:
   discrete_tomography_cardinality_message(Chirality _c) : c(_c) {}

   template<typename LEFT_FACTOR, typename MSG>
   void RepamLeft(LEFT_FACTOR& f_left, const MSG& msg)
   {
      assert(f_left.min_conv.size() == msg.size());
      // nothing needs to be done, reparametrization is shared
   }

   template<typename RIGHT_FACTOR, typename MSG>
   void RepamRight(RIGHT_FACTOR& f_right, const MSG& msg)
   {
      if(c == Chirality::left) {
         assert(f_right.left.size() == msg.size());
         for(INDEX i=0; i<f_right.left.size(); ++i) {
            f_right.left[i] += normalize(msg[i]);
         }
      } else {
         assert(c == Chirality::right);
         assert(f_right.right.size() == msg.size());
         for(INDEX i=0; i<f_right.right.size(); ++i) {
            f_right.right[i] += normalize(msg[i]);
         }
      }
   }

   template<typename LEFT_FACTOR, typename MSG>
   void send_message_to_right(const LEFT_FACTOR& f_left, MSG& msg, const REAL omega)
   {
      auto& left = f_left.left;
      auto& right = f_left.right;
      const auto& min_conv = f_left.min_conv;
      vector<REAL> msg_val(min_conv.size());
      tropical_convolution::min_conv(left.begin(), left.end(), right.begin(), right.end(), msg_val.begin(), msg_val.end());
      assert(msg_val.size() == min_conv.size());

      for(INDEX i=0; i<min_conv.size(); ++i) {
         msg_val[i] = normalize(msg_val[i] - min_conv[i]);
      }

      msg -= omega*msg_val;
   }

   template<typename RIGHT_FACTOR, typename MSG>
   void send_message_to_left(const RIGHT_FACTOR& f_right, MSG& msg, const REAL omega)
   {
      assert(false);
   }

   template<typename LEFT_FACTOR, typename RIGHT_FACTOR>
   bool ComputeRightFromLeftPrimal(const LEFT_FACTOR& l, RIGHT_FACTOR& r)
   {
      const INDEX left_val = l.primal_sum; 
      INDEX& right_val = c == Chirality::left ? r.primal_left : r.primal_right;
      if(left_val < l.size() && left_val != right_val) {
         right_val = left_val;
         return true;
      } 
      return false;
   }

   template<typename LEFT_FACTOR, typename RIGHT_FACTOR>
   bool ComputeLeftFromRightPrimal(LEFT_FACTOR& l, const RIGHT_FACTOR& r)
   {
      const INDEX right_val = c == Chirality::left ? r.primal_left : r.primal_right;
      INDEX& left_val = l.primal_sum; 
      if(right_val != std::numeric_limits<INDEX>::max() && right_val != left_val) {
         left_val = right_val;
         return true;
      }
      return false;
   }

   template<typename LEFT_FACTOR, typename RIGHT_FACTOR>
   bool CheckPrimalConsistency(const LEFT_FACTOR& l, const RIGHT_FACTOR& r) const
   {
      if(c == Chirality::left) {
         return l.primal_sum == r.primal_left;
      } else {
         assert(c == Chirality::right);
         return l.primal_sum == r.primal_right;
      }
   }

private:
   Chirality c;
};

// left is unary simplex factor, right is cardinality factor
class unary_simplex_discrete_tomography_cardinality_message {
public:
   unary_simplex_discrete_tomography_cardinality_message(Chirality _c) : c(_c) {}

   template<typename LEFT_FACTOR, typename MSG>
   void RepamLeft(LEFT_FACTOR& f_left, const MSG& msg)
   {
      assert(f_left.size() == msg.size());
      for(INDEX i=0; i<f_left.size(); ++i) {
         f_left[i] += normalize(msg[i]);
      }
   }

   template<typename RIGHT_FACTOR, typename MSG>
   void RepamRight(RIGHT_FACTOR& f_right, const MSG& msg)
   {
      if(c == Chirality::left) {
         assert(msg.size() == f_right.left.size());
         for(INDEX i=0; i<f_right.left.size(); ++i) {
            f_right.left[i] += normalize(msg[i]);
         }
      } else {
         assert(c == Chirality::right);
         assert(msg.size() == f_right.right.size());
         for(INDEX i=0; i<f_right.right.size(); ++i) {
            f_right.right[i] += normalize(msg[i]);
         }
      }
   }

   template<typename LEFT_FACTOR, typename MSG>
   void send_message_to_right(const LEFT_FACTOR& f_left, MSG& msg, const REAL omega)
   {
      msg -= omega*f_left;
   }

   template<typename RIGHT_FACTOR, typename MSG>
   void send_message_to_left(const RIGHT_FACTOR& f_right, MSG& msg, const REAL omega)
   {
      assert(false);
   }

   template<typename LEFT_FACTOR, typename RIGHT_FACTOR>
   bool ComputeRightFromLeftPrimal(const LEFT_FACTOR& l, RIGHT_FACTOR& r)
   {
      const INDEX left_val = l.primal(); 
      INDEX& right_val = c == Chirality::left ? r.primal_left : r.primal_right;
      if(left_val < l.size() && left_val != right_val) {
         right_val = left_val;
         return true;
      } 
      return false;
   }

   template<typename LEFT_FACTOR, typename RIGHT_FACTOR>
   bool ComputeLeftFromRightPrimal(LEFT_FACTOR& l, const RIGHT_FACTOR& r)
   {
      INDEX& left_val = l.primal(); 
      const INDEX right_val = c == Chirality::left ? r.primal_left : r.primal_right;
      if(right_val != std::numeric_limits<INDEX>::max() && right_val != left_val) {
         left_val = right_val;
         return true;
      }
      return false;
   }


   template<typename LEFT_FACTOR, typename RIGHT_FACTOR>
   bool CheckPrimalConsistency(const LEFT_FACTOR& l, const RIGHT_FACTOR& r) const
   {
      if(c == Chirality::left) {
         return l.primal() == r.primal_left;
      } else {
         assert(c == Chirality::right);
         return l.primal() == r.primal_right;
      }
   }
private:
   Chirality c;
};

} // end namespace LP_MP

#endif // LP_MP_DISCRETE_TOMOGRAPHY_CARDINALITY_HXX
