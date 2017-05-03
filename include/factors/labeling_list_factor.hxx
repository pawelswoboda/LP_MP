#ifndef LP_MP_LABELING_LIST_FACTOR_HXX
#define LP_MP_LABELING_LIST_FACTOR_HXX

#include <array>
#include <bitset>
#include "vector.hxx"
#include "config.hxx"

// factors and messages for simplex over compile-time specified labelings.
// to do: make _impl functions private
// to do: possilby code can be made nicer with boost::hana

namespace LP_MP {

template<INDEX... LABELS>
class labeling {
public:

   template<INDEX LABEL_NO>
   constexpr static INDEX label()
   { 
      static_assert(LABEL_NO < sizeof...(LABELS), "label number must be smaller than number of labels");
      return label_impl<LABEL_NO, LABELS...>();
   } 

   constexpr static INDEX no_labels()
   {
      return sizeof...(LABELS);
   }

   static bool matches(const std::bitset< sizeof...(LABELS)>& l)
   {
      return matches_impl<0, LABELS...>(l); 
   }

   static bool partial_matches(const std::bitset< sizeof...(LABELS)>& l, const std::bitset<sizeof...(LABELS)>& mask)
   {
      return partial_matches_impl<0, LABELS...>(l, mask);
   }

private:

   template<INDEX LABEL_NO, INDEX LABEL, INDEX... LABELS_REST>
   constexpr static
   typename std::enable_if<LABEL_NO == 0,INDEX>::type label_impl()
   {
      return LABEL;
   }

   template<INDEX LABEL_NO, INDEX LABEL, INDEX... LABELS_REST>
   constexpr static
   typename std::enable_if<(LABEL_NO > 0),INDEX>::type label_impl()
   {
      return label_impl<LABEL_NO-1, LABELS_REST...>();
   }

   template<INDEX I, INDEX... LABELS_REST>
   static typename std::enable_if<(I == sizeof...(LABELS)),bool>::type
   matches_impl(const std::bitset< sizeof...(LABELS)>& l)
   {
      return true;
   }

   template<INDEX I, INDEX LABEL, INDEX... LABELS_REST>
   static typename std::enable_if<(I < sizeof...(LABELS)),bool>::type
   matches_impl(const std::bitset< sizeof...(LABELS)>& l)
   {
      if(l[I] != LABEL) {
         return false;
      } else {
         return matches_impl<I+1, LABELS_REST...>(l);
      }
   } 

   template<INDEX I, INDEX... LABELS_REST>
   static typename std::enable_if<(I == sizeof...(LABELS)),bool>::type
   partial_matches_impl(const std::bitset< sizeof...(LABELS)>& l, const std::bitset< sizeof...(LABELS)>& mask)
   {
      return true;
   }

   template<INDEX I, INDEX LABEL, INDEX... LABELS_REST>
   static typename std::enable_if<(I < sizeof...(LABELS)),bool>::type
   partial_matches_impl(const std::bitset< sizeof...(LABELS)>& l, const std::bitset< sizeof...(LABELS)>& mask)
   {
      if((mask[I] && l[I] == LABEL) || !mask[I]) {
         return partial_matches_impl<I+1, LABELS_REST...>(l, mask);
      } else {
         assert(mask[I] && l[I] != LABEL);
         return false;
      }
   } 

};

template<typename... LABELINGS> // all labels must be instances of labeling
class labelings
{
public:
   ~labelings()
   {
      static_assert(sizeof...(LABELINGS) > 0, "at least one labeling must be present");
      // to do: check whether each label occurs at most once and each labeling has same number of labels.
   }

   template<INDEX LABELING_NO, INDEX LABEL_NO, typename LABELING, typename... LABELINGS_REST>
   constexpr static typename std::enable_if<LABELING_NO == 0,INDEX>::type 
   get_label()
   {
      return LABELING::template label<LABEL_NO>();
   }

   template<INDEX LABELING_NO, INDEX LABEL_NO, typename LABELING, typename... LABELINGS_REST>
   constexpr static typename std::enable_if<(LABELING_NO > 0),INDEX>::type 
   get_label()
   {
      return get_label<LABELING_NO-1, LABEL_NO, LABELINGS_REST...>();
   } 

   template<INDEX LABELING_NO, INDEX LABEL_NO>
   constexpr static INDEX label()
   {
      static_assert(LABELING_NO < sizeof...(LABELINGS), "labeling number must be smaller than number of labelings");
      return get_label<LABELING_NO,LABEL_NO,LABELINGS...>();
   } 

   constexpr static INDEX no_labels()
   {
      return no_labels_impl<LABELINGS...>();
   }

   template<typename LABELING, typename... LABELINGS_REST>
   constexpr static INDEX no_labels_impl()
   {
      return LABELING::no_labels();
   }

   constexpr static INDEX no_labelings()
   {
      return sizeof...(LABELINGS);
   }

   template<INDEX I, typename... LABELINGS_REST>
   static typename std::enable_if<(I >= sizeof...(LABELINGS)), INDEX>::type
   matching_labeling_impl(const std::bitset< no_labels()>& l)
   {
      return no_labelings();
   }

   template<INDEX I, typename LABELING, typename... LABELINGS_REST>
   static typename std::enable_if<(I < sizeof...(LABELINGS)), INDEX>::type
   matching_labeling_impl(const std::bitset< no_labels()>& l)
   {
      if(LABELING::matches(l)) {
         return I;
      } else {
         return matching_labeling_impl<I+1,LABELINGS_REST...>(l);
      }
   } 

   template<INDEX LABELING_NO, INDEX I, typename... LABELINGS_REST>
   static typename std::enable_if<(I >= sizeof...(LABELINGS)), INDEX>::type
   partial_match_impl(const std::bitset< no_labels()>& l, const std::bitset<no_labels()>& mask)
   {
      return false;
   }

   template<INDEX LABELING_NO, INDEX I, typename LABELING, typename... LABELINGS_REST>
   static typename std::enable_if<(I < sizeof...(LABELINGS)), INDEX>::type
   partial_match_impl(const std::bitset< no_labels()>& l, const std::bitset<no_labels()>& mask)
   {
      if(LABELING_NO == I) {
         return LABELING::partial_matches(l, mask);
      } else {
         return partial_match_impl<LABELING_NO,I+1,LABELINGS_REST...>(l, mask);
      }
   } 
   
   static INDEX matching_labeling(const std::bitset<no_labels()>& l)
   {
      return matching_labeling_impl<0, LABELINGS...>(l); 
   }

   template<INDEX LABELING_NO>
   static bool partial_match(const std::bitset<no_labels()>& l, const std::bitset<no_labels()>& mask)
   {
      return partial_match_impl<LABELING_NO, 0, LABELINGS...>(l, mask); 
   }

};

template<>
class labelings<>
{
public:
   template<INDEX LABELING_NO, INDEX LABEL_NO>
   constexpr static INDEX label()
   {
      return 42;
   }

   constexpr static INDEX no_labelings()
   {
      return 0;
   }

   constexpr static INDEX no_labels()
   {
      return 0;
   }

   template<typename VEC>
   static INDEX matching_labeling(const VEC& l)
   {
      return 0;
   }
};

// MP_ROUNDING = true enables additional unknown label and informs labeling_message to also provide ReceiveRestrictedMessageFromRight
template<typename LABELINGS, bool IMPLICIT_ORIGIN, bool MP_ROUNDING = false>
class labeling_factor : public array<REAL, LABELINGS::no_labelings()>
{
public:
   using labelings = LABELINGS;

   constexpr static bool supports_mp_rounding() { return MP_ROUNDING; }

   labeling_factor() 
   {
      std::fill(this->begin(), this->end(), 0.0);
   }
   ~labeling_factor()
   {
      static_assert(IMPLICIT_ORIGIN, "implicit zero label obligatory yet");
   }

   constexpr static bool has_implicit_origin() { return IMPLICIT_ORIGIN; } // means zero label has cost 0 and is not recorded.

   constexpr static INDEX size() {
      return LABELINGS::no_labelings();
   }
   
   constexpr static INDEX primal_size() {
      return LABELINGS::no_labels(); 
   }

   REAL LowerBound() const
   {
      if(IMPLICIT_ORIGIN) {
         return std::min(0.0, *std::min_element(this->begin(), this->end()));
      } else {
         return *std::min_element(this->begin(), this->end());
      }
   }

   REAL EvaluatePrimal() const
   {
      const INDEX labeling_no = LABELINGS::matching_labeling(primal_);
      if(labeling_no < size()) {
         return (*this)[labeling_no];
      }
      // check for zero labeling
      if(has_implicit_origin() && primal_.count() == 0) {
         return 0.0;
      }
      return std::numeric_limits<REAL>::infinity();
   }

   template<INDEX LABELING_NO>
   bool partial_match() const
   {
      return LABELINGS::template partial_match<LABELING_NO>(primal_, labelled_mask_);
   }

   auto& primal() { return primal_; }
   const auto& primal() const { return primal_; }

   auto& mask() { return labelled_mask_; }
   const auto& mask() const { return labelled_mask_; }

   void init_primal() 
   {
      //if(MP_ROUNDING) {
         labelled_mask_.reset();
         for(INDEX i=0; i<labelled_mask_.size(); ++i) {
            assert(labelled_mask_[i] == 0);
         }
      //} 
   }
   template<typename ARCHIVE> void serialize_dual(ARCHIVE& ar) { ar( *static_cast<array<REAL,size()>*>(this) ); }
   template<typename ARCHIVE> void serialize_primal(ARCHIVE& ar) { ar( primal_ ); }

private:
   // possible use std::bitset for more space efficiency
   // for MP_ROUNDING use two bitsets, to encode unknown values as well
   std::bitset<primal_size()> primal_;
   std::bitset<primal_size()> labelled_mask_;
};

// this is not so nice: We only want to make std::array behave sanely with vector expression like - and multiplying it with a constant
template<typename T, INDEX N>
class array_expression : public std::array<T,N>, public vector_expression<T, array_expression<T,N>>
{
public:
   INDEX size() const { return std::array<T,N>::size(); } 
   T& operator[](const INDEX i) { return static_cast<std::array<T,N>*>(this)->operator[](i); }
   T operator[](const INDEX i) const { return static_cast<const std::array<T,N>*>(this)->operator[](i); }
};

// we assume that LEFT_LABELING contains sublageings of RIGHT_LABELING, where we INDICES indicate i-th entry of LEFT_LABELING is mapped to INDICES[i]-th entry of right labeling
//template<typename LEFT_LABELINGS, typename RIGHT_LABELINGS, INDEX... INDICES>
template<typename LEFT_FACTOR, typename RIGHT_FACTOR, INDEX... INDICES>
class labeling_message {
using type = labeling_message<LEFT_FACTOR, RIGHT_FACTOR, INDICES...>;
using LEFT_LABELINGS = typename LEFT_FACTOR::labelings;
using RIGHT_LABELINGS = typename RIGHT_FACTOR::labelings;
using msg_val_type = array_expression<REAL, LEFT_LABELINGS::no_labelings()>;

public:
   ~labeling_message() 
   {
      static_assert(sizeof...(INDICES) == LEFT_LABELINGS::no_labels(), "each left label must be matched to a right one");
   }

   template<typename LEFT_LABELING, typename RIGHT_LABELING, INDEX LEFT_INDEX, INDEX... I_REST>
   constexpr static typename std::enable_if<(LEFT_INDEX >= LEFT_LABELING::no_labels()),bool>::type
   matches_impl()
   {
      return true;
   }
   template<typename LEFT_LABELING, typename RIGHT_LABELING, INDEX LEFT_INDEX, INDEX I, INDEX... I_REST>
   constexpr static typename std::enable_if<(LEFT_INDEX < LEFT_LABELING::no_labels()),bool>::type
   matches_impl()
   {
      if(LEFT_LABELING::template label<LEFT_INDEX>() == RIGHT_LABELING::template label<I>()) {
         return type::matches_impl<LEFT_LABELING, RIGHT_LABELING, LEFT_INDEX+1, I_REST...>();
      } else {
         return false;
      }
   }

   template<typename LEFT_LABELING, typename RIGHT_LABELING>
   constexpr static bool matches()
   {
      return matches_impl<LEFT_LABELING, RIGHT_LABELING, 0, INDICES...>();
   }

   template<typename RIGHT_LABELING, INDEX I, typename... LEFT_LABELINGS_REST>
   constexpr static typename std::enable_if<(I >= LEFT_LABELINGS::no_labelings()),INDEX >::type
   matching_left_labeling_impl(labelings<LEFT_LABELINGS_REST...>)
   {
      return I; // return 1 + number of left labelings;
   }
   template<typename RIGHT_LABELING, INDEX I, typename LEFT_LABELING, typename... LEFT_LABELINGS_REST>
   constexpr static typename std::enable_if<(I < LEFT_LABELINGS::no_labelings()),INDEX>::type
   matching_left_labeling_impl(labelings<LEFT_LABELING, LEFT_LABELINGS_REST...>)
   {
      if(matches<LEFT_LABELING, RIGHT_LABELING>()) {
         return I;
      } else {
         return matching_left_labeling_impl<RIGHT_LABELING, I+1>(labelings<LEFT_LABELINGS_REST...>{});
      }
   }

   // we assume that there is as most one matching left labeling
   template<typename RIGHT_LABELING>
   constexpr static INDEX matching_left_labeling()
   {
      return matching_left_labeling_impl<RIGHT_LABELING,0>(LEFT_LABELINGS{});
   } 

  template<INDEX I, typename... RIGHT_LABELINGS_REST>
  typename std::enable_if<(I >= RIGHT_LABELINGS::no_labelings())>::type 
  compute_msg_impl(msg_val_type& msg_val, const RIGHT_FACTOR& r, REAL& min_of_labels_not_taken, labelings<RIGHT_LABELINGS_REST...>)
  {
     return;
  }

  template<INDEX I, typename RIGHT_LABELING, typename... RIGHT_LABELINGS_REST>
  typename std::enable_if<(I < RIGHT_LABELINGS::no_labelings())>::type 
  compute_msg_impl(msg_val_type& msg_val, const RIGHT_FACTOR& r, REAL& min_of_labels_not_taken, labelings<RIGHT_LABELING, RIGHT_LABELINGS_REST...>)
  {
     INDEX left_label_number = matching_left_labeling<RIGHT_LABELING>(); // note: we should be able to qualify with constexpr! Is this an llvm bug?
     if(left_label_number < msg_val.size()) {
        msg_val[left_label_number] = std::min(msg_val[left_label_number], r[I]);
     } else {
        assert(left_label_number == msg_val.size());
        min_of_labels_not_taken = std::min(min_of_labels_not_taken, r[I]); 
     }
     compute_msg_impl<I+1, RIGHT_LABELINGS_REST...>(msg_val, r, min_of_labels_not_taken, labelings<RIGHT_LABELINGS_REST...>{});
  }

  //template<typename RIGHT_FACTOR>
  void compute_msg(msg_val_type& msg_val, const RIGHT_FACTOR& r)
  {
     REAL min_of_labels_not_taken;
     if(r.has_implicit_origin()) {
        min_of_labels_not_taken = 0.0;
     } else {
        min_of_labels_not_taken = std::numeric_limits<REAL>::infinity();
     }
     std::fill(msg_val.begin(), msg_val.end(), std::numeric_limits<REAL>::infinity());
     compute_msg_impl<0>(msg_val, r, min_of_labels_not_taken, RIGHT_LABELINGS{});
     // note: this is possibly wrong, if r.has_implicit_origin() is false
     for(auto& v : msg_val) {
        v -= min_of_labels_not_taken;
     }
  }

  template<INDEX I, typename... RIGHT_LABELINGS_REST>
  typename std::enable_if<(I >= RIGHT_LABELINGS::no_labelings())>::type 
  compute_restricted_msg_impl(msg_val_type& msg_val, const RIGHT_FACTOR& r, REAL& min_of_labels_not_taken, labelings<RIGHT_LABELINGS_REST...>)
  {
     return;
  }

  template<INDEX I, typename RIGHT_LABELING, typename... RIGHT_LABELINGS_REST>
  typename std::enable_if<(I < RIGHT_LABELINGS::no_labelings())>::type 
  compute_restricted_msg_impl(msg_val_type& msg_val, const RIGHT_FACTOR& r, REAL& min_of_labels_not_taken, labelings<RIGHT_LABELING, RIGHT_LABELINGS_REST...>)
  {
     INDEX left_label_number = matching_left_labeling<RIGHT_LABELING>(); // note: we should be able to qualify with constexpr! Is this an llvm bug?
     const bool partial_match = r.template partial_match<I>();

     if(partial_match) {
        if(left_label_number < msg_val.size()) {
           msg_val[left_label_number] = std::min(msg_val[left_label_number], r[I]);
        } else {
           assert(left_label_number == msg_val.size());
           min_of_labels_not_taken = std::min(min_of_labels_not_taken, r[I]); 
        }
     }
     
     compute_restricted_msg_impl<I+1, RIGHT_LABELINGS_REST...>(msg_val, r, min_of_labels_not_taken, labelings<RIGHT_LABELINGS_REST...>{});
  }

  //template<typename RIGHT_FACTOR>
  void compute_restricted_msg(msg_val_type& msg_val, const RIGHT_FACTOR& r)
  {
     std::fill(msg_val.begin(), msg_val.end(), std::numeric_limits<REAL>::infinity());
     REAL min_of_labels_not_taken;
     if(r.has_implicit_origin()) {
        min_of_labels_not_taken = 0.0;
     } else {
        min_of_labels_not_taken = std::numeric_limits<REAL>::infinity();
     }

     compute_restricted_msg_impl<0>(msg_val, r, min_of_labels_not_taken, RIGHT_LABELINGS{});
     // note: this is possibly wrong is hass_implicit_origin is false()
     for(auto& v : msg_val) {
        v -= min_of_labels_not_taken;
     }
     for(auto v : msg_val) {
        assert(!std::isnan(v));
     }
  }

   template<typename MSG, INDEX I, typename... RIGHT_LABELINGS_REST>
   typename std::enable_if<(I >= RIGHT_LABELINGS::no_labelings())>::type 
   repam_right_impl(RIGHT_FACTOR& r, const MSG& msg, labelings<RIGHT_LABELINGS_REST...>)
   {
      return;
   }

   template<typename MSG, INDEX I, typename RIGHT_LABELING, typename... RIGHT_LABELINGS_REST>
   typename std::enable_if<(I < RIGHT_LABELINGS::no_labelings())>::type 
   repam_right_impl(RIGHT_FACTOR& r, const MSG& msg, labelings<RIGHT_LABELING, RIGHT_LABELINGS_REST...>)
   {
     INDEX left_label_number = matching_left_labeling<RIGHT_LABELING>(); // note: we should be able to qualify with constexpr! Is this an llvm bug?
     if(left_label_number < msg.size()) {
        r[I] += msg[left_label_number];
     }
     repam_right_impl<MSG, I+1, RIGHT_LABELINGS_REST...>(r, msg, labelings<RIGHT_LABELINGS_REST...>{});
   }

   template<typename MSG>
   void RepamRight(RIGHT_FACTOR& r, const MSG& msg)
   {
      // msg has dimension equal to number of left labelings;
      // go over all right labelings, find corresponding left labeling (if there is any) and if so, add msg value
      repam_right_impl<MSG, 0>(r, msg, RIGHT_LABELINGS{}); 
   }

   template<typename MSG>
   void ReceiveMessageFromRight(const RIGHT_FACTOR& r, MSG& msg)
   {
      // msg has dimension equal to number of left labelings;
      // go over all right labelings. Then find left labeling corresponding to it, and compute minimum
      msg_val_type msg_val; // additional last entry is minimum over unmatched right labelings
      compute_msg(msg_val, r);
      msg -= msg_val;
   }

   template<typename MSG, bool ENABLE = LEFT_FACTOR::supports_mp_rounding()>
   typename std::enable_if<ENABLE>::type 
   ReceiveRestrictedMessageFromRight(const RIGHT_FACTOR& r, MSG& msg)
   {
      msg_val_type msg_val;
      compute_restricted_msg(msg_val, r);
      msg -= msg_val;
   }

   template<typename MSG>
   void RepamLeft(LEFT_FACTOR& l, const MSG& msg)
   {
      for(INDEX i=0; i<LEFT_LABELINGS::no_labelings(); ++i) {
         l[i] += msg[i];
      }
   }

   template<typename MSG>
   void SendMessageToRight(const LEFT_FACTOR& l, MSG& msg, const REAL omega)
   {
      msg -= omega*l;
   }

   /*
   template<INDEX LEFT_INDEX, INDEX... RIGHT_INDICES_REST>
   typename std::enable_if<(LEFT_INDEX >= LEFT_LABELINGS::no_labels())>::type
   compute_right_from_left_primal_impl(const std::bitset<LEFT_LABELINGS::no_labels()>& l, std::bitset< RIGHT_LABELINGS::no_labels()>& r)
   {}
   template<INDEX LEFT_INDEX, INDEX RIGHT_INDEX, INDEX... RIGHT_INDICES_REST>
   typename std::enable_if<(LEFT_INDEX < LEFT_LABELINGS::no_labels())>::type
   compute_right_from_left_primal_impl(const std::bitset<LEFT_LABELINGS::no_labels()>& l, std::bitset< RIGHT_LABELINGS::no_labels()>& r)
   {
      r[RIGHT_INDEX] = l[LEFT_INDEX];
      compute_right_from_left_primal_impl<LEFT_INDEX+1, RIGHT_INDICES_REST...>(l,r); 
   }

   //template<typename LEFT_FACTOR, typename RIGHT_FACTOR>
   void ComputeRightFromLeftPrimal(const LEFT_FACTOR& l, RIGHT_FACTOR& r)
   {
      compute_right_from_left_primal_impl<0, INDICES...>(l.primal(), r.primal()); 
   } 
   */

   template<INDEX LEFT_INDEX, INDEX... RIGHT_INDICES_REST>
   typename std::enable_if<(LEFT_INDEX >= LEFT_LABELINGS::no_labels())>::type
   compute_right_from_left_primal_impl(
         const std::bitset<LEFT_LABELINGS::no_labels()>& l, std::bitset< RIGHT_LABELINGS::no_labels()>& r,
         const std::bitset<LEFT_LABELINGS::no_labels()>& l_mask, std::bitset< RIGHT_LABELINGS::no_labels()>& r_mask
         )
   {}
   template<INDEX LEFT_INDEX, INDEX RIGHT_INDEX, INDEX... RIGHT_INDICES_REST>
   typename std::enable_if<(LEFT_INDEX < LEFT_LABELINGS::no_labels())>::type
   compute_right_from_left_primal_impl(
         const std::bitset<LEFT_LABELINGS::no_labels()>& l, std::bitset< RIGHT_LABELINGS::no_labels()>& r,
         const std::bitset<LEFT_LABELINGS::no_labels()>& l_mask, std::bitset< RIGHT_LABELINGS::no_labels()>& r_mask
         )
   {
      r[RIGHT_INDEX] = l[LEFT_INDEX];
      r_mask[RIGHT_INDEX] = true;
      compute_right_from_left_primal_impl<LEFT_INDEX+1, RIGHT_INDICES_REST...>(l,r, l_mask, r_mask); 
   }

   //template<typename LEFT_FACTOR, typename RIGHT_FACTOR>
   void ComputeRightFromLeftPrimal(const LEFT_FACTOR& l, RIGHT_FACTOR& r)
   {
      compute_right_from_left_primal_impl<0, INDICES...>(l.primal(), r.primal(), l.mask(), r.mask()); 
   } 


   template<INDEX RIGHT_INDEX, typename... RIGHT_LABELINGS_REST>
   static typename std::enable_if<(RIGHT_INDEX >= RIGHT_LABELINGS::no_labelings())>::type
   print_matching_impl(labelings<RIGHT_LABELINGS_REST...>)
   {}
   template<INDEX RIGHT_INDEX, typename RIGHT_LABELING, typename... RIGHT_LABELINGS_REST>
   static typename std::enable_if<(RIGHT_INDEX < RIGHT_LABELINGS::no_labelings())>::type
   print_matching_impl(labelings<RIGHT_LABELING, RIGHT_LABELINGS_REST...>)
   {
      INDEX left_label_number = matching_left_labeling<RIGHT_LABELING>(); // note: we should be able to qualify with constexpr! Is this an llvm bug?
      std::cout << left_label_number << "," << RIGHT_INDEX << "\n";
      print_matching_impl<RIGHT_INDEX+1>(labelings<RIGHT_LABELINGS_REST...>{});
   }
   static void print_matching()
   {
      // for each right labeling, print left one that is matched (for debugging purposes)
      print_matching_impl<0>(RIGHT_LABELINGS{});
   }

private:
};

} // end namespace LP_MP

#endif //  LP_MP_LABELING_LIST_FACTOR_HXX
