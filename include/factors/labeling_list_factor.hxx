#ifndef LP_MP_LABELING_LIST_FACTOR_HXX
#define LP_MP_LABELING_LIST_FACTOR_HXX

#include <array>
#include <bitset>
#include "vector.hxx"
#include "config.hxx"
//#include "cereal/types/bitset.hpp"
//#include "cereal/archives/binary.hpp"

#ifdef WITH_SAT
#include "sat_interface.hxx"
#endif

// to do: make _impl functions private
// make functions static whenever they can be

namespace LP_MP {

template<INDEX... LABELS>
struct labeling {

   template<INDEX LABEL_NO, INDEX LABEL, INDEX... LABELS_REST>
   constexpr static
   typename std::enable_if<LABEL_NO == 0,INDEX>::type label_impl()
   {
     static_assert(LABEL == 0 || LABEL == 1,"");
      return LABEL;
   }

   template<INDEX LABEL_NO, INDEX LABEL, INDEX... LABELS_REST>
   constexpr static
   typename std::enable_if<(LABEL_NO > 0),INDEX>::type label_impl()
   {
      return label_impl<LABEL_NO-1, LABELS_REST...>();
   }

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

   static bool matches(const std::bitset< sizeof...(LABELS)>& l)
   {
      return matches_impl<0, LABELS...>(l); 
   }

   template<INDEX I, INDEX... LABELS_REST>
   static typename std::enable_if<(I == sizeof...(LABELS))>::type
   get_labeling_impl(const std::bitset< sizeof...(LABELS)>& l)
   {
      return;
   }

   template<INDEX I, INDEX LABEL, INDEX... LABELS_REST>
   static typename std::enable_if<(I < sizeof...(LABELS))>::type
   get_labeling_impl(std::bitset< sizeof...(LABELS)>& l)
   {
     l[I] = (LABEL == 0) ? false : true;
     return get_labeling_impl<I+1, LABELS_REST...>(l);
   }

   static std::bitset<no_labels()> get_labeling()
   {
     std::bitset<no_labels()> l;
     get_labeling_impl<0, LABELS...>(l);
     return l;
   }
};

template<typename... LABELINGS> // all labels must be instances of labeling
struct labelings
{
   /* destructor makes labelings a non-constant type
   ~labelings()
   {
      static_assert(sizeof...(LABELINGS) > 0, "at least one labeling must be present");
      // to do: check whether each label occurs at most once and each labeling has same number of labels.
   }
   */

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
      static_assert(sizeof...(LABELINGS) > 0, "at least one labeling must be present");
      // to do: check whether each label occurs at most once and each labeling has same number of labels.

      static_assert(LABELING_NO < sizeof...(LABELINGS), "labeling number must be smaller than number of labelings");
      return get_label<LABELING_NO,LABEL_NO,LABELINGS...>();
   }

   constexpr static INDEX no_labelings()
   {
      return sizeof...(LABELINGS);
   }

   template<typename LABELING, typename... LABELINGS_REST>
   constexpr static INDEX no_labels_impl()
   {
      return LABELING::no_labels();
   }
   constexpr static INDEX no_labels()
   {
      return no_labels_impl<LABELINGS...>();
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

   static INDEX matching_labeling(const std::bitset< no_labels()>& l)
   {
      return matching_labeling_impl<0, LABELINGS...>(l); 
   }

   template<INDEX I, typename... LABELINGS_REST>
   static typename std::enable_if<(I >= sizeof...(LABELINGS)), std::bitset<no_labels()>>::type
   labeling_impl(const INDEX no)
   {
     assert(false);
      return std::bitset<no_labels()>();
   }

   template<INDEX I, typename LABELING, typename... LABELINGS_REST>
   static typename std::enable_if<(I < sizeof...(LABELINGS)), std::bitset<no_labels()>>::type
   labeling_impl(const INDEX no)
   {
     if(I == no) {
       return LABELING::get_labeling();
     } else {
       return labeling_impl<I+1,LABELINGS_REST...>(no);
     }
   }

   static std::bitset<no_labels()> labeling(const INDEX no)
   {
     return labeling_impl<0, LABELINGS...>(no);
   }

};

template<>
class labelings<>
{
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

template<typename LABELINGS, bool IMPLICIT_ORIGIN>
class labeling_factor : public array<REAL, LABELINGS::no_labelings()>
{
public:
   labeling_factor() 
   {
      std::fill(this->begin(), this->end(), 0.0);
   }
   ~labeling_factor()
   {}

   constexpr static bool has_implicit_origin() { return IMPLICIT_ORIGIN; } // means zero label has cost 0 and is not recorded.

   constexpr static INDEX size() {
      return LABELINGS::no_labelings();
   }
   
   constexpr static INDEX primal_size() {
      return LABELINGS::no_labels(); 
   }

   REAL LowerBound() const
   {

      /*
      static_assert(std::is_same<REAL, double>::value, "");
         // for this to work, array must be aligned
      REAL min;
      auto it = this->begin();
      // first compute minimum entry with SIMD
      if(this->size() >= 4) {
         simdpp::float64<4> min_vec = simdpp::load(it);
         for(it+=4; it+4<this->end(); it+=4) {
            simdpp::float64<4> tmp = simdpp::load( it );
            min_vec = simdpp::min(min_vec, tmp);
         }

         min = simdpp::reduce_min(min_vec);
      } else {
         min = *it;
         ++it;
      }
      for(; it<this->end(); ++it) {
         min = std::min(min, *it);
      } 
      */

      if(has_implicit_origin()) {
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

   // return two possible variable states
   // branch on current primal vs. not current primal
   void branch_left()
   {
      // set cost of label associated not with primal to infinity, i.e. current label should always be taken
      const INDEX labeling_no = LABELINGS::matching_labeling(primal_);
      assert(labeling_no < this->size());
      for(INDEX i=0; i<this->size(); ++i) {
         (*this)[i] = std::numeric_limits<REAL>::infinity(); 
      }
      // also the zero labeling must be forbidden. How to do? IMPLICIT_ORIGIN must be dropped for this.
      assert(false);
      assert(!has_implicit_origin());
   }

   void branch_right()
   {
      // set cost of primal label to infinity
      assert(EvaluatePrimal() < std::numeric_limits<REAL>::infinity());
      const INDEX labeling_no = LABELINGS::matching_labeling(primal_);
      assert(labeling_no < this->size());
      (*this)[labeling_no] = std::numeric_limits<REAL>::infinity();
      assert(LowerBound() < std::numeric_limits<REAL>::infinity());
   }

   auto& primal() { return primal_; }
   const auto& primal() const { return primal_; }

   void init_primal() {}
   template<typename ARCHIVE> void serialize_dual(ARCHIVE& ar) { ar( binary_data<REAL>(&(*this)[0], size()) ); }//*static_cast<array<REAL,size()>*>(this) ); }
   template<typename ARCHIVE> void serialize_primal(ARCHIVE& ar) { ar( primal_ ); }

#ifdef WITH_SAT
   template<typename SAT_SOLVER>
   void construct_sat_clauses(SAT_SOLVER& s) const
   {
      const INDEX no_vars = size() + (has_implicit_origin() ? 1 : 0);
      auto vars = create_sat_variables(s, no_vars);
      add_simplex_constraint_sat(s, vars.begin(), vars.end());
   }

   template<typename VEC>
   void reduce_sat(VEC& assumptions, const REAL th, sat_var begin) const
   {
      const REAL lb = LowerBound();
      for(INDEX i=0; i<this->size(); ++i) {
         if((*this)[i] > lb + th) { 
            assumptions.push_back(-to_literal(begin+i));
         }
      } 
   }

   template<typename SAT_SOLVER>
   void convert_primal(SAT_SOLVER& s, sat_var first)
   {
      primal_.reset();
      for(INDEX i=first; i<first+this->size(); ++i) {
         if(lglderef(s,to_literal(i)) == 1) {
           primal_ = LABELINGS::labeling(i); 
         }
      }
   }
#endif

private:
   std::bitset<primal_size()> primal_;
};

// we assume that LEFT_LABELING contains sublageings of RIGHT_LABELING, where we INDICES indicate i-th entry of LEFT_LABELING is mapped to INDICES[i]-th entry of right labeling
template<typename LEFT_LABELINGS, typename RIGHT_LABELINGS, INDEX... INDICES>
class labeling_message {
using type = labeling_message<LEFT_LABELINGS, RIGHT_LABELINGS, INDICES...>;
using msg_val_type = array<REAL, LEFT_LABELINGS::no_labelings()>;

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

  template<typename RIGHT_FACTOR, INDEX I, typename... RIGHT_LABELINGS_REST>
  typename std::enable_if<(I >= RIGHT_LABELINGS::no_labelings())>::type 
  compute_msg_impl(msg_val_type& msg_val, const RIGHT_FACTOR& r, REAL& min_of_labels_not_taken, labelings<RIGHT_LABELINGS_REST...>)
  {
     return;
  }

  template<typename RIGHT_FACTOR, INDEX I, typename RIGHT_LABELING, typename... RIGHT_LABELINGS_REST>
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
     compute_msg_impl<RIGHT_FACTOR, I+1, RIGHT_LABELINGS_REST...>(msg_val, r, min_of_labels_not_taken, labelings<RIGHT_LABELINGS_REST...>{});
  }

  template<typename RIGHT_FACTOR>
  void compute_msg(msg_val_type& msg_val, const RIGHT_FACTOR& r)
  {
     REAL min_of_labels_not_taken;
     if(r.has_implicit_origin()) {
        min_of_labels_not_taken = 0.0;
     } else {
        min_of_labels_not_taken = std::numeric_limits<REAL>::infinity();
     }
     std::fill(msg_val.begin(), msg_val.end(), std::numeric_limits<REAL>::infinity());
     compute_msg_impl<RIGHT_FACTOR, 0>(msg_val, r, min_of_labels_not_taken, RIGHT_LABELINGS{});
     // note: this is possibly wrong, if r.has_implicit_origin() is false
     for(auto& v : msg_val) {
        v -= min_of_labels_not_taken;
     }
     for(INDEX i=0; i<msg_val.size(); ++i) {
        assert(!std::isnan(msg_val[i]));
     }
  }

   template<typename RIGHT_FACTOR, typename MSG, INDEX I, typename... RIGHT_LABELINGS_REST>
   typename std::enable_if<(I >= RIGHT_LABELINGS::no_labelings())>::type 
   repam_right_impl(RIGHT_FACTOR& r, const MSG& msg, labelings<RIGHT_LABELINGS_REST...>)
   {
      return;
   }

   template<typename RIGHT_FACTOR, typename MSG, INDEX I, typename RIGHT_LABELING, typename... RIGHT_LABELINGS_REST>
   typename std::enable_if<(I < RIGHT_LABELINGS::no_labelings())>::type 
   repam_right_impl(RIGHT_FACTOR& r, const MSG& msg, labelings<RIGHT_LABELING, RIGHT_LABELINGS_REST...>)
   {
     INDEX left_label_number = matching_left_labeling<RIGHT_LABELING>(); // note: we should be able to qualify with constexpr! Is this an llvm bug?
     if(left_label_number < msg.size()) {
        r[I] += msg[left_label_number];
     }
     repam_right_impl<RIGHT_FACTOR, MSG, I+1, RIGHT_LABELINGS_REST...>(r, msg, labelings<RIGHT_LABELINGS_REST...>{});
   }

   template<typename RIGHT_FACTOR, typename MSG>
   void RepamRight(RIGHT_FACTOR& r, const MSG& msg)
   {
      // msg has dimension equal to number of left labelings;
      // go over all right labelings, find corresponding left labeling (if there is any) and if so, add msg value
      for(INDEX i=0; i<LEFT_LABELINGS::no_labelings(); ++i) {
         assert(!std::isnan(msg[i]));
         assert(std::isfinite(msg[i]));
      } 
      repam_right_impl<RIGHT_FACTOR, MSG, 0>(r, msg, RIGHT_LABELINGS{}); 
      for(INDEX i=0; i<r.size(); ++i) {
         assert(!std::isnan(r[i]));
         assert(std::isfinite(r[i]));
      }
   }

   template<typename RIGHT_FACTOR, typename MSG>
   void send_message_to_left(const RIGHT_FACTOR& r, MSG& msg, const REAL omega)
   {
      // msg has dimension equal to number of left labelings;
      // go over all right labelings. Then find left labeling corresponding to it, and compute minimum
      msg_val_type msg_val;
      compute_msg(msg_val, r);
      msg -= omega*msg_val;
   }

   template<typename LEFT_FACTOR, typename MSG>
   void RepamLeft(LEFT_FACTOR& l, const MSG& msg)
   {
      for(INDEX i=0; i<LEFT_LABELINGS::no_labelings(); ++i) {
         assert(std::isfinite(msg[i]));
         assert(!std::isnan(msg[i]));
         l[i] += msg[i];
         assert(!std::isnan(l[i]));
         assert(std::isfinite(l[i]));
      }
   }

   template<typename LEFT_FACTOR, typename MSG>
   void send_message_to_right(const LEFT_FACTOR& l, MSG& msg, const REAL omega)
   {
      for(INDEX i=0; i<l.size(); ++i) { assert(!std::isnan(l[i])); }
      msg -= omega*l;
   }

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

   template<typename LEFT_FACTOR, typename RIGHT_FACTOR>
   void ComputeRightFromLeftPrimal(const LEFT_FACTOR& l, RIGHT_FACTOR& r)
   {
      compute_right_from_left_primal_impl<0, INDICES...>(l.primal(), r.primal()); 
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

#ifdef WITH_SAT
   template<INDEX LEFT_LABELING_NO, typename... RIGHT_LABELINGS_REST>
   constexpr static std::size_t no_corresponding_labelings_impl(labelings<RIGHT_LABELINGS_REST...>)
   {
      return 0;
   }
   template<INDEX LEFT_LABELING_NO, typename RIGHT_LABELING, typename... RIGHT_LABELINGS_REST>
   constexpr static std::size_t no_corresponding_labelings_impl(labelings<RIGHT_LABELING, RIGHT_LABELINGS_REST...>)
   {
      if(matching_left_labeling<RIGHT_LABELING>() == LEFT_LABELING_NO) {
         return 1 + no_corresponding_labelings_impl<LEFT_LABELING_NO>(labelings<RIGHT_LABELINGS_REST...>{});
      } else {
         return no_corresponding_labelings_impl<LEFT_LABELING_NO>(labelings<RIGHT_LABELINGS_REST...>{});
      }
   }
   template<INDEX LEFT_LABELING_NO>
   constexpr static std::size_t no_corresponding_labelings()
   {
      return no_corresponding_labelings_impl<LEFT_LABELING_NO>(RIGHT_LABELINGS{}); 
   }


   template<INDEX LEFT_LABELING_NO, INDEX RIGHT_LABELING_IDX, typename ARRAY_IT, typename... RIGHT_LABELINGS_REST>
   void corresponding_labelings_impl(ARRAY_IT idx, labelings<RIGHT_LABELINGS_REST...>) const
   {
      static_assert(RIGHT_LABELING_IDX == RIGHT_LABELINGS::no_labelings(), "");
      return;
   }
   template<INDEX LEFT_LABELING_NO, INDEX RIGHT_LABELING_IDX, typename ARRAY_IT, typename RIGHT_LABELING, typename... RIGHT_LABELINGS_REST>
   void corresponding_labelings_impl(ARRAY_IT it, labelings<RIGHT_LABELING, RIGHT_LABELINGS_REST...>) const
   {
      INDEX left_label_number = matching_left_labeling<RIGHT_LABELING>(); // note: we should be able to qualify with constexpr!
      if(left_label_number == LEFT_LABELING_NO) {
         *it = RIGHT_LABELING_IDX;
         corresponding_labelings_impl<LEFT_LABELING_NO, RIGHT_LABELING_IDX+1>(it+1, labelings<RIGHT_LABELINGS_REST...>{});
      } else {
         corresponding_labelings_impl<LEFT_LABELING_NO, RIGHT_LABELING_IDX+1>(it, labelings<RIGHT_LABELINGS_REST...>{});
      } 
   }
   template<INDEX LEFT_LABELING_NO>
   std::array< sat_var, no_corresponding_labelings<LEFT_LABELING_NO>() > corresponding_labelings() const
   {
      std::array< sat_var, no_corresponding_labelings<LEFT_LABELING_NO>() > idx;
      corresponding_labelings_impl<LEFT_LABELING_NO, 0>(idx.begin(), RIGHT_LABELINGS{});
      return idx;
   }


   template<INDEX LEFT_LABELING_NO, typename SAT_SOLVER>
   typename std::enable_if<(LEFT_LABELING_NO >= LEFT_LABELINGS::no_labelings())>::type
   construct_sat_clauses_impl(SAT_SOLVER& s, const sat_var left_begin, const sat_var right_begin) const
   {}
   template<INDEX LEFT_LABELING_NO, typename SAT_SOLVER>
   typename std::enable_if<(LEFT_LABELING_NO < LEFT_LABELINGS::no_labelings())>::type
   construct_sat_clauses_impl(SAT_SOLVER& s, const sat_var left_begin, const sat_var right_begin) const
   {
      auto right_idx = corresponding_labelings<LEFT_LABELING_NO>();
      for(auto& idx : right_idx) {
         idx += right_begin;
      }
      auto one_active = one_active_indicator_sat(s, right_idx.begin(), right_idx.end());
      make_sat_var_equal(s, to_literal(left_begin + LEFT_LABELING_NO), to_literal(one_active));
      construct_sat_clauses_impl<LEFT_LABELING_NO+1>(s, left_begin, right_begin);
   }
   template<typename SAT_SOLVER, typename LEFT_FACTOR, typename RIGHT_FACTOR>
   void construct_sat_clauses(SAT_SOLVER& s, const LEFT_FACTOR& l, const RIGHT_FACTOR& r, const sat_var left_begin, const sat_var right_begin) const
   {
      construct_sat_clauses_impl<0>(s, left_begin, right_begin);
   }

#endif
private:
};

} // end namespace LP_MP

#endif //  LP_MP_LABELING_LIST_FACTOR_HXX
