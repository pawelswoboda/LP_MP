#ifndef LP_MP_MULTICUT_FACTORS_MESSAGES_HXX
#define LP_MP_MULTICUT_FACTORS_MESSAGES_HXX

#include "LP_MP.h"
#include "factors/labeling_list_factor.hxx"

namespace LP_MP {

using multicut_edge_labelings = labelings< labeling<1> >;
using multicut_edge_factor = labeling_factor< multicut_edge_labelings, true >;

using multicut_triplet_labelings = labelings<
   labeling<0,1,1>,
   labeling<1,0,1>,
   labeling<1,1,0>,
   labeling<1,1,1>
      >;
using multicut_triplet_factor = labeling_factor< multicut_triplet_labelings, true >;

using multicut_edge_triplet_message_0 = labeling_message< multicut_edge_labelings, multicut_triplet_labelings, 0 >;
using multicut_edge_triplet_message_1 = labeling_message< multicut_edge_labelings, multicut_triplet_labelings, 1 >;
using multicut_edge_triplet_message_2 = labeling_message< multicut_edge_labelings, multicut_triplet_labelings, 2 >;

using multicut_odd_3_wheel_labelings = labelings<
   // one node separated from rest
   labeling<1,1,1,0,0,0>,
   labeling<1,0,1,1,0,0>,
   labeling<0,1,1,0,0,1>,
   labeling<0,0,0,1,1,1>,

   // two components of two nodes each
   labeling<0,1,1,1,1,0>,
   labeling<1,0,1,1,0,1>,
   labeling<1,1,0,0,1,1>,

   // three components
   labeling<0,1,1,1,1,1>,
   labeling<1,0,1,1,1,1>,
   labeling<1,1,0,1,1,1>,
   labeling<1,1,1,0,1,1>,
   labeling<1,1,1,1,0,1>,
   labeling<1,1,1,1,1,0>,

   // four components
   labeling<1,1,1,1,1,1>
   >;
using multicut_odd_3_wheel_factor = labeling_factor< multicut_odd_3_wheel_labelings, true >;

using multicut_triplet_odd_3_wheel_message_012 = labeling_message< multicut_triplet_labelings, multicut_odd_3_wheel_labelings, 0,1,2>;
using multicut_triplet_odd_3_wheel_message_013 = labeling_message< multicut_triplet_labelings, multicut_odd_3_wheel_labelings, 0,3,4>;
using multicut_triplet_odd_3_wheel_message_023 = labeling_message< multicut_triplet_labelings, multicut_odd_3_wheel_labelings, 1,3,5>;
using multicut_triplet_odd_3_wheel_message_123 = labeling_message< multicut_triplet_labelings, multicut_odd_3_wheel_labelings, 2,4,5>;

/*
using multicut_triplet_plus_spoke_labelings = labelings<
   labeling<0,1,1,0>,
   labeling<1,0,1,0>,
   labeling<1,1,0,0>,
   labeling<1,1,1,0>,
   labeling<0,1,1,1>,
   labeling<1,0,1,1>,
   labeling<1,1,0,1>,
   labeling<1,1,1,1>,
   labeling<0,0,0,1>
      >;
using multicut_triplet_plus_spoke_factor = labeling_factor< multicut_triplet_plus_spoke_labelings, true >; 

using multicut_triplet_plus_spoke_cover_message = labeling_message< multicut_triplet_labelings, multicut_triplet_plus_spoke_labelings, 0,1,2 >;
using multicut_triplet_plus_spoke_message = labeling_message< multicut_triplet_labelings, multicut_triplet_plus_spoke_labelings, 0,1,2 >;
*/

} // end namespace LP_MP 

#endif // LP_MP_MULTICUT_FACTORS_MESSAGES_HXX
