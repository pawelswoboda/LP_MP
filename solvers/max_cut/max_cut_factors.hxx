#ifndef LP_MP_MAX_CUT_FACTORS
#define LP_MP_MAX_CUT_FACTORS

#include "LP_MP.h"
#include "factors/labeling_list_factor.hxx"

namespace LP_MP {

using max_cut_edge_labelings = labelings< labeling<1> >;
using max_cut_edge_factor = labeling_factor< max_cut_edge_labelings, true >;

using max_cut_triplet_labelings = labelings<
   labeling<1,1,0>, // 100
   labeling<1,0,1>, // 010
   labeling<0,1,1>  // 001
   >;
using max_cut_triplet_factor = labeling_factor< max_cut_triplet_labelings, true >;

using max_cut_edge_triplet_message_0 = labeling_message< max_cut_edge_labelings, max_cut_triplet_labelings, 0 >;
using max_cut_edge_triplet_message_1 = labeling_message< max_cut_edge_labelings, max_cut_triplet_labelings, 1 >;
using max_cut_edge_triplet_message_2 = labeling_message< max_cut_edge_labelings, max_cut_triplet_labelings, 2 >;

// do we also need factors on 4 nodes for triangulation purposes for odd bicycle wheels? I think so!
using max_cut_4_clique_labelings = labelings <
   labeling<1,1,0,1,0,0>, // 1000
   labeling<1,0,1,0,1,0>, // 0100
   labeling<0,1,1,0,0,1>, // 0010
   labeling<0,0,0,1,1,1>, // 0001

   labeling<0,1,1,1,1,0>, // 1100
   labeling<1,0,1,1,0,1>, // 1010
   labeling<1,1,0,0,1,1>  // 1001
      >;

using max_cut_4_clique_factor = labeling_factor< max_cut_4_clique_labelings, true >;

using max_cut_triplet_4_clique_message_012 = labeling_message< max_cut_triplet_labelings, max_cut_4_clique_labelings, 0,1,2 >;
using max_cut_triplet_4_clique_message_013 = labeling_message< max_cut_triplet_labelings, max_cut_4_clique_labelings, 0,3,4 >;
using max_cut_triplet_4_clique_message_023 = labeling_message< max_cut_triplet_labelings, max_cut_4_clique_labelings, 1,3,5 >;
using max_cut_triplet_4_clique_message_123 = labeling_message< max_cut_triplet_labelings, max_cut_4_clique_labelings, 2,4,5 >;

/*
                       i1



               
                    i3_
                       \__
                          i4




   i0                                       i2
   */
using max_cut_5_clique_labelings = labelings <
// one node, four nodes
labeling<1,1,0,1,1,0,0,0,0,0>, // 10000
labeling<1,0,1,0,0,1,1,0,0,0>, // 01000
labeling<0,1,1,0,0,0,0,1,1,0>, // 00100
labeling<0,0,0,1,0,1,0,1,0,1>, // 00010
labeling<0,0,0,0,1,0,1,0,1,1>, // 00001
// two nodes, three nodes
labeling<0,1,1,1,1,1,1,0,0,0>, // 11000
labeling<1,0,1,1,1,0,0,1,1,0>, // 10100
labeling<1,1,0,0,1,1,0,1,0,1>, // 10010
labeling<1,1,0,1,0,0,1,0,1,1>, // 10001
labeling<1,1,0,0,0,1,1,1,1,0>, // 01100
labeling<1,0,1,1,0,0,1,1,0,1>, // 01010
labeling<1,0,1,0,1,1,0,0,1,1>, // 01001
labeling<0,1,1,1,0,1,0,0,1,1>, // 00110
labeling<0,1,1,0,1,0,1,1,0,1>, // 00101
labeling<0,0,0,1,1,1,1,1,1,0>  // 00011
/*
   labeling<0,1,1,0,0,0,0,1,1,0>, // 00100
   labeling<1,0,1,1,1,0,0,0,0,0>, // 10000 
   labeling<1,1,0,0,0,1,1,0,0,0>, // 01000 

   labeling<0,1,1,1,0,1,0,0,1,1>, // 00110 
   labeling<1,0,1,0,1,1,0,1,0,1>, // 10010 
   labeling<1,1,0,1,0,0,1,1,0,1>, // 01010 

   labeling<0,1,1,1,0,1,0,0,1,1>, // 00101 
   labeling<1,0,1,1,0,0,1,0,1,1>, // 10001 
   labeling<1,1,0,0,1,1,0,0,1,1>, // 01001 

   labeling<0,1,1,1,1,1,1,0,0,0>, // 11000 
   labeling<1,0,1,0,0,1,1,1,1,0>, // 01100 
   labeling<1,1,0,1,1,0,0,1,1,0>, // 10100 

   labeling<0,0,0,1,0,1,0,1,0,1>, // 00010 
   labeling<0,0,0,0,1,0,1,0,1,1>, // 00001 

   labeling<0,0,0,1,1,1,1,1,1,0>  // 00011 
   */
   >;

using max_cut_5_clique_factor = labeling_factor< max_cut_5_clique_labelings, true >;

using max_cut_4_5_clique_message_0123 = labeling_message< max_cut_4_clique_labelings, max_cut_5_clique_labelings, 0,1,2,3,5,7 >; // 01->01, 02->02, 12->12, 03->03, 13->13, 23->23
using max_cut_4_5_clique_message_0124 = labeling_message< max_cut_4_clique_labelings, max_cut_5_clique_labelings, 0,1,2,4,6,8 >; // 01->01, 02->02, 12->12, 03->04, 13->14, 23->24
using max_cut_4_5_clique_message_0134 = labeling_message< max_cut_4_clique_labelings, max_cut_5_clique_labelings, 0,3,5,4,6,9 >; // 01->01, 02->03, 13->13, 03->04, 13->14, 23->34
using max_cut_4_5_clique_message_0234 = labeling_message< max_cut_4_clique_labelings, max_cut_5_clique_labelings, 1,3,7,4,8,9 >; // 01->02, 02->03, 12->23, 03->04, 13->24, 23->34
using max_cut_4_5_clique_message_1234 = labeling_message< max_cut_4_clique_labelings, max_cut_5_clique_labelings, 2,5,7,6,8,9 >; // 01->12, 02->13, 12->23, 03->14, 13->24, 23->34

/*
using max_cut_triplet_odd_bicycle_3_wheel_message_012 = labeling_message< max_cut_triplet_labelings, max_cut_odd_bicycle_3_wheel_labelings, 0,1,2 >;
using max_cut_triplet_odd_bicycle_3_wheel_message_013 = labeling_message< max_cut_triplet_labelings, max_cut_odd_bicycle_3_wheel_labelings, 0,3,5 >;
using max_cut_triplet_odd_bicycle_3_wheel_message_014 = labeling_message< max_cut_triplet_labelings, max_cut_odd_bicycle_3_wheel_labelings, 0,4,6 >;
using max_cut_triplet_odd_bicycle_3_wheel_message_023 = labeling_message< max_cut_triplet_labelings, max_cut_odd_bicycle_3_wheel_labelings, 1,3,7 >;
using max_cut_triplet_odd_bicycle_3_wheel_message_024 = labeling_message< max_cut_triplet_labelings, max_cut_odd_bicycle_3_wheel_labelings, 1,4,8 >;
using max_cut_triplet_odd_bicycle_3_wheel_message_034 = labeling_message< max_cut_triplet_labelings, max_cut_odd_bicycle_3_wheel_labelings, 3,4,9 >;
using max_cut_triplet_odd_bicycle_3_wheel_message_123 = labeling_message< max_cut_triplet_labelings, max_cut_odd_bicycle_3_wheel_labelings, 2,5,7 >;
using max_cut_triplet_odd_bicycle_3_wheel_message_124 = labeling_message< max_cut_triplet_labelings, max_cut_odd_bicycle_3_wheel_labelings, 2,6,8 >;
using max_cut_triplet_odd_bicycle_3_wheel_message_134 = labeling_message< max_cut_triplet_labelings, max_cut_odd_bicycle_3_wheel_labelings, 5,6,9 >;
using max_cut_triplet_odd_bicycle_3_wheel_message_234 = labeling_message< max_cut_triplet_labelings, max_cut_odd_bicycle_3_wheel_labelings, 7,8,9 >;
*/

} // end namespace LP_MP

#endif // LP_MP_MAX_CUT_FACTORS
