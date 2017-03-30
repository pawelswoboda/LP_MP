#ifndef LP_MP_MLT_FACTORS_MESSAGES_HXX
#define LP_MP_MLT_FACTORS_MESSAGES_HXX

#include "config.hxx"
#include "factors/labeling_list_factor.hxx"

namespace LP_MP {

// states: last edge is edge from t+1 frame
//	011
//	101
//	110
//	111
//	001
using mlt_triplet_labelings = labelings<
	labeling<0,1,1>,
	labeling<1,0,1>,
	labeling<1,1,0>,
	labeling<1,1,1>,
	labeling<0,0,1>
>;
using mlt_triplet_factor = labeling_factor<mlt_triplet_labelings, true>;

using edge_mlt_triplet_message_0 = labeling_message<multicut_edge_labelings, mlt_triplet_labelings,0>;
using edge_mlt_triplet_message_1 = labeling_message<multicut_edge_labelings, mlt_triplet_labelings,1>;
using edge_mlt_triplet_message_2 = labeling_message<multicut_edge_labelings, mlt_triplet_labelings,2>;


// states:
// first inter-frame edges, then intra-frame edges for t+1

// 000 011 | 0
// 000 101 | 1
// 000 110 | 2

// 100 101 | 3
// 100 111 | 4
// 010 110 | 5
// 010 111 | 6
// 001 011 | 7
// 001 111 | 8

// 011 101 | 9
// 011 111 | 10
// 101 110 | 11
// 101 111 | 12
// 011 011 | 13
// 011 111 | 14

// 111 000 | 15
// 111 011 | 16
// 111 101 | 17
// 111 110 | 18
// 111 111 | 19
using mlt_bifurcation_labelings = labelings<
	labeling<0,0,0, 0,1,1>,
	labeling<0,0,0, 0,1,1>,
	labeling<0,0,0, 0,1,1>,

	labeling<1,0,0, 1,0,1>,
	labeling<1,0,0, 1,1,1>,
	labeling<0,1,0, 1,1,0>,
	labeling<0,1,0, 1,1,1>,
	labeling<0,0,1, 0,1,1>,
	labeling<0,0,1, 1,1,1>,

	labeling<0,1,1, 1,0,1>,
	labeling<0,1,1, 1,1,1>,
	labeling<1,0,1, 1,1,0>,
	labeling<1,0,1, 1,1,1>,
	labeling<0,1,1, 0,1,1>,
	labeling<0,1,1, 1,1,1>,

	labeling<1,1,1, 0,0,0>,
	labeling<1,1,1, 0,1,1>,
	labeling<1,1,1, 1,0,1>,
	labeling<1,1,1, 1,1,0>,
	labeling<1,1,1, 1,1,1>
>;
using mlt_bifurcation_factor = labeling_factor<mlt_bifurcation_labelings,true>;

using mlt_triplet_bifurcation_message_0 = labeling_message<mlt_triplet_labelings, mlt_bifurcation_labelings, 0,1,3>;
using mlt_triplet_bifurcation_message_1 = labeling_message<mlt_triplet_labelings, mlt_bifurcation_labelings, 1,2,4>;
using mlt_triplet_bifurcation_message_2 = labeling_message<mlt_triplet_labelings, mlt_bifurcation_labelings, 0,2,5>;
using mc_triplet_bifurcation_message = labeling_message<multicut_triplet_labelings, mlt_bifurcation_labelings, 3,5,4>;

}

#endif