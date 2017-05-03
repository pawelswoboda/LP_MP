#ifndef LP_MP_MULTICUT_FACTORS_MESSAGES_HXX
#define LP_MP_MULTICUT_FACTORS_MESSAGES_HXX

#include "LP_MP.h"
#include "factors/labeling_list_factor.hxx"
#include "graph.hxx"
#include "union_find.hxx"

namespace LP_MP {

using multicut_edge_labelings = labelings< labeling<1> >;
using multicut_edge_factor = labeling_factor< multicut_edge_labelings, true, true >;

using multicut_triplet_labelings = labelings<
   labeling<0,1,1>,
   labeling<1,0,1>,
   labeling<1,1,0>,
   labeling<1,1,1>
      >;
using multicut_triplet_factor = labeling_factor< multicut_triplet_labelings, true >;

using multicut_edge_triplet_message_0 = labeling_message< multicut_edge_factor, multicut_triplet_factor, 0 >;
using multicut_edge_triplet_message_1 = labeling_message< multicut_edge_factor, multicut_triplet_factor, 1 >;
using multicut_edge_triplet_message_2 = labeling_message< multicut_edge_factor, multicut_triplet_factor, 2 >;

/*

   ______i1_____
  /      |      \
 /       |       \
/   _____i3____   \
|  /           \  |
| /             \ |
i0_______________i2

   */
using multicut_odd_3_wheel_labelings = labelings<
   // two components
   // one node separated from the rest
   labeling<1,1,0,1,0,0>, // 1000
   labeling<1,0,1,0,1,0>, // 0100
   labeling<0,1,1,0,0,1>, // 0010
   labeling<0,0,0,1,1,1>, // 0001

   // two components of two nodes each
   labeling<0,1,1,1,1,0>, // 1100 // 4
   labeling<1,0,1,1,0,1>, // 0101 // 5
   labeling<1,1,0,0,1,1>, // 1001 // 6

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

using multicut_triplet_odd_3_wheel_message_012 = labeling_message< multicut_triplet_factor, multicut_odd_3_wheel_factor, 0,1,2>;
using multicut_triplet_odd_3_wheel_message_013 = labeling_message< multicut_triplet_factor, multicut_odd_3_wheel_factor, 0,3,4>;
using multicut_triplet_odd_3_wheel_message_023 = labeling_message< multicut_triplet_factor, multicut_odd_3_wheel_factor, 1,3,5>;
using multicut_triplet_odd_3_wheel_message_123 = labeling_message< multicut_triplet_factor, multicut_odd_3_wheel_factor, 2,4,5>;

/*
                       i1



               
                    i3_
                       \__
                          i4




   i0                                       i2
   */
// is also the full graph on 5 nodes
using multicut_odd_bicycle_3_wheel_labelings = labelings<
// two components
// one node separated from the rest
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
labeling<0,1,1,0,1,1,0,0,1,1>, // 00101
labeling<0,0,0,1,1,1,1,1,1,0>, // 00011

// three conponents
// two singleton components
labeling<1,1,1,1,1,1,1,0,0,0>, // 01222
labeling<1,1,1,1,1,0,0,1,1,0>, // 02122
labeling<1,1,0,1,1,1,0,1,0,1>, // 02212
labeling<1,1,0,1,1,0,1,0,1,1>, // 02221
labeling<1,1,1,0,0,1,1,1,1,0>, // 20122
labeling<1,0,1,1,0,1,1,1,0,1>, // 20212
labeling<1,0,1,0,1,1,1,0,1,1>, // 20221
labeling<0,1,1,1,0,1,0,1,1,1>, // 22012
labeling<0,1,1,0,1,0,1,1,1,1>, // 22021
labeling<0,0,0,1,1,1,1,1,1,1>, // 22201
// one singleton component, two two-element components
labeling<1,1,0,1,1,1,1,1,1,0>, // 01122
labeling<1,1,1,1,1,0,1,1,0,1>, // 01212
labeling<1,1,1,1,1,1,0,1,0,1>, // 01221
labeling<1,0,1,1,1,1,1,1,1,0>, // 10122
labeling<1,1,1,0,1,1,1,1,0,1>, // 10212
labeling<1,1,1,1,0,1,1,1,0,1>, // 10221
labeling<0,1,1,1,1,1,1,1,1,0>, // 11022
labeling<1,1,1,0,1,1,0,1,1,1>, // 12012
labeling<1,1,1,1,0,0,1,1,1,0>, // 12021
labeling<0,1,1,1,1,1,1,1,0,1>, // 11202
labeling<1,0,1,1,1,1,0,1,1,1>, // 12102
labeling<1,1,0,1,0,1,1,1,1,1>, // 12201
labeling<0,1,1,1,1,1,1,0,1,1>, // 11220
labeling<1,0,1,1,1,0,1,1,1,1>, // 12120
labeling<1,1,0,0,1,1,1,1,1,1>, // 12210

// four components
labeling<0,1,1,1,1,1,1,1,1,1>,
labeling<1,0,1,1,1,1,1,1,1,1>,
labeling<1,1,0,1,1,1,1,1,1,1>,
labeling<1,1,1,0,1,1,1,1,1,1>,
labeling<1,1,1,1,0,1,1,1,1,1>,
labeling<1,1,1,1,1,0,1,1,1,1>,
labeling<1,1,1,1,1,1,0,1,1,1>,
labeling<1,1,1,1,1,1,1,0,1,1>,
labeling<1,1,1,1,1,1,1,1,0,1>,
labeling<1,1,1,1,1,1,1,1,1,0>,
// five components
labeling<1,1,1,1,1,1,1,1,1,1>
   >;
using multicut_odd_bicycle_3_wheel_factor = labeling_factor< multicut_odd_bicycle_3_wheel_labelings, true >;

using multicut_odd_3_wheel_odd_bicycle_message0123 = labeling_message< multicut_odd_3_wheel_factor, multicut_odd_bicycle_3_wheel_factor, 0,1,3,2,5,7 >;
using multicut_odd_3_wheel_odd_bicycle_message0124 = labeling_message< multicut_odd_3_wheel_factor, multicut_odd_bicycle_3_wheel_factor, 0,1,4,2,6,8 >;
using multicut_odd_3_wheel_odd_bicycle_message0134 = labeling_message< multicut_odd_3_wheel_factor, multicut_odd_bicycle_3_wheel_factor, 0,3,4,5,6,9 >;
using multicut_odd_3_wheel_odd_bicycle_message0234 = labeling_message< multicut_odd_3_wheel_factor, multicut_odd_bicycle_3_wheel_factor, 1,3,4,7,8,9 >;
using multicut_odd_3_wheel_odd_bicycle_message1234 = labeling_message< multicut_odd_3_wheel_factor, multicut_odd_bicycle_3_wheel_factor, 2,5,6,7,8,9 >;

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

// class for holding necessary information such that a given edge knows whether it can be locally set to 0 or 1 given partial labeling
enum class multicut_possible_labels {arbitrary, cut, non_cut};
class multicut_rounding {
public:
   multicut_rounding(const INDEX no_nodes, const INDEX no_arcs, const std::vector<INDEX>& no_outgoing_arcs) :
      uf(no_nodes), g(no_nodes, no_arcs, no_outgoing_arcs), bfs(no_nodes)
   {}

   multicut_possible_labels possible_labels(const INDEX i, const INDEX j) 
   {
      if(uf.connected(i,j)) {
         return multicut_possible_labels::non_cut;
      } else if(bfs.path_exists(i,j,g)) { // check whether distance w.r.t. cut edges is 1 
         return multicut_possible_labels::cut;
      } else {
         return multicut_possible_labels::arbitrary;
      }
   }

   void set(const INDEX i, const INDEX j, const bool cut) 
   {
      if(cut) {
         g.add_edge(i,j,1.0);
      } else {
         uf.merge(i,j);
         g.add_edge(i,j,0.0);
      } 
   }
   private:
   UnionFind uf; // records whether two nodes are in the same component
   Graph g;
   multicut_path_search bfs; 
};

// multicut edge factor additionally holding pointer to edge value which it records when doing rounding. This will be the edge cost for the subsequent primal rounding with KL+j
class rounding_multicut_edge_factor : public multicut_edge_factor 
{
public:
   rounding_multicut_edge_factor(const INDEX _i, const INDEX _j) :
      multicut_edge_factor(),
      i(_i), j(_j), mc_rounding(nullptr)
   {
      assert(i < j);
   }

   void set_rounding(multicut_rounding* m)
   {
      mc_rounding = m;
   }

   void MaximizePotentialAndComputePrimal()
   {
      auto r = mc_rounding->possible_labels(i,j);
      if(r == multicut_possible_labels::cut) {
         this->primal()[0] = true;
      } else if(r == multicut_possible_labels::non_cut) {
         this->primal()[0] = false;
      } else {
         assert(r == multicut_possible_labels::arbitrary);
         this->primal()[0] =  (*this)[0] < 0 ? true : false;
      }
      mc_rounding->set(i,j,this->primal()[0]);
   }

   //REAL get_rounding_cost() const { return edge_val_; }
   REAL get_rounding_cost() const { return (*this)[0]; }
private:
   multicut_rounding* mc_rounding;
   const INDEX i;
   const INDEX j;
};

} // end namespace LP_MP 

#endif // LP_MP_MULTICUT_FACTORS_MESSAGES_HXX
