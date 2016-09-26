/*
 *
 *  Created by Amir Globerson and David Sontag on 8/10/08.
 *  Updated by David Sontag, Do Kook Choe, and Yitao Li in 2012.
 *  Copyright 2008 MIT, 2012 NYU. All rights reserved.
 *
 *  Adapted to LP_MP by Paul Swoboda (only 2012 code)
 *
 */
#ifndef LP_MP_CYCLE_INEQUALITIES_HXX
#define LP_MP_CYCLE_INEQUALITIES_HXX

#include <stdlib.h>
#include <vector>
#include <algorithm>
#include <set>
#include <list>
#include <map>
#include <queue>

#include "marray.hxx" // check how to replace this
#include "config.hxx"

#define DEBUG_MODE 1

#define Inf 9999999999.9

namespace LP_MP {

   // do zrobinia: 
   // - we minimize instead of maximizing. Replace max by min everywhere? Currently, negative potentials are taken.
   // - LP_MP tightening takes fixed epsilon, and searches for violated inequalities regarding this value. Sontag et al code does not do this, but searches for violated inequalities over all epsilon. Change this as well.
   // - use hash instead of std::map
   // - graph data for cycle search are held in manually allocated memory. Replace by std containers whenever possible. Possibility of memory leaks in code

template<typename MRF_CONSTRUCTOR>
class Cycle 
{
public:
   //typedef std::map<std::pair<int, int>, int> mapType;
   typedef std::vector<std::vector<std::pair<int, REAL> > > adj_type;

   // Structure for storing a candidate triplet cluster for tightening
   struct TripletCluster {
      REAL bound;
      int i,j,k;
      int ij_intersect_loc, jk_intersect_loc, ki_intersect_loc;
      std::vector<std::set<SIGNED_INDEX> > labels;

      bool operator <(const TripletCluster & rhs) const {
         return bound < rhs.bound;
      }
   };

   Cycle (const MRF_CONSTRUCTOR& mrf);
   ~Cycle() {};

   // UAI 2008
   // ADD_TRIPLET_FUNCTION must be function of type
   // <void(const SIGNED_INDEX,const SIGNED_INDEX, const SIGNED_INDEX, const std::vector<SIGNED_INDEX>, const std::vector<SIGNED_INDEX>, const std::vector<SIGNED_INDEX>)
   template<typename ADD_TRIPLET_FUNCTION>
   int TightenTriplet(
         ADD_TRIPLET_FUNCTION addTripletFun,
         int nclus_to_add_max, 
         const REAL epsilon,
         std::map<std::vector<int>, bool >& triplet_set);

   // UAI 2012
   template<typename ADD_TRIPLET_FUNCTION>
   int TightenCycle(
         ADD_TRIPLET_FUNCTION addTripletFun,
         const int nclus_to_add,  
         std::vector<int>& projection_imap_var,
         std::vector<std::vector<int> >& partition_imap,
         std::vector<std::list<int> >& cycle_set,
         const int method);


private:
   // UAI 2008
   template<typename BELIEF_ARRAY>
   REAL maximizeIndependently(const BELIEF_ARRAY & beliefs) const;
   template<typename BELIEF_ARRAY>
   REAL maximizeCycle(const BELIEF_ARRAY& beliefs) const;//, std::vector<std::set<SIGNED_INDEX> >& cyclePi) const;
   template<typename BELIEF_ARRAY>
   REAL maximizeTriangle(const BELIEF_ARRAY& beliefs) const;

   // UAI 2012
   static bool edge_sorting(std::list<int> i, std::list<int> j);
   static bool edge_sorting3(std::pair<std::pair<int, int>, REAL> i, std::pair<std::pair<int, int>, REAL> j);
   void delete_projection_graph(int num_vars, std::vector<std::vector<int> > &projection_map, std::vector<int> &projection_imap_var, adj_type &projection_adjacency_list, REAL* &array_of_sij);
   REAL find_smn( const std::vector<int>& partition_i, const std::vector<int>& partition_j, const int factorId);
   //REAL find_smn_state_i(const SIGNED_INDEX single_i, const std::vector<int>& partition_j, const int factorId, std::vector<std::vector<REAL> >& max_i_bij_not_xi);
   //REAL find_smn_state_j(const SIGNED_INDEX single_j, const std::vector<int>& partition_j, const int factorId, std::vector<std::vector<REAL> >& max_j_bij_not_xj);
   REAL find_smn_state_i(const SIGNED_INDEX single_i, const std::vector<int>& partition_j, const int factorId, const marray::Marray<REAL>& max_i_bij_not_xi);
   REAL find_smn_state_j(const SIGNED_INDEX single_j, const std::vector<int>& partition_j, const int factorId, const marray::Marray<REAL>& max_j_bij_not_xj);
   void find_partition(std::vector<std::map<std::vector<int>, int> >& partition_set,  const int factorId, int index_i, int index_j);
   std::tuple<marray::Marray<REAL>, marray::Marray<REAL>, marray::Marray<REAL>>
      calculate_cond_min_beliefs(const INDEX factor_id);
   int create_expanded_projection_graph(
      std::vector<int>& projection_imap_var,
      std::vector<std::vector<std::pair<int, REAL> > >& projection_adjacency_list, 
      std::map<std::pair<int, int>, REAL>&projection_edge_weights,
      REAL* &array_of_sij,
      int& array_of_sij_size,
      std::vector<std::vector<int> >& partition_imap);
   void create_k_projection_graph(
      std::vector<std::vector<int> > &projection_map, 
      int& num_projection_nodes, 
      std::vector<int> &projection_imap_var, 
      std::vector<std::vector<int> > &partition_imap,
      std::map<std::pair<int, int>, REAL>& projection_edge_weights,
      adj_type &projection_adjacency_list,
      REAL* &array_of_sij,
      int& array_of_sij_size);
   REAL find_optimal_R(adj_type &projection_adjacency_list, REAL* &array_of_sij, int array_of_sij_size);
   // do zrobienia: not nice: use C++ random shuffle
   int* random_permutation(int n);
   void FindCycles(std::vector<std::list<int> > &cycle_set, REAL optimal_R, int ncycles_to_add, adj_type &projection_adjacency_list);
   void shortcut(std::list<int> &cycle, std::vector<int> &projection_imap_var, std::map<std::pair<int, int>, REAL>& projection_edge_weights, int& num_projection_nodes);

   template<typename ADD_TRIPLET_FUNCTION>
   int add_cycle( 
         ADD_TRIPLET_FUNCTION addTripletFun,
         //std::function<void(const int,const int, const int, const std::vector<SIGNED_INDEX>, const std::vector<SIGNED_INDEX>, const std::vector<SIGNED_INDEX>)> addTripletFun,
         std::list<int> &cycle, 
         std::vector<int> &projection_imap_var, 
         std::vector<std::vector<int> >& partition_imap,
         std::map<std::vector<int>, bool >& triplet_set, 
         int& num_projection_nodes);


   const MRF_CONSTRUCTOR& gm_;
};

template<class MRF_CONSTRUCTOR>
Cycle<MRF_CONSTRUCTOR>::Cycle
( const MRF_CONSTRUCTOR& gm)
   : gm_(gm)
{}


/////////////////////////////////////////////////////////////////////////////////
// Code for union-find data structure (used by FindPartition)
// do zrobienia: use existing union find class

struct Node { //node for union-find
  int bit; //10, 01, or 11
  int rank; 
  Node* parent;
  int i_size; //number of elements in i 
  int j_size; //number of elements in j 
  
  Node(int b) : bit(b), rank(0), parent(this) {
    if (b == 2) { //represent elements of i
      i_size = 1;
      j_size = 0;
    }
    else { //reprsent elements of j
      i_size = 0;
      j_size = 1;
    }
  }
};
  
inline Node* find(Node* n) {
	if (n != n->parent) {
		n->parent = find(n->parent);
	}
	return n->parent;
} 


inline Node* merge(Node* x, Node* y) {
  Node* root_x = find(x);
  Node* root_y = find(y);
  if (root_x == root_y) { //x and y have the same head
    return root_x;
  }
  
  if (root_x->rank > root_y->rank) {
    root_y->parent = root_x; //head x
    root_x->i_size += root_y->i_size; //add number of elements of i to the head node
    root_x->j_size += root_y->j_size; //add number of elements of j to the head node    
    return root_x;
  }
  else {
    root_x->parent = root_y; //head y
    if (root_x->rank == root_y->rank) {
      ++(root_y->rank);
    }
    root_y->i_size += root_x->i_size; //add number of elements of i to the head node 
    root_y->j_size += root_x->j_size; //add number of elements of j to the head node
    return root_y;
  }  
}

/////////////////////////////////////////////////////////////////////////////////
// Code to evaluate how good a cycle cluster is
 
template<class MRF_CONSTRUCTOR>
template<typename BELIEF_ARRAY>
REAL
Cycle<MRF_CONSTRUCTOR>::maximizeIndependently(const BELIEF_ARRAY& beliefs) const
{
	REAL sum=0.0;
	for(int i=0; i < beliefs.size(); i++) {
		sum += *std::max_element(beliefs[i].begin(), beliefs[i].end());
	}
	return sum;
}

/*
REAL getValCycle(std::vector<MulDimArr*> & beliefs, std::vector<bool> & b_transpose, std::vector<int> & assignments)
{
	REAL sum=0.0;
	std::vector<int> inds; inds.push_back(-1); inds.push_back(-1); // temp

	// All except the last edge
	for(int i=0; i < beliefs.size()-1; i++)
	{
		inds[b_transpose[i]?1:0] = assignments[i];
		inds[b_transpose[i]?0:1] = assignments[i+1];
		sum += beliefs[i]->GetVal(inds);
	}

	// Now do last edge
	inds[b_transpose[beliefs.size()-1]?1:0] = assignments[beliefs.size()-1];
	inds[b_transpose[beliefs.size()-1]?0:1] = assignments[0];
	sum += beliefs[beliefs.size()-1]->GetVal(inds);

	return sum;
}
*/

// faster alternative to maximizeCycle for three elements
template<typename MRF_CONSTRUCTOR>
template<typename BELIEF_ARRAY>
REAL 
Cycle<MRF_CONSTRUCTOR>::maximizeTriangle(const BELIEF_ARRAY& beliefs) const
{
   assert(beliefs.size() == 3);
   REAL max_val = -std::numeric_limits<REAL>::infinity();

   // Fix value of the first variable
   const INDEX dim1 = beliefs[0].shape(0);
   const INDEX dim2 = beliefs[0].shape(1);
   const INDEX dim3 = beliefs[1].shape(1);
   assert(beliefs[0].shape(1) == beliefs[1].shape(0) && beliefs[1].shape(1) == beliefs[2].shape(0) && beliefs[2].shape(1) == beliefs[0].shape(0));
   for(INDEX i1=0; i1<dim1; ++i1) {
      for(INDEX i2=0; i2<dim2; ++i2) {
         for(INDEX i3=0; i3<dim3; ++i3) {
            max_val = std::max(max_val, beliefs[0](i1,i2) + beliefs[1](i2,i3) + beliefs[2](i3,i1));
         }
      }
   }
   return max_val;
}


// Choose partition such that all labelings which have value > independent bound are included in the partiton cyclePi
template<typename MRF_CONSTRUCTOR>
template<typename BELIEF_ARRAY>
REAL
Cycle<MRF_CONSTRUCTOR>::maximizeCycle(const BELIEF_ARRAY& beliefs) const
{
   double max_val = -std::numeric_limits<REAL>::max();//-Inf;

   // Fix value of the first variable
   int first_var_size = beliefs[0].shape(0);
   int second_var_size = beliefs[0].shape(1);
   for(int vo=0; vo < first_var_size; vo++)
   {
      std::array<SIGNED_INDEX,2> inds;
      inds[0] = vo;

      // Do first edge (construct initial field)
      std::vector<REAL> field;
      for(int v2=0; v2 < second_var_size; v2++)
      {
         inds[1] = v2;
         field.push_back(beliefs[0](inds[0],inds[1]));
      }

      // Go over rest of edges, except last (which has to be treated specially)
      for(int i=1; i < beliefs.size()-1; i++)
      {
         std::vector<REAL> new_field;
         for(int v2=0; v2 < beliefs[i].shape(1); v2++) //->m_base_sizes[b_transpose[i]?0:1]; v2++)
         {
            inds[0] = -1; inds[1] = -1;
            inds[1] = v2;

            // Take max
            double tmp_max_val = -Inf;
            for(int v1=0; v1 < field.size(); v1++)
            {
               //inds[b_transpose[i]?1:0] = v1;
               inds[0] = v1;
               tmp_max_val = std::max(tmp_max_val, field[v1] + beliefs[i](inds[0],inds[1]));//->GetVal(inds));
            }
            new_field.push_back(tmp_max_val);
         }
         // do zrobienia: std::swap
         field.clear(); // necessary?
         field = new_field;
      }

      // Do last edge (fix endpoint value to vo)
      inds[0] = -1;
      inds[1] = vo;

      // Take max
      double tmp_max_val = -Inf;
      for(int v1=0; v1 < field.size(); v1++)
      {
         inds[0] = v1;
         tmp_max_val = std::max(tmp_max_val, field[v1] + beliefs[beliefs.size()-1](inds[0],inds[1]));//->GetVal(inds));
      }

      max_val = std::max(max_val, tmp_max_val);
   }
   return max_val;
}

/////////////////////////////////////////////////////////////////////////////////





// do zrobienia: add eps
// do zrobienia> templatize for beliefs array, such that std::array is possible
// do zrobienia: use special functions for cycle of length 3 or 4
/*
template<typename MRF_CONSTRUCTOR>
REAL
Cycle<MRF_CONSTRUCTOR>::maximizeCycle(const std::vector<marray::View<REAL,true> > & beliefs, std::vector<std::set<SIGNED_INDEX> >& cyclePi) const
{
   const REAL bound_indep = maximizeIndependently(beliefs);

   REAL max_val = -Inf;
   cyclePi.clear();
   cyclePi.resize(beliefs.size()); // return labels on cycle which go beyond the independent bound

   // Fix value of the first variable
   const int first_var_size = beliefs[0].shape(0);
   const int second_var_size = beliefs[0].shape(1);
   for(int vo=0; vo < first_var_size; vo++) {
      std::vector<SIGNED_INDEX> cycleLabelingTmp(beliefs.size()); // do zrobienia: better array type

      // Do first edge (construct initial field)
      std::vector<std::vector<REAL> > field(1); // do zrobieni: template array type
      for(int v2=0; v2 < second_var_size; v2++) {
         field[0].push_back(beliefs[0](vo,v2));
      }

      // Go over rest of edges, except last (which has to be treated specially)
      for(int i=1; i < beliefs.size()-1; i++) {
         std::vector<REAL> new_field(beliefs[i].shape(1));
         for(int v2=0; v2 < beliefs[i].shape(1); v2++) {
            // Take max
            REAL tmp_max_val = -Inf;
            for(int v1=0; v1 < field.back().size(); v1++) {
               tmp_max_val = std::max(tmp_max_val, field[i-1][v1]+beliefs[i](v1,v2));
            }
            new_field[v2] = tmp_max_val;
         }
         field.push_back(new_field);
      }

      // Do last edge (fix endpoint value to vo)
      // Take max
      REAL tmp_max_val = -Inf;
      for(int v1=0; v1 < field.back().size(); v1++) {
         REAL cur_val = field.back()[v1]+beliefs.back()(v1,vo);
         if(tmp_max_val < cur_val) {
            tmp_max_val = std::max(tmp_max_val, cur_val);
         }
         if(cur_val + eps < bound_indep) {
            cyclePi.back().insert(v1);
         }
         field.back()[v1] = cur_val;
      }

      if(max_val < tmp_max_val) {
         max_val = std::max(max_val, tmp_max_val);
      }
      }
      return max_val;
      {

      if(bound_indep > eps + tmp_max_val) {
         cyclePi[0].insert(vo);

         // go backwards and include all labels for which there is a labeling with value > bound_indep
         for(int i=beliefs.size()-2; i>0; i--) {

            std::vector<REAL> new_field(beliefs[i].shape(0),-Inf);
            for(int v1=0; v1 < beliefs[i].shape(0); v1++) {
               // Take max
               REAL tmp_max_val = -Inf;
               for(int v2=0; v2 < beliefs[i].shape(1); v2++) {
                  tmp_max_val = std::max(tmp_max_val, field[i][v2]+beliefs[i+1](v1,v2));
                  new_field[v1] = std::max(new_field[v1], field[i][v2] - beliefs[i+1](v1,v2));
               }
            }
            field[i-1] = new_field;

            for(int v1=0; v1<field[i].size(); v1++) {
               for(int v2=0; v2 < beliefs[i].shape(1); v2++) {
                  const REAL cur_val = field[i-1][v1] + beliefs[i](v1,v2) + field[i][v2];
                  if(cur_val + eps < bound_indep) {
                     cyclePi[i].insert(v1);
                  }
               } 
            }
         }
         for(size_t i=0; i<cyclePi[i].size(); i++) {
            assert( cyclePi[i].size() > 0 );
         }
      }

   }


   if(beliefs.size() == 3) {
      std::vector<std::set<SIGNED_INDEX> > cyclePiTest(3);
      for(size_t x_i=0; x_i<beliefs[0].shape(0); x_i++) {
         for(size_t x_j=0; x_j<beliefs[1].shape(0); x_j++) {
            for(size_t x_k=0; x_k<beliefs[2].shape(0); x_k++) {
               const REAL cost = beliefs[0](x_i,x_j) + beliefs[1](x_j,x_k) + beliefs[2](x_k,x_i);
               if(cost + eps < bound_indep) {
                  cyclePiTest[0].insert(x_i);
                  cyclePiTest[1].insert(x_j);
                  cyclePiTest[2].insert(x_k);
               }
            }
         }
      }
      assert(cyclePiTest == cyclePi); // possibly this is due to large values in pairwise potentials coming from infty on diagonals in matching problems.
   }
   return bound_indep - max_val;
}
*/


/////////////////////////////////////////////////////////////////////////////////
// Implementation of UAI 2008 algorithm (just for triplets; square functionality removed)
// do zrobienia: add square functionality again?
//               only store triplets with bound > CLUSTER_THR, not all! Only those are added, the rest is discarded either way.
template<typename MRF_CONSTRUCTOR>
template<typename ADD_TRIPLET_FUNCTION>
int
Cycle<MRF_CONSTRUCTOR>::TightenTriplet(
         //std::function<void(const SIGNED_INDEX,const SIGNED_INDEX, const SIGNED_INDEX, const std::vector<SIGNED_INDEX>, const std::vector<SIGNED_INDEX>, const std::vector<SIGNED_INDEX>)> addTripletFun,
         ADD_TRIPLET_FUNCTION addTripletFun,
         int nclus_to_add_max, 
         const REAL epsilon,
         std::map<std::vector<int>, bool >& triplet_set) 
{
   int nClustersAdded = 0;
   int nNewClusters = 0;

   if(DEBUG_MODE)
      std::cout << "Doing pre-processing for adding triplet clusters." << std::endl;

   // Initialize adjacency list (filled in later) TODO: only do this when needed
   //std::vector<int> adjacency_list[mplp.m_var_sizes.size()];
   std::vector<std::vector<int> > adjacency_list(gm_.GetNumberOfVariables());

   // Construct adjacency list for the graph
   // Iterate over all of the edges (we do this by looking at the edge intersection sets)
   for(size_t factorId=0; factorId<gm_.GetNumberOfPairwiseFactors(); factorId++) {
      auto vars = gm_.GetPairwiseVariables(factorId);
      // Get the two nodes i & j
      const size_t i=std::get<0>(vars);
      const size_t j=std::get<1>(vars);
      assert(i<j);
      adjacency_list[i].push_back(j);
      adjacency_list[j].push_back(i);
   }

   // Sort the adjacency list, for fast intersections later
   for(int i=0; i < adjacency_list.size(); i++) {
      std::sort(adjacency_list[i].begin(), adjacency_list[i].end());
   }

  // Count the number of triangles
  std::vector<int>::iterator intersects_iter_end;
  std::vector<int> commonNodes(gm_.GetNumberOfVariables());
  for(size_t factorId=0; factorId<gm_.GetNumberOfPairwiseFactors(); factorId++) {
     auto vars = gm_.GetPairwiseVariables(factorId);
     // Get the two nodes i & j
     const size_t i=std::get<0>(vars);
     const size_t j=std::get<1>(vars);

     // Now find all neighbors of both i and j to see where the triangles are
     intersects_iter_end = set_intersection(adjacency_list[i].begin(), adjacency_list[i].end(), adjacency_list[j].begin(), adjacency_list[j].end(), commonNodes.begin());

     for(std::vector<int>::const_iterator n=commonNodes.begin(); n != intersects_iter_end; ++n) {
        // Since a triplet shows up three times as an edge plus
        // a node, we only consider it for the case when i<j<k
        if(j < *n)
           nNewClusters++;
     }
  }

  if(nNewClusters == 0) {
    if(DEBUG_MODE)
      std::cout << "nNewClusters = 0. Returning." << std::endl;
    return 0;
  }

  if(DEBUG_MODE)
    std::cout << "Looking for triangle clusters to add (" << nNewClusters << " triplets) " << std::endl;

  // TODO: put this elsewhere so that the space isn't re-allocated continuously?
  // Enumerate over all of the edges
  std::vector<TripletCluster> newCluster(nNewClusters);
			
  int index=0;
			
  // Iterate over all of the edge intersection sets
  // do zrobienia: parallelize
  for(size_t factorId=0; factorId<gm_.GetNumberOfPairwiseFactors(); factorId++) {
     auto vars = gm_.GetPairwiseVariables(factorId);
     // Get the two nodes i & j
     const INDEX i=std::get<0>(vars);
     const INDEX j=std::get<1>(vars);

     // Now find all neighbors of both i and j to see where the triangles are
     // TEMP TEMP -- fails at i=0, j=1, on i==3.
     intersects_iter_end = set_intersection(adjacency_list[i].begin(), adjacency_list[i].end(), adjacency_list[j].begin(), adjacency_list[j].end(), commonNodes.begin());

     for(std::vector<int>::const_iterator n=commonNodes.begin(); n != intersects_iter_end; ++n) {
        INDEX k = *n;

        // Since a triplet shows up three times as an edge plus
        // a node, we only consider it for the case when i<j<k 
        if(!(j<k))
           continue;

        std::array<INDEX,3> inds = {{i,j,k}};
        newCluster[index].i = inds[0];
        newCluster[index].j = inds[1];
        newCluster[index].k = inds[2];

        // Find the intersection sets for this triangle
        newCluster[index].ij_intersect_loc = gm_.GetPairwiseFactorId(inds[0],inds[1]);
        newCluster[index].jk_intersect_loc = gm_.GetPairwiseFactorId(inds[1],inds[2]);
        newCluster[index].ki_intersect_loc = gm_.GetPairwiseFactorId(inds[0],inds[2]);

        // Construct the beliefs for each edge, which will be maximized below
        std::array<marray::Marray<REAL>,3> beliefs;

        SIGNED_INDEX shape_ij[] = {SIGNED_INDEX(gm_.GetNumberOfLabels(inds[0])), SIGNED_INDEX(gm_.GetNumberOfLabels(inds[1]))};
        beliefs[0] = marray::Marray<REAL>(shape_ij, shape_ij+2);
        //gm_[ newCluster[index].ij_intersect_loc ].copyValues(beliefs[0].begin());
        //gm_.CopyPairwiseValues( newCluster[index].ij_intersect_loc, beliefs[0].begin() );
        for(INDEX label_i=0; label_i<gm_.GetNumberOfLabels(i); ++label_i) {
           for(INDEX label_j=0; label_j<gm_.GetNumberOfLabels(j); ++label_j) {
              beliefs[0](label_i,label_j) = -gm_.GetPairwiseValue(newCluster[index].ij_intersect_loc, label_i, label_j);
           }
        }

        SIGNED_INDEX shape_jk[] = {SIGNED_INDEX(gm_.GetNumberOfLabels(inds[1])), SIGNED_INDEX(gm_.GetNumberOfLabels(inds[2]))};
        beliefs[1] = marray::Marray<REAL>(shape_jk, shape_jk+2);
        //gm_[ newCluster[index].jk_intersect_loc ].copyValues(beliefs[1].begin());
        //gm_.CopyPairwiseValues( newCluster[index].jk_intersect_loc, beliefs[1].begin() );
        for(INDEX label_j=0; label_j<gm_.GetNumberOfLabels(j); ++label_j) {
           for(INDEX label_k=0; label_k<gm_.GetNumberOfLabels(k); ++label_k) {
              beliefs[1](label_j,label_k) = -gm_.GetPairwiseValue(newCluster[index].jk_intersect_loc, label_j, label_k);
           }
        }

        SIGNED_INDEX shape_ik[] = {SIGNED_INDEX(gm_.GetNumberOfLabels(inds[0])), SIGNED_INDEX(gm_.GetNumberOfLabels(inds[2]))};
        beliefs[2] = marray::Marray<REAL>(shape_ik, shape_ik+2);
        //gm_[ newCluster[index].ki_intersect_loc ].copyValues(beliefs[2].begin());
        //gm_.CopyPairwiseValues( newCluster[index].ki_intersect_loc, beliefs[2].begin() );
        for(INDEX label_i=0; label_i<gm_.GetNumberOfLabels(i); ++label_i) {
           for(INDEX label_k=0; label_k<gm_.GetNumberOfLabels(k); ++label_k) {
              beliefs[2](label_i,label_k) = -gm_.GetPairwiseValue(newCluster[index].ki_intersect_loc, label_i, label_k);
           }
        }

        // this is due to transpose being a private member of beliefs type
        std::array<marray::View<REAL,true>,3> b; 
        b[0] = beliefs[0];
        b[1] = beliefs[1];
        b[2] = beliefs[2];
        b[2].transpose();
        //b[0].transpose();
        //b[1].transpose();
        //beliefs[3].transpose();

        //std::vector<std::set<SIGNED_INDEX> > cyclePi;
        const REAL boundIndep = maximizeIndependently(b);
        //const REAL boundCycle = maximizeCycle(b);
        const REAL boundCycle = maximizeTriangle(b);
        newCluster[index].bound = boundIndep - boundCycle;
        assert(newCluster[index].bound >=  - eps);
        index++;
     }
  }

  // TODO opt: have a class for a cluster, so we can have different types and sort by bound,
  //       choosing best one by bound.
  //       Make the sorting and adding independent of the type of graph...
			
  // Sort the clusters by the bound
  if(DEBUG_MODE)
    std::cout << "Optimized over all triangles, now add those with greatest bound\n";

  std::sort(newCluster.begin(), newCluster.end());

  if(DEBUG_MODE)
     std::cout << " -- Considered " << nNewClusters << " clusters, smallest bound " << newCluster[std::max(nNewClusters-nclus_to_add_max, 0)].bound << ", largest bound " << newCluster[nNewClusters-1].bound << std::endl;

  // Add the top nclus_to_add clusters to the relaxation
  for(int clusterId = nNewClusters-1; clusterId >= 0 && nClustersAdded < nclus_to_add_max; clusterId--) {
     std::cout << newCluster[clusterId].bound << " = bound\n";

     // Now add cluster ijk
     std::array<int,3> ijk_inds;
     //ijk_inds.push_back(newCluster[clusterId].i); ijk_inds.push_back(newCluster[clusterId].j); ijk_inds.push_back(newCluster[clusterId].k);
     ijk_inds[0] = newCluster[clusterId].i;
     ijk_inds[1] = newCluster[clusterId].j;
     ijk_inds[2] = newCluster[clusterId].k;

     bool tripletAdded = addTripletFun(ijk_inds[0], ijk_inds[1], ijk_inds[2]);
     if(tripletAdded) {
        nClustersAdded++;
        if(DEBUG_MODE) {
           std::cout << "Cluster added on nodes " << newCluster[clusterId].i << ", " << newCluster[clusterId].j << ", " << newCluster[clusterId].k << std::endl;
        }
     }

     if(newCluster[clusterId].bound < epsilon) {
        std::cout << "cluster with too small improvement found, abort\n";
        break;
     }
  }

  return nClustersAdded;
}



/////////////////////////////////////////////////////////////////////////////////
// Everything that follows is for the UAI 2012 cycle finding algorithm which
// can find arbitrary-length cycles to use in tightening the relaxation.

template<typename MRF_CONSTRUCTOR>
bool 
Cycle<MRF_CONSTRUCTOR>::edge_sorting(std::list<int> i, std::list<int> j) {
	return i.back() > j.back();
}

template<typename MRF_CONSTRUCTOR>
bool 
Cycle<MRF_CONSTRUCTOR>::edge_sorting3(std::pair<std::pair<int, int>, REAL> i, std::pair<std::pair<int, int>, REAL> j) {
  return i.second > j.second;
}

// De-allocate memory relating to the projection graph (TODO: finish writing this)
template<typename MRF_CONSTRUCTOR>
void 
Cycle<MRF_CONSTRUCTOR>::delete_projection_graph(int num_vars, std::vector<std::vector<int> > &projection_map, std::vector<int> &projection_imap_var, adj_type &projection_adjacency_list, REAL* &array_of_sij) {

  //for(int node=0; node < num_vars; node++)
  //  delete []projection_map[node];
  //delete []projection_map;
  
  //delete projection_imap_var;
  
  //delete []projection_adjacency_list;

  delete []array_of_sij;
}

// Compute a single edge weight in the projection graph
template<typename MRF_CONSTRUCTOR>
REAL 
Cycle<MRF_CONSTRUCTOR>::find_smn(
      const std::vector<int>& partition_i, 
      const std::vector<int>& partition_j, 
      const int factorId)
{
   const int var_i = std::get<0>(gm_.GetPairwiseVariables(factorId));
   const int var_j = std::get<1>(gm_.GetPairwiseVariables(factorId));
   const int var_i_size = gm_.GetNumberOfLabels(var_i);
   const int var_j_size = gm_.GetNumberOfLabels(var_j);

  int* whole_i = new int[var_i_size];
  int* whole_j = new int[var_j_size];  
  for (int i = 0; i < var_i_size; i++) {
    whole_i[i] = 0;
  }
  for (int i = 0; i < var_j_size; i++) {
    whole_j[i] = 0;
  }

  REAL smn = -Inf;
  for (int i = 0; i < partition_i.size(); i++) {
    whole_i[partition_i[i]] = 1;
    for (int j = 0; j < partition_j.size(); j++) {
      whole_j[partition_j[j]] = 1;
    }
  }
  
  REAL sec_max = -Inf;
  // do zrobienia: or exchange the two loops?
  for (int i = 0; i < var_i_size; i++) {    
     for (int j = 0; j < var_j_size; j++) {
      if (whole_i[i] == whole_j[j]) {
        std::array<int,2> inds = {{i,j}};// inds.push_back(i); inds.push_back(j);
        //REAL temp_val = gm_[factorId](inds.begin());//edge_belief->GetVal(inds);
         REAL temp_val = -gm_.GetPairwiseValue(factorId,inds[0],inds[1]);

        if (smn < temp_val) smn = temp_val;
      }
      else {
        std::array<int,2> inds = {{i,j}};// inds.push_back(i); inds.push_back(j);
        //REAL temp_val = gm_[factorId](inds.begin());//edge_belief->GetVal(inds);
         REAL temp_val = -gm_.GetPairwiseValue(factorId,inds[0],inds[1]);
        if (sec_max < temp_val) sec_max = temp_val;
      }
    }
  }
  smn -= sec_max;

  delete [] whole_i;
  delete [] whole_j;

  return smn;
}

// Compute a single edge weight in the projection graph (more efficiently)
template<typename MRF_CONSTRUCTOR>
REAL 
Cycle<MRF_CONSTRUCTOR>::find_smn_state_i(
      const SIGNED_INDEX single_i,
      const std::vector<int>& partition_j,
      const int factorId,
      const marray::Marray<REAL>& max_i_bij_not_xi) 
      //int single_i,
      //int var_i_size,
      //int var_j_size,
      //MulDimArr* edge_belief,
      //std::vector<std::vector<REAL> >& max_i_bij_not_xi) 
{
	REAL max, sec_max;
	max = sec_max = -Inf;
   const int var_i = std::get<0>(gm_.GetPairwiseVariables(factorId));
   const int var_j = std::get<1>(gm_.GetPairwiseVariables(factorId));
   const int var_i_size = gm_.GetNumberOfLabels(var_i);
   const int var_j_size = gm_.GetNumberOfLabels(var_j);
   //std::vector<REAL> whole_i(var_i_size,0.0);
   //std::vector<REAL> whole_j(var_j_size,0.0);
   // do zrobienia: intialize more efficiently
   std::vector<REAL> whole_i(var_i_size);
   std::vector<REAL> whole_j(var_j_size);
	for (int i = 0; i < var_i_size; i++) {
		whole_i[i] = 0;
	}
	for (int j = 0; j < var_j_size; j++) {
		whole_j[j] = 0;
	}
	whole_i[single_i] = 1;
	for (int j = 0; j < partition_j.size(); j++) { 
		//brute force max_{pi(x_i)=pi(x_j)=1}bij(x_i,x_j)
		whole_j[partition_j[j]] = 1;
		std::array<int,2> inds = {{single_i,partition_j[j]}};// inds.push_back(single_i); inds.push_back(partition_j[j]);
		//REAL tmp = gm_[factorId](inds.begin());//edge_belief->GetVal(inds);
		REAL tmp = -gm_.GetPairwiseValue(factorId,inds[0],inds[1]);
		if (max < tmp) max = tmp;
	}

	for (int j = 0; j < var_j_size; j++) {
		if (whole_j[j] == 0) {
			//bruteforce max_{pi(x_i)=pi(x_j)=0}bij(x_i,x_j)
			//std::vector<int> inds; inds.push_back(single_i); inds.push_back(j);
         std::array<int,2> inds = {{single_i,j}};
			//REAL tmp = gm_[factorId](inds.begin());//edge_belief->GetVal(inds);
         REAL tmp = -gm_.GetPairwiseValue(factorId,inds[0],inds[1]);
			if (sec_max < tmp) sec_max = tmp; 

			if (max < max_i_bij_not_xi(single_i,j)) max = max_i_bij_not_xi(single_i,j);
		}
		else {
			if (sec_max < max_i_bij_not_xi(single_i,j)) sec_max = max_i_bij_not_xi(single_i,j);			
		}
	}	
	
	return max - sec_max;
}

// Compute a single edge weight in the projection graph (more efficiently)
template<typename MRF_CONSTRUCTOR>
REAL 
Cycle<MRF_CONSTRUCTOR>::find_smn_state_j(
      const SIGNED_INDEX single_j,
      const std::vector<int>& partition_i,
      const int factorId,
      const marray::Marray<REAL>& max_j_bij_not_xj) 
      //const std::vector<int>& partition_i, int var_i_size, int single_j, int var_j_size, MulDimArr* edge_belief, std::vector<std::vector<REAL> >& max_j_bij_not_xj) 
{
	REAL max, sec_max;
	max = sec_max = -Inf;
   const int var_i = std::get<0>(gm_.GetPairwiseVariables(factorId));
   const int var_j = std::get<1>(gm_.GetPairwiseVariables(factorId));
   const int var_i_size = gm_.GetNumberOfLabels(var_i);
   const int var_j_size = gm_.GetNumberOfLabels(var_j);
   //std::vector<REAL> whole_i(var_i_size,0.0);
   //std::vector<REAL> whole_j(var_j_size,0.0);
   std::vector<REAL> whole_i(var_i_size);
   std::vector<REAL> whole_j(var_j_size);
	for (int i = 0; i < var_i_size; i++) {
		whole_i[i] = 0;
	}
	for (int j = 0; j < var_j_size; j++) {
		whole_j[j] = 0;
	}
	whole_j[single_j] = 1;
	for (int i = 0; i < partition_i.size(); i++) { 
		//brute force max_{pi(x_i)=pi(x_j)=1}bij(x_i,x_j)
		whole_i[partition_i[i]] = 1;
		//std::vector<int> inds; inds.push_back(partition_i[i]); inds.push_back(single_j);
		std::array<int,2> inds = {{partition_i[i], single_j}};
		//REAL tmp = gm_[factorId](inds.begin());//edge_belief->GetVal(inds);
		REAL tmp = -gm_.GetPairwiseValue(factorId,inds[0],inds[1]);
		if (max < tmp) max = tmp;
	}

	for (int i = 0; i < var_i_size; i++) {
		if (whole_i[i] == 0) {
			//bruteforce max_{pi(x_i)=pi(x_j)=0}bij(x_i,x_j)
			//std::vector<int> inds; inds.push_back(i); inds.push_back(single_j);
			std::array<int,2> inds = {{i, single_j}};
         //REAL tmp = gm_[factorId](inds.begin());//edge_belief->GetVal(inds);
         REAL tmp = -gm_.GetPairwiseValue(factorId,inds[0],inds[1]);
			if (sec_max < tmp) sec_max = tmp; 

			if (max < max_j_bij_not_xj(i,single_j)) max = max_j_bij_not_xj(i,single_j);
		}
		else {
			if (sec_max < max_j_bij_not_xj(i,single_j)) sec_max = max_j_bij_not_xj(i,single_j);			
		}
	}	

	return max - sec_max;
}


// This will find the best partitioning of the states of a single edge
// Implements the FindPartition algorithm described in supplementary material of UAI 2012 paper.
template<typename MRF_CONSTRUCTOR>
void 
Cycle<MRF_CONSTRUCTOR>::find_partition(
      std::vector<std::map<std::vector<int>, int> >& partition_set, 
      const int factorId,
      int index_i,  // do zrobienia: index_i and index_j not needed
      int index_j) 
{
   const int var_i = std::get<0>(gm_.GetPairwiseVariables(factorId));
   const int var_j = std::get<1>(gm_.GetPairwiseVariables(factorId));
   int size_i = gm_.GetNumberOfLabels(var_i);
   int size_j = gm_.GetNumberOfLabels(var_j);

  //int size_i = gm_[factorId].shape(0);//mplp.m_var_sizes[index_i];
  //int size_j = gm_[factorId].shape(1);//mplp.m_var_sizes[index_j];

  std::vector<Node*> x_i; //partition for node i
  std::vector<Node*> x_j; //partition for node j

  for (int i = 0; i < size_i; i++) {
    Node* t = new Node(2);
    x_i.push_back(t);
  }
  for (int j = 0; j < size_j; j++) {
    Node* t = new Node(1);
    x_j.push_back(t);
  }

  // sort edges by their weights in decreasing order
  std::vector<std::pair<std::pair<int, int>, REAL> > sorted_edge_set;
  for (int i = 0; i < size_i; i++) {
    for (int j = 0; j < size_j; j++) {
      REAL temp_val = -gm_.GetPairwiseValue(factorId,i,j);
      sorted_edge_set.push_back(std::pair<std::pair<int, int>, REAL>(std::make_pair(std::make_pair(i, j), temp_val)));
    }
  }
  std::sort(sorted_edge_set.begin(), sorted_edge_set.end(), edge_sorting3);

  typename std::vector<std::pair<std::pair<int, int>, REAL> >::iterator it = sorted_edge_set.begin();
  REAL val_s = it->second; //val_s = bij(x_i^*, x_j^*);

  // union find step;
  for ( ; it != sorted_edge_set.end(); it++) {
    int ind_i, ind_j;
    ind_i = it->first.first; //state of x_i
    ind_j = it->first.second; //staet of x_j
    Node* i_root = find(x_i[ind_i]); //find head of ind_i
    Node* j_root = find(x_j[ind_j]); //find head of ind_j   
    if (i_root == j_root) continue; // if i, j belong to the same partition already, nothing to do
    
    int bit_i = i_root->bit; //bit of ind_i node
    int bit_j = j_root->bit; //bit of ind_j node   
    if (bit_i == 2 && bit_j == 1) { //two singletons
      ;
    }
    if (bit_i == 2 && bit_j == 3) { 
      if (size_i == 2) { //if number of partition is 2, then break
        val_s -= it->second; //val_s = bij(x_i^*, x_j^*) - max_{pi(x_i) != pi(x_j)} bij(x_i, x_j)
        break;
      }
      size_i--;
    }
    if (bit_i == 3 && bit_j == 1) {
      if (size_j == 2) {
        val_s -= it->second;        
        break;
      }
      size_j--;
    }
    if (bit_i == 3 && bit_j == 3) {
      if (size_i == 2 || size_j == 2) {
        val_s -= it->second;
        break;
      }
      size_i--; size_j--;
    }    
    Node* new_root = merge(i_root, j_root); //if # of partiton > 2, merge two two partitions
    new_root->bit = 3;
  }

  if (size_i != 2 && size_j != 2) {
     std::cout << "something is wrong" << std::endl;
  }
  
  if (size_i == 1 || size_j == 1) {
     std::cout << "should not happen" << std::endl;
    return;
  }
  std::vector<int> part_i, part_j; 

  //find the smallest partition and use it as an index
  // If there's a tie, use the partition with the smallest index. Keep it sorted.
  //int min = gm_[factorId].shape(0) + 1;//mplp.m_var_sizes[index_i] + 1; 
  int min = size_i + 1;//mplp.m_var_sizes[index_i] + 1; 
  Node* head = NULL;
  for (int i = 0; i < size_i; i++) {
    Node* t = find(x_i[i]);
    if ((size_i == 2 || t->bit == 3) && min > t->i_size) {
      min = t->i_size;
      head = t;
    }
  }
  for (int i = 0; i < size_i; i++) {
    Node* t = find(x_i[i]);
    if (head == t)
      part_i.push_back(i);
  }

  //min = gm_[factorId].shape(1) + 1;//mplp.m_var_sizes[index_j] + 1;
  min = size_j + 1;//mplp.m_var_sizes[index_j] + 1;
  for (int j = 0; j < size_j; j++) {
    Node* t = find(x_j[j]);
    if ((size_j == 2 || t->bit == 3) && min > t->j_size) {
      min = t->j_size;
      head = t;
    }
  }
  for (int j = 0; j < size_j; j++) {
    Node* t = find(x_j[j]);
    if (head == t)
      part_j.push_back(j);
  } 

  // We have one std::map per variable.
  
  if (val_s != 0) {
    // Don't add a partition that already exists
    if (partition_set[index_i].find(part_i) == partition_set[index_i].end()) {

      // Add this partition
      // Keep track of the number of times that this partition is used
      partition_set[index_i][part_i] = 1;
    }
    else
      // Keep track of the number of times that this partition is used
      partition_set[index_i][part_i]++;

    if (partition_set[index_j].find(part_j) == partition_set[index_j].end()) {

      // Add this partition. Keep track of the number of times that this partition is used
      partition_set[index_j][part_j] = 1;
    }
    else
      partition_set[index_j][part_j]++;
  }

  for (int i = 0; i < size_i; i++) {
    delete [] x_i[i];
  }
  for (int j = 0; j < size_j; j++) {
    delete [] x_j[j];
  }
}

template<typename MRF_CONSTRUCTOR>
std::tuple<marray::Marray<REAL>, marray::Marray<REAL>, marray::Marray<REAL>>
Cycle<MRF_CONSTRUCTOR>::calculate_cond_min_beliefs(const INDEX factor_id)
{
   const int i = std::get<0>(gm_.GetPairwiseVariables(factor_id));
   const int j = std::get<1>(gm_.GetPairwiseVariables(factor_id));

   // Do some pre-processing for speed.
   SIGNED_INDEX shape[] = {SIGNED_INDEX(gm_.GetNumberOfLabels(i)), SIGNED_INDEX(gm_.GetNumberOfLabels(j))};
   // do zrobienia: rename min and change signs etc.
   marray::Marray<REAL> max_j_bij_not_xj(shape, shape+2);
   marray::Marray<REAL> max_i_bij_not_xi(shape, shape+2);

   for(int state1=0; state1 < gm_.GetNumberOfLabels(i); state1++) {

      // Find max over state2
      REAL largest_val = -gm_.GetPairwiseValue(factor_id,state1,0);
      int largest_ind = 0;
      for(int state2=1; state2 < gm_.GetNumberOfLabels(j); state2++) {
         REAL tmp_val = -gm_.GetPairwiseValue(factor_id,state1,state2);
         if(tmp_val > largest_val) {
            largest_val = tmp_val;
            largest_ind = state2;
         }
      }

      // Find second largest val over state2
      REAL sec_largest_val = -gm_.GetPairwiseValue(factor_id,state1,0);
      for(int state2=1; state2 < gm_.GetNumberOfLabels(j); state2++) {

         REAL tmp_val = -gm_.GetPairwiseValue(factor_id,state1,state2);

         if(tmp_val > sec_largest_val && state2 != largest_ind) {
            sec_largest_val = tmp_val;
         }
      }

      // Assign values
      for(int state2=0; state2 < gm_.GetNumberOfLabels(j); state2++)
         max_j_bij_not_xj(state1,state2) = largest_val;
      max_j_bij_not_xj(state1,largest_ind) = sec_largest_val;
   }


   for(int state1=0; state1 < gm_.GetNumberOfLabels(j); state1++) {

      // Find max over state2
      REAL largest_val = -gm_.GetPairwiseValue(factor_id,0,state1);
      int largest_ind = 0;
      for(int state2=1; state2 < gm_.GetNumberOfLabels(i); state2++) {

         REAL tmp_val = -gm_.GetPairwiseValue(factor_id,state2,state1);//gm_[factorId](inds.begin());//edge_belief->GetVal(inds);

         if(tmp_val > largest_val) {
            largest_val = tmp_val;
            largest_ind = state2;
         }
      }

      // Find second largest val over state2
      REAL sec_largest_val = -gm_.GetPairwiseValue(factor_id,0,state1);
      for(int state2=1; state2 < gm_.GetNumberOfLabels(i); state2++) {

         REAL tmp_val = -gm_.GetPairwiseValue(factor_id,state2,state1);//gm_[factorId](inds.begin());//edge_belief->GetVal(inds);

         if(tmp_val > sec_largest_val && state2 != largest_ind) {
            sec_largest_val = tmp_val;
         }
      }

      // Assign values
      for(int state2=0; state2 < gm_.GetNumberOfLabels(i); state2++)
         max_i_bij_not_xi(state2,state1) = largest_val;
      max_i_bij_not_xi(largest_ind,state1) = sec_largest_val;
   }

   // decompose: max_{x_j!=x_j}max_{x_i != x_i'}. Then, use the above computations.
   marray::Marray<REAL> max_ij_bij_not_xi_xj(shape, shape+2);
   for(int state1=0; state1 < gm_.GetNumberOfLabels(i); state1++) {

      // Find max over state2
      REAL largest_val = -std::numeric_limits<REAL>::infinity();
      int largest_ind = 0;
      for(int state2=1; state2 < gm_.GetNumberOfLabels(j); state2++) {
         REAL tmp_val = max_i_bij_not_xi(state1,state2);

         if(tmp_val > largest_val) {
            largest_val = tmp_val;
            largest_ind = state2;
         }
      }

      // Find second largest val over state2
      REAL sec_largest_val = -std::numeric_limits<REAL>::infinity();
      for(int state2=0; state2 < gm_.GetNumberOfLabels(j); state2++) {
         REAL tmp_val = max_i_bij_not_xi(state1,state2);

         if(tmp_val > sec_largest_val && state2 != largest_ind) {
            sec_largest_val = tmp_val;
         }
      }

      // Assign values
      for(int state2=0; state2 < gm_.GetNumberOfLabels(j); state2++)
         max_ij_bij_not_xi_xj(state1,state2) = largest_val;
      max_ij_bij_not_xi_xj(state1,largest_ind) = sec_largest_val;
   }

   return std::move(std::make_tuple( std::move(max_j_bij_not_xj), std::move(max_i_bij_not_xi), std::move(max_ij_bij_not_xi_xj)));
}

// Create the expanded projection graph by including all singleton partitions and also
// all partitions found by calling FindPartition on all edges.
template<typename MRF_CONSTRUCTOR>
int 
Cycle<MRF_CONSTRUCTOR>::create_expanded_projection_graph(
      std::vector<int>& projection_imap_var,
      std::vector<std::vector<std::pair<int, REAL> > >& projection_adjacency_list, 
      std::map<std::pair<int, int>, REAL>&projection_edge_weights,
      REAL* &array_of_sij,
      int& array_of_sij_size,
      std::vector<std::vector<int> >& partition_imap) 
{
  // projection_imap_var std::maps from projection node to variable
  // partition_imap std::maps from projection node to std::vector of states

  int num_of_vars = gm_.GetNumberOfVariables();
  std::vector<std::map<std::vector<int>,int> > partition_set;
  std::set<REAL> set_of_sij;
  int num_projection_nodes = 0;

  for (int i = 0; i < num_of_vars; i++) {
    std::map<std::vector<int>,int> p;
    partition_set.push_back(p);

    // push all singleton partitions
    for (int i_state = 0; i_state < gm_.GetNumberOfLabels(i); i_state++) {
      std::vector<int> single;
      single.push_back(i_state);
      partition_set[i][single] = 1;
    }
  }

  bool partition = true;
  clock_t find_partition_start_time;
  if (partition) {
     find_partition_start_time = clock();
     //THIS IS THE NEW PARTITIONING ALGORITHM

     //for (std::mapType::const_iterator it = mplp.m_intersect_std::map.begin(); it != mplp.m_intersect_std::map.end(); ++it) {    
     for(size_t factorId=0; factorId<gm_.GetNumberOfPairwiseFactors(); factorId++) {
        const int i = std::get<0>(gm_.GetPairwiseVariables(factorId));//gm_.variableOfFactor(factorId,0);
        const int j = std::get<1>(gm_.GetPairwiseVariables(factorId));//gm_.variableOfFactor(factorId,1);
        // Check to see if i and j have at least two states each -- otherwise, cannot be part of any frustrated edge
        if(gm_.GetNumberOfLabels(i) <= 1 || gm_.GetNumberOfLabels(j) <= 1)
           continue;

        find_partition(partition_set, factorId, i, j);
     }
  }
  clock_t find_partition_end_time = clock();
  REAL find_partition_total_time = (REAL)(find_partition_end_time - find_partition_start_time)/CLOCKS_PER_SEC;
  if (DEBUG_MODE) {
     std::cout << " -- find_partition. Took " << find_partition_total_time << " seconds" << std::endl; 
  }
  //}

  // Create one projection node for each remaining partition. Overload value of partition_set
  // to now refer to the node number in the projection graph.
  for(int i=0; i < partition_set.size(); i++) {
    for(typename std::map<std::vector<int>, int>::iterator it = partition_set[i].begin(); it != partition_set[i].end(); it++) {
      it->second = num_projection_nodes++;
      projection_imap_var.push_back(i);
      partition_imap.push_back(it->first);

      std::vector<std::pair<int, REAL> > temp;
      projection_adjacency_list.push_back(temp);
    }
  }

  clock_t find_smn_start_time = clock();

  // Create projection graph edges for each edge of original graph and each std::pair of partitions
  for(size_t factorId=0; factorId<gm_.GetNumberOfPairwiseFactors(); factorId++) {
        const int i = std::get<0>(gm_.GetPairwiseVariables(factorId));//gm_.variableOfFactor(factorId,0);
        const int j = std::get<1>(gm_.GetPairwiseVariables(factorId));//gm_.variableOfFactor(factorId,1);
        //int ij_intersect_loc = it->second;

        // Do some pre-processing to make the case of single state partitions very fast

        /*
        MulDimArr* edge_belief = &mplp.m_sum_into_intersects[ij_intersect_loc];
        if (mplp.m_all_intersects[ij_intersect_loc][0] != i) {
           int tmp_i = i;
           i = j;
           j = tmp_i;
        }
        */

        if(gm_.GetNumberOfLabels(i) <= 1 || gm_.GetNumberOfLabels(j) <= 1) continue;

        auto max_bij_cond = calculate_cond_min_beliefs(factorId);
        const auto& max_j_bij_not_xj = std::get<0>(max_bij_cond);
        const auto& max_i_bij_not_xi = std::get<1>(max_bij_cond);
        const auto& max_ij_bij_not_xi_xj = std::get<2>(max_bij_cond);

        // Now, for each partition of node i and each partition of node j, compute edge weights
        // If the edge weight is non-zero, insert edge into adjacency list

        for(typename std::map<std::vector<int>, int>::iterator it_i = partition_set[i].begin(); it_i != partition_set[i].end(); it_i++) {
           int n = it_i->second;

           for(typename std::map<std::vector<int>, int>::iterator it_j = partition_set[j].begin(); it_j != partition_set[j].end(); it_j++) {
              int m = it_j->second;

              REAL smn = 0;
              if (it_i->first.size() == 1 && it_j->first.size() == 1) {
                 int xi = it_i->first[0]; int xj = it_j->first[0];
                 //std::vector<int> inds; inds.push_back(xi); inds.push_back(xj);
                 std::array<int,2> inds = {{xi,xj}};
                 //REAL tmp_val = gm_[factorId](inds.begin());//edge_belief->GetVal(inds);
                 REAL tmp_val = -gm_.GetPairwiseValue(factorId,inds[0],inds[1]);//gm_[factorId](inds.begin());//edge_belief->GetVal(inds);
                 smn = std::max(tmp_val, max_ij_bij_not_xi_xj(xi,xj)) - std::max(max_i_bij_not_xi(xi,xj), max_j_bij_not_xj(xi,xj));
              }
              else if (it_i->first.size() == 1) {
                 int xi = it_i->first[0];
                 smn = find_smn_state_i(xi, it_j->first, factorId, max_i_bij_not_xi);
              }
              else if (it_j->first.size() == 1) {
                 int xj = it_j->first[0];
                 smn = find_smn_state_j(xj, it_i->first, factorId, max_j_bij_not_xj);
              }
              else {
                 // This computes smn
                 smn = find_smn(it_i->first, it_j->first, factorId);
              }

              if (smn != 0) {
                 set_of_sij.insert(std::abs(smn));

                 projection_adjacency_list[m].push_back(std::make_pair(n, smn));
                 projection_adjacency_list[n].push_back(std::make_pair(m, smn));

                 projection_edge_weights.insert(std::pair<std::pair<int, int>, REAL>(std::pair<int, int>(n, m), smn));
                 projection_edge_weights.insert(std::pair<std::pair<int, int>, REAL>(std::pair<int, int>(m, n), smn));
              }

           }
        }
  }
  clock_t find_smn_end_time = clock();
  REAL find_smn_total_time = (REAL)(find_smn_end_time - find_smn_start_time)/CLOCKS_PER_SEC;
  if (DEBUG_MODE) {
     std::cout << " -- find_smn. Took " <<  find_smn_total_time << " seconds" << std::endl;
  }
  
  array_of_sij = new REAL[set_of_sij.size()];
  array_of_sij_size = 0;
  
  for (typename std::set<REAL>::iterator set_iter = set_of_sij.begin(); set_iter != set_of_sij.end(); set_iter++) {  //NOTE: PASSED AS FUNCTION ARG, DO NOT DELETE HERE!!
    array_of_sij[array_of_sij_size++] = *set_iter;
  }

  return num_projection_nodes;
}


// Creates the k-projection graph (just a single partition per variable)
// only include variables which are not locally excluded
template<typename MRF_CONSTRUCTOR>
void 
Cycle<MRF_CONSTRUCTOR>::create_k_projection_graph(
      std::vector<std::vector<int> > &projection_map, 
      int& num_projection_nodes, 
      std::vector<int> &projection_imap_var, 
      std::vector<std::vector<int> > &partition_imap,
      std::map<std::pair<int, int>, REAL>& projection_edge_weights,
      adj_type &projection_adjacency_list,
      REAL* &array_of_sij,
      int& array_of_sij_size) 
{

  // TODO: make sure for binary variables that there is only one node per variable (rather than 2).
  // TODO: projection_edge_weights can likely be removed from this function and elsewhere.
  // TODO: most of these ints can be changed to be unsigned and/or fewer bits. Look into memory allocation.
  
  // only look at locally optimal variables
   //std::vector<std::vector<bool> > locallyOptimal = LocallyMinimal(gm_); // do zrobienia

  // Initialize the projection graph
  int projection_node_iter = 0;
  for(int node=0; node < gm_.GetNumberOfVariables(); node++) {
    std::vector<int> i;
    
    for(int state=0; state < gm_.GetNumberOfLabels(node); state++) {
       std::vector<std::pair<int, REAL> > temp;
       projection_adjacency_list.push_back(temp);
       i.push_back(projection_node_iter++);
    }
    projection_map.push_back(i);
  }
  // Construct inverse std::map
  num_projection_nodes = projection_node_iter;
			
  projection_node_iter = 0;
  for(int node=0; node < gm_.GetNumberOfVariables(); node++) {
    for(int state=0; state < gm_.GetNumberOfLabels(node); state++) {
       projection_imap_var.push_back(node);

       std::vector<int> tmp_state_vector;
       tmp_state_vector.push_back(state);
       partition_imap.push_back(tmp_state_vector);

       projection_node_iter++;
    }
  }
  
  // Iterate over all of the edges (we do this by looking at the edge intersection sets)
  std::set<REAL> set_of_sij;
  for(size_t factorId=0; factorId<gm_.GetNumberOfPairwiseFactors(); factorId++) {
        // Get the two nodes i & j and the edge intersection set. Put in right order.
        const int i = std::get<0>(gm_.GetPairwiseVariables(factorId));// gm_.variableOfFactor(factorId,0);
        const int j = std::get<1>(gm_.GetPairwiseVariables(factorId));//gm_.variableOfFactor(factorId,1);
        //int ij_intersect_loc = it->second;

        /*
        MulDimArr* edge_belief = &mplp.m_sum_into_intersects[ij_intersect_loc];
        if(mplp.m_all_intersects[ij_intersect_loc][0] != i)  { // swap i and j
           int tmp_i = i;
           i = j;
           j = tmp_i;
        }
        */

        // Check to see if i and j have at least two states each -- otherwise, cannot be part of any frustrated edge
        if(gm_.GetNumberOfLabels(i) <= 1 || gm_.GetNumberOfLabels(j) <= 1)
           continue;

        auto max_bij_cond = calculate_cond_min_beliefs(factorId);
        const auto& max_j_bij_not_xj = std::get<0>(max_bij_cond);
        const auto& max_i_bij_not_xi = std::get<1>(max_bij_cond);
        const auto& max_ij_bij_not_xi_xj = std::get<2>(max_bij_cond);

        // For each of their states
        for(int xi=0; xi < gm_.GetNumberOfLabels(i); xi++) {
           int m = projection_map[i][xi];

           for(int xj=0; xj < gm_.GetNumberOfLabels(j); xj++) {
              const int n = projection_map[j][xj];

              //std::vector<int> inds; inds.push_back(xi); inds.push_back(xj);
              const REAL tmp_val = -gm_.GetPairwiseValue(factorId,xi,xj);//gm_[factorId](inds.begin());//edge_belief->GetVal(inds);

              // Compute s_mn for this edge
              REAL val_s = std::max(tmp_val, max_ij_bij_not_xi_xj(xi,xj)) - std::max(max_i_bij_not_xi(xi,xj), max_j_bij_not_xj(xi,xj));
              //if(std::max(tmp_val, max_ij_bij_not_xi_xj(xi,xj)) == std::numeric_limits<REAL>::infinity()
              //      && std::max(max_i_bij_not_xi(xi,xj), max_j_bij_not_xj(xi,xj)) == std::numeric_limits<REAL>::infinity()) {
              //      val_s = 0.0;
              //      }

              //assert(std::isnan(val_s) == false);

              // TODO: use threshold here, to make next stage faster
              if(val_s != 0 && !std::isnan(val_s)) {          
                 projection_adjacency_list[m].push_back(std::make_pair(n, val_s));
                 projection_adjacency_list[n].push_back(std::make_pair(m, val_s));
                 set_of_sij.insert(std::abs(val_s));
              }

              // Insert into edge weight std::map
              projection_edge_weights.insert(std::pair<std::pair<int,int>,REAL>(std::pair<int,int>(n, m), val_s));
              projection_edge_weights.insert(std::pair<std::pair<int,int>,REAL>(std::pair<int,int>(m, n), val_s));
           }
        }

  }
  
  // Sort list_of_sij and remove duplicates
  array_of_sij = new REAL[set_of_sij.size()]; 
  array_of_sij_size = 0;
	
  for(typename std::set<REAL>::iterator set_iter = set_of_sij.begin(); set_iter != set_of_sij.end(); ++set_iter) {
    array_of_sij[array_of_sij_size++] = *set_iter;
  }
}


// Does binary search over the edge weights to find largest edge weight
// such that there is an odd-signed cycle.
template<typename MRF_CONSTRUCTOR>
REAL 
Cycle<MRF_CONSTRUCTOR>::find_optimal_R(adj_type &projection_adjacency_list, REAL* &array_of_sij, int array_of_sij_size) {

  // Do binary search over sij
  int bin_search_lower_bound = 0;
  int bin_search_upper_bound = array_of_sij_size;
  REAL sij_min = -1;
  int num_projection_nodes = projection_adjacency_list.size();

  // do zrobienia: czy to prawidlowo?
  while(bin_search_lower_bound <= bin_search_upper_bound) {
				
    // Compute mid-point
    int R_pos = floor((bin_search_lower_bound + bin_search_upper_bound)/2);
    REAL R = array_of_sij[R_pos];
    std::cout << "Test: " << bin_search_lower_bound << " < " << R_pos << " < " << bin_search_upper_bound << ", R = " << R << std::endl;
				
    // Does there exist an odd signed cycle using just edges with sij >= R? If so, go up. If not, go down.
    bool found_odd_signed_cycle = false;
				
    // Initialize
    // TODO: do not allocate memory for this every time! Just re-initialize.
    
    std::vector<int> node_sign(num_projection_nodes);
    std::fill(node_sign.begin(), node_sign.end(), 0); // denotes "not yet seen"      
    //for(int m=0; m<num_projection_nodes; m++) {
    //  node_sign[m] = 0; // denotes "not yet seen"      
    //}
    
    // Graph may be disconnected, so check from all nodes    
    for(int i = 0; i < num_projection_nodes && !found_odd_signed_cycle; i++) {
      if(node_sign[i] == 0) {
        node_sign[i] = 1; // root node
        std::queue<int> q;
        q.push(i);
        
        while (!q.empty() && !found_odd_signed_cycle) {
          const int current = q.front();
          assert(node_sign[current] != 0); // do zrobienia: correct?
          q.pop();
          
          for (int j = 0; j < projection_adjacency_list[current].size(); j++) {
            const REAL smn = projection_adjacency_list[current][j].second;
            // Ignore edges with weight less than R
            if (std::abs(smn) < R)
              continue;
            
            const int next = projection_adjacency_list[current][j].first;
            const int sign_of_smn = (smn > 0) - (smn < 0);
            
            if (node_sign[next] == 0) {
              node_sign[next] = node_sign[current] * sign_of_smn;
              assert(node_sign[next] != 0);
              q.push(next);
            }
            else if(node_sign[next] == -node_sign[current] * sign_of_smn) {
              // Found an odd-signed cycle! Can quit.
              found_odd_signed_cycle = true;
              break;
            }
          }
        }
      }
    }

    if(found_odd_signed_cycle) {
      sij_min = R;
      bin_search_lower_bound = R_pos+1;
    }
    else
      bin_search_upper_bound = R_pos-1;
  }
  //assert(sij_min >= 0);

  return sij_min;
}


// Returns an array which is a random permutation of the numbers 0 through n-1.
template<typename MRF_CONSTRUCTOR>
int* 
Cycle<MRF_CONSTRUCTOR>::random_permutation(int n) {
  int *p = new int[n];
  for (int i = 0; i < n; ++i) {
    int j = rand() % (i + 1);
    p[i] = p[j];
    p[j] = i;
  }
  return p;
}


// Given an undirected graph, finds an odd-signed cycle.
// This works by breath-first search. Better might be to find minimal depth tree.
// NOTE: this function depends on the random seed becase it creates a random spanning tree.
template<typename MRF_CONSTRUCTOR>
void 
Cycle<MRF_CONSTRUCTOR>::FindCycles(std::vector<std::list<int> > &cycle_set, const REAL optimal_R, const int ncycles_to_add, adj_type &projection_adjacency_list) {

  REAL R = optimal_R;
  int num_projection_nodes = projection_adjacency_list.size();
  
  // Initialize  (NOTE: uses heap allocation)
  int *node_sign = new int[num_projection_nodes], *node_parent = new int[num_projection_nodes], *node_depth = new int[num_projection_nodes];
  
  for(int i = 0; i < num_projection_nodes; i++) {
    node_sign[i] = 0; // denotes "not yet seen"
    node_depth[i] = -1;
    node_parent[i] = -1;
  }

  // construct the rooted spanning tree(s) -- randomly choose a root!
  int* random_projection_node = random_permutation(num_projection_nodes);
  //  if(DEBUG_MODE) {
  //    cout << "First random node is: " << random_projection_node[0] << std::endl;
  //  }
  
  for(int ri = 0; ri < num_projection_nodes; ri++) {
    int i = random_projection_node[ri];    
    if(node_sign[i] == 0) {
      node_sign[i] = 1; // root node
      node_parent[i] = i;
      node_depth[i] = 0;
      
      std::queue<int> q;
      q.push(i);

      while (!q.empty()) {
        int current = q.front();        
        q.pop();
        
        // Randomize the adjacency list        
        int* random_nbrs = random_permutation(projection_adjacency_list[current].size());
        for (int rj = 0; rj < projection_adjacency_list[current].size(); rj++) {
          int j = random_nbrs[rj];
          int next = projection_adjacency_list[current][j].first;
          if (node_sign[next] == 0) {
            REAL smn = projection_adjacency_list[current][j].second;
            if (std::abs(smn) < R) {
              ;
            }
            else {
              int sign_of_smn = (smn > 0) - (smn < 0);
              node_sign[next] = node_sign[current] * sign_of_smn;
              node_parent[next] = current;									
              node_depth[next] = node_depth[current] + 1;
              q.push(next);
            }
          }
        }
        delete [] random_nbrs;
      }
    }					
  }
  delete [] random_projection_node;

  // construct edge set that contains edges that are not parts of the tree
  // TODO: This can be made faster! Can be performed in CONSTANT time.
  std::vector<std::list<int> > edge_set;
  std::map<std::pair<int, int>, int> edge_map;
  
  for (int i = 0; i < num_projection_nodes; i++) {
    for (int j = 0; j < projection_adjacency_list[i].size(); j++) {
      if (node_parent[i] == j || node_parent[j] == i) {							
        continue;
      }
      REAL smn = projection_adjacency_list[i][j].second;

      if (std::abs(smn) < R) {
        continue;
      }
      int jj = projection_adjacency_list[i][j].first;
      int sign_of_smn = (smn > 0) - (smn < 0);

      if (node_sign[i] == -node_sign[jj] * sign_of_smn) { //cycle found
        int depth_of_i, depth_of_j, temp_di, temp_dj;
        depth_of_i = temp_di = node_depth[i];
        depth_of_j = temp_dj = node_depth[jj];
        
        int anc_i, anc_j;
        anc_i = i;
        anc_j = jj;

        while (temp_dj > temp_di) {
          anc_j = node_parent[anc_j];
          temp_dj--;
        }
        while (temp_di > temp_dj) {
          anc_i = node_parent[anc_i];
          temp_di--;
        }
        while (temp_di >= 0) {
          if (anc_i == anc_j) { //least common ancestor found
             std::list<int> temp;
            temp.push_back(i);
            temp.push_back(jj);
            temp.push_back(depth_of_i - temp_di + depth_of_j - temp_dj);
            assert( depth_of_i - temp_di + depth_of_j - temp_dj >= 0 );
            typename std::map<std::pair<int, int>, int>::iterator m_it = edge_map.find(std::make_pair(i, jj));
            if (m_it == edge_map.end()) {
              edge_map[std::pair<int, int>(i, jj)] = int(depth_of_i - temp_di + depth_of_j - temp_dj);
              edge_set.push_back(temp);		
            }								
            break;
          }
          anc_i = node_parent[anc_i];
          anc_j = node_parent[anc_j];		
          temp_di--;
          temp_dj--;								
        }
      }
    }
  }
  
  // sort the edge by the distance between two nodes in the tree
  std::sort(edge_set.begin(), edge_set.end(), edge_sorting);	

  // Note: here we do not check for duplicate variables 
  for (int i = 0; i < edge_set.size() && cycle_set.size() < ncycles_to_add; i++) {

    // Get next edge (and depth) from sorted list
    std::list<int> temp = edge_set.back();
    edge_set.pop_back();

    // Find least common ancestor
    int left, right, left_d, right_d, ancestor, left_anc, right_anc;
    left = temp.front();
    temp.pop_front();
    right = temp.front();
    left_d = node_depth[left];
    right_d = node_depth[right];
    left_anc = left;
    right_anc = right;

    // find the first common ancestor
    while (left_d > right_d) {
      left_anc = node_parent[left_anc];
      left_d--;
    }
    while (right_d > left_d) {
      right_anc = node_parent[right_anc];
      right_d--;
    }
    while (left_anc != right_anc) {
      left_anc = node_parent[left_anc];
      right_anc = node_parent[right_anc];
    }
    ancestor = left_anc;

    std::list<int> cycle;

    // backtrace the found cycle
    while (node_parent[left] != ancestor) {
      cycle.push_back(left);
      left = node_parent[left];
    }
    cycle.push_back(left);
    cycle.push_back(node_parent[left]);
    while (node_parent[right] != ancestor) {
      cycle.push_front(right);
      right = node_parent[right];
    }
    cycle.push_front(right);
    cycle_set.push_back(cycle);
  }
  delete [] node_sign;
  delete [] node_parent;
  delete [] node_depth;
}


// Check to see if there are any duplicate variables, and, if so, shortcut
template<typename MRF_CONSTRUCTOR>
void 
Cycle<MRF_CONSTRUCTOR>::shortcut(
      std::list<int> &cycle,
      std::vector<int> &projection_imap_var,
      std::map<std::pair<int, int>, REAL>& projection_edge_weights,
      int& num_projection_nodes) 
{
   assert(cycle.size() >= 3);
  bool exist_duplicates = true;
  std::vector<int> cycle_array(cycle.size()); int tmp_ind = 0;
  std::vector<int> cycle_sign(cycle.size());

  for (typename std::list<int>::iterator it=cycle.begin(); it!=cycle.end(); ++it) {
    cycle_array[tmp_ind++] = *it;
  }

  int cycle_start = 0;
  int cycle_end = cycle.size()-1;
				
  // first, shorten the cycle
  while(exist_duplicates) {
    exist_duplicates = false;
    
    // This std::map allows us to quickly find duplicates and their locations within the cycle
    std::map<int, int> duplicates_map;
					
    // Must keep track of sign.
    cycle_sign[cycle_start] = 1;
    
    duplicates_map.insert(std::pair<int,int>(projection_imap_var[cycle_array[cycle_start]], cycle_start));
    for(int i=cycle_start+1; i <= cycle_end; i++) {
      
      // Get edge weight
      typename std::map< std::pair<int,int>, REAL>::iterator weights_iter = projection_edge_weights.find(std::make_pair(cycle_array[i-1], cycle_array[i]));
      
      REAL smn = weights_iter->second;
      int sign_of_smn = (smn > 0) - (smn < 0);
      cycle_sign[i] = cycle_sign[i-1]*sign_of_smn;
			
      // Does node i already exist in the std::map?
      typename std::map<int, int>::iterator dups_iter = duplicates_map.find(projection_imap_var[cycle_array[i]]);
      if( dups_iter != duplicates_map.end() ) { // Duplicate found
        int first_occurence_index = dups_iter->second;
        
        // Look to see if the first half of the cycle is violated.
        weights_iter = projection_edge_weights.find(std::make_pair(cycle_array[first_occurence_index], cycle_array[i-1]));        
        REAL edge_weight = weights_iter->second;
        int sign_of_edge = (edge_weight > 0) - (edge_weight < 0);
        
        // NOTE: Bound decrease could be LESS than promised!
        if(cycle_sign[first_occurence_index] == -cycle_sign[i-1] * sign_of_edge) {
          
          // Found a violated cycle! Can quit
          cycle_start = first_occurence_index;
          cycle_end = i-1;
        }
        else {
          
          // Otherwise, cut out the part of the cycle from first_occurence_index to cycle_end, and continue
          for(int pos=first_occurence_index-1; pos>= cycle_start; pos--)
            cycle_array[pos+i-first_occurence_index] = cycle_array[pos];
          cycle_start += i-first_occurence_index;
          
          exist_duplicates = true;
        }
        
        break;
      }
      else {
        // Insert into the std::map
        duplicates_map.insert(std::pair<int,int>(projection_imap_var[cycle_array[i]], i));
      }
    }
  }

  // Modify the cycle
  cycle.clear();
  for(int i=cycle_start; i <= cycle_end; i++)
    cycle.push_back(cycle_array[i]);
}


// Triangulate cycles and add sparse triplets to the graphical model
template<typename MRF_CONSTRUCTOR>
template<typename ADD_TRIPLET_FUNCTION>
int
Cycle<MRF_CONSTRUCTOR>::add_cycle(
      //std::function<void(const int,const int, const int, const std::vector<SIGNED_INDEX>, const std::vector<SIGNED_INDEX>, const std::vector<SIGNED_INDEX>)> addTripletFun,
      ADD_TRIPLET_FUNCTION addTripletFun,
      std::list<int> &cycle, 
      std::vector<int> &projection_imap_var, 
      std::vector<std::vector<int> >& partition_imap,
      std::map<std::vector<int>, bool >& triplet_set, 
      int& num_projection_nodes) 
{
   assert(cycle.size() >= 3);
  int nClustersAdded = 0;

  // Number of clusters we're adding is length_cycle - 2
  int nNewClusters = cycle.size() - 2;
  std::vector<TripletCluster> newCluster(nNewClusters);
  int cluster_index = 0;

  // do zrobienia: partitioning can be done less redundant
  // construct partition of label space for cycle
  std::vector< std::vector<SIGNED_INDEX> > pi(0);
  for (typename std::list<int>::iterator it=cycle.begin(); it!=cycle.end(); ++it) {
     std::vector<int> temp = partition_imap[*it];          
     pi.push_back( std::vector<SIGNED_INDEX>(gm_.GetNumberOfLabels( projection_imap_var[*it] ), 0) );
     for(size_t l=0; l<temp.size(); l++) {
        pi.back()[ temp[l] ] = 1;
     }
  }



  // Convert cycle to array
  std::vector<int> cycle_array(cycle.size()); int tmp_ind = 0;
  for (typename std::list<int>::iterator it=cycle.begin(); it!=cycle.end(); ++it) { // this is unnecesary
    //    if (DEBUG_MODE) {
    //      if (*it >= num_projection_nodes) {
    //        cout << "Cycle uses non-trivial partitionining." << std::endl;
    //      }
    //    }
    cycle_array[tmp_ind++] = *it;
  }
		
  std::vector<std::vector<std::vector<SIGNED_INDEX> > > pi_c(nNewClusters);
  // Found violated cycle, now triangulate and add to the relaxation!
  for(int i=0; i+1 < cycle.size()-2-i; i++) {
					
    // Add projection_imap_var applied to [i, i+1, cycle.size()-2-i]
    newCluster[cluster_index].i = projection_imap_var[cycle_array[i]];
    newCluster[cluster_index].j = projection_imap_var[cycle_array[i+1]];
    newCluster[cluster_index].k = projection_imap_var[cycle_array[cycle.size()-2-i]];

    pi_c[cluster_index].resize(3);
    pi_c[cluster_index][0] = pi[i];
    pi_c[cluster_index][1] = pi[i+1];
    pi_c[cluster_index][2] = pi[cycle.size()-2-i];
    
    cluster_index++;
  }

  for(int i=cycle.size()-1; i-1 > cycle.size()-1-i; i--) {
      
    // Add projection_imap_var applied to [i, i-1, cycle.size()-1-i]
    newCluster[cluster_index].i = projection_imap_var[cycle_array[i]];
    newCluster[cluster_index].j = projection_imap_var[cycle_array[i-1]];
    newCluster[cluster_index].k = projection_imap_var[cycle_array[cycle.size()-1-i]];
      
    pi_c[cluster_index].resize(3);
    pi_c[cluster_index][0] = pi[i];
    pi_c[cluster_index][1] = pi[i-1];
    pi_c[cluster_index][2] = pi[cycle.size()-1-i];
    
    cluster_index++;
  }

  // Add the top nclus_to_add clusters to the relaxation
  for(int clusterId = 0; clusterId < nNewClusters; clusterId++) {
    // Check that these clusters and intersection sets haven't already been added
    std::vector<int> temp;
    temp.push_back(newCluster[clusterId].i);
    temp.push_back(newCluster[clusterId].j);
    temp.push_back(newCluster[clusterId].k);
    sort(temp.begin(), temp.end());

    // Check to see if this cluster involves two of the same variables
    // (could happen because we didn't shortcut)
    if(temp[0] == temp[1] || temp[0] == temp[2] || temp[1] == temp[2]) {
      //      if(DEBUG_MODE)
      //	cout << "skipping this triplet because it is an edge." << std::endl;
      continue;
    }
    
    typename std::map<std::vector<int>, bool >::iterator t_itr = triplet_set.find(temp);
    if (t_itr == triplet_set.end())
      triplet_set.insert(std::pair<std::vector<int>, bool >(temp, true));
    else {
      //	if(DEBUG_MODE)
      //	  cout << "   Triplet was already present. Skipping..." << std::endl;
      continue; 
    }
				
    // add cluster ijk
    //std::vector<int> ijk_inds;
    //ijk_inds.push_back(newCluster[clusterId].i); ijk_inds.push_back(newCluster[clusterId].j); ijk_inds.push_back(newCluster[clusterId].k);
    std::array<int,3> ijk_inds = {{newCluster[clusterId].i,newCluster[clusterId].j,newCluster[clusterId].k}};
    
    //std::vector<int> ijk_intersect_inds;
    //ijk_intersect_inds.push_back(newCluster[clusterId].ij_intersect_loc);
    //ijk_intersect_inds.push_back(newCluster[clusterId].jk_intersect_loc);
    //ijk_intersect_inds.push_back(newCluster[clusterId].ki_intersect_loc);
	
    // do zrobienia: use sorting routine
    // sort indices
    if(ijk_inds[0] > ijk_inds[1]) {
      std::swap(ijk_inds[0],ijk_inds[1]);
      std::swap(pi_c[clusterId][0], pi_c[clusterId][1]);
    }
    if(ijk_inds[1] > ijk_inds[2]) {
      std::swap(ijk_inds[1],ijk_inds[2]);
      std::swap(pi_c[clusterId][1], pi_c[clusterId][2]);
    }
    if(ijk_inds[0] > ijk_inds[1]) {
      std::swap(ijk_inds[0],ijk_inds[1]);
      std::swap(pi_c[clusterId][0], pi_c[clusterId][1]);
    }
    //std::sort(ijk_inds.begin(), ijk_inds.end());
    assert(ijk_inds[0] < ijk_inds[1] && ijk_inds[1] < ijk_inds[2]);

    if(addTripletFun(ijk_inds[0], ijk_inds[1], ijk_inds[2])) {//, pi_c[clusterId][0], pi_c[clusterId][1], pi_c[clusterId][2])) {
       nClustersAdded++;
    }
    // TODO: log which clusters are chosen...
		
  }
  
  return nClustersAdded;
}

/**
 * Main function implementing UAI 2012 cycle tightening algorithm.
 * 
 * method=1: use create_k_projection_graph
 * method=2: use create_expanded_projection_graph
 */
template<typename MRF_CONSTRUCTOR>
template<typename ADD_TRIPLET_FUNCTION>
int 
Cycle<MRF_CONSTRUCTOR>::TightenCycle(
      ADD_TRIPLET_FUNCTION addTripletFun,
      const int nclus_to_add,  
      std::vector<int>& projection_imap_var,
      std::vector<std::vector<int> >& partition_imap,
      std::vector<std::list<int> >& cycle_set,
      const int method) 
{
   std::map<std::vector<int>, bool > triplet_set;

   int nClustersAdded = 0;
   //int nNewClusters;

   if (DEBUG_MODE) std::cout << "Finding the most violated cycle...." << std::endl;

   // This std::map allows us to quickly look up the edge weights
   std::map<std::pair<int, int>, REAL> projection_edge_weights;
   int num_projection_nodes;
   std::vector<std::vector<int> > projection_map;
   //std::vector<int> projection_imap_var;  
   std::vector<std::vector<std::pair<int, REAL> > > projection_adjacency_list;
   //std::vector<std::vector<int> > partition_imap;
   REAL* array_of_sij; int array_of_sij_size;

   // Define the projection graph and all edge weights
   if(method == 2)
      num_projection_nodes = create_expanded_projection_graph(projection_imap_var, projection_adjacency_list, projection_edge_weights, array_of_sij, array_of_sij_size, partition_imap);
   else if(method == 1)
      create_k_projection_graph(projection_map, num_projection_nodes, projection_imap_var, partition_imap, projection_edge_weights, projection_adjacency_list, array_of_sij, array_of_sij_size);
   else {
      std::cout << "ERROR: method not defined." << std::endl;
      return 0;
   }

   //std::vector<std::list<int> > cycle_set;
   REAL optimal_R = find_optimal_R(projection_adjacency_list, array_of_sij, array_of_sij_size);
   if (DEBUG_MODE) std::cout << "R_optimal = " << optimal_R << std::endl;

   //promised_bound = optimal_R;

   // Look for cycles. Some will be discarded.
   clock_t start_time = clock();

   if (optimal_R > 0) {
      // TODO: this is almost certainly doing more computation than necessary. Might want to change
      // nclus_to_add*10 to nclus_to_add, and comment out all but the top 3.
      FindCycles(cycle_set, optimal_R, nclus_to_add*10, projection_adjacency_list);
      FindCycles(cycle_set, optimal_R/2, nclus_to_add*10, projection_adjacency_list);
      FindCycles(cycle_set, optimal_R/4, nclus_to_add*10, projection_adjacency_list);
      FindCycles(cycle_set, optimal_R/8, nclus_to_add*10, projection_adjacency_list);
      FindCycles(cycle_set, optimal_R/16, nclus_to_add*10, projection_adjacency_list);
      FindCycles(cycle_set, optimal_R/32, nclus_to_add*10, projection_adjacency_list);
      FindCycles(cycle_set, optimal_R/64, nclus_to_add*10, projection_adjacency_list);
      FindCycles(cycle_set, optimal_R/128, nclus_to_add*10, projection_adjacency_list);
   }

   clock_t end_time = clock();
   REAL total_time = (REAL)(end_time - start_time)/CLOCKS_PER_SEC;
  if (DEBUG_MODE) {
     std::cout << " -- FindCycles. Took " << total_time << " seconds" << std::endl;
     std::cout << " Added " << cycle_set.size() << " cycles." << std::endl;
  } 

  
  // Add all cycles that we've found to the relaxation!
  //start_time = clock();
  for (int z = 0; z < cycle_set.size() && nClustersAdded < nclus_to_add; z++) {
     assert(cycle_set[z].size() >= 3);

     // Check to see if there are any duplicate nodes, and, if so, shortcut
     shortcut(cycle_set[z], projection_imap_var, projection_edge_weights, num_projection_nodes);

     // Output the cycle
     if (DEBUG_MODE){
        for (typename std::list<int>::iterator it=cycle_set[z].begin(); it!=cycle_set[z].end(); ++it) {
           std::cout << projection_imap_var[*it] << "(";
           std::vector<int> temp = partition_imap[*it];          
           for(int s=0; s < temp.size()-1; s++)
              std::cout << temp[s] << ",";
           std::cout << temp[temp.size()-1] << "), ";          
        }
        std::cout << std::endl;
     }

     // Add cycle to the relaxation by creating triplets
     nClustersAdded += add_cycle(addTripletFun, cycle_set[z], projection_imap_var, partition_imap, triplet_set, num_projection_nodes);
  }

  //end_time = clock();
  //total_time = (REAL)(end_time - start_time)/CLOCKS_PER_SEC;
  //if (DEBUG_MODE) {
  //   std::cout << " -- shortcut + add_cycles. Took " << total_time << " seconds" << std::endl;
  //} 
  
  delete_projection_graph(gm_.GetNumberOfVariables(), projection_map, projection_imap_var, projection_adjacency_list, array_of_sij);
  return nClustersAdded;
}

} // end namespace LP_MP

#endif // LP_MP_CYCLE_INEQUALITIES_HXX

