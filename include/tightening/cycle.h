/*
 *
 *  Created by Amir Globerson and David Sontag on 8/10/08.
 *  Updated by David Sontag, Do Kook Choe, and Yitao Li in 2012.
 *  Copyright 2008 MIT, 2012 NYU. All rights reserved.
 *
 */
#ifndef CYCLE_H
#define CYCLE_H

#include <stdlib.h>
#include <vector>
#include <algorithm>
#include <set>
#include <list>
#include <map>
#include <queue>
#include "mplp_alg.h"

#define DEBUG_MODE 1

#define Inf 9999999999.9
#define CYCLE_THRESH .00001  // TODO: Figure out how to set this
#define CLUSTER_THR .0000001

typedef map<pair<int, int>, int> mapType;
typedef vector<vector<pair<int, double> > > adj_type;

// Structure for storing a candidate triplet cluster for tightening
struct TripletCluster {
   double bound;
   int i,j,k;
   int ij_intersect_loc, jk_intersect_loc, ki_intersect_loc;

   bool operator <(const TripletCluster & rhs) const {
      return bound < rhs.bound;
   }
};


/////////////////////////////////////////////////////////////////////////////////
// Code for union-find data structure (used by FindPartition)

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

Node* find(Node* n) {
   if (n != n->parent) {
      n->parent = find(n->parent);
   }
   return n->parent;
} 


Node* merge(Node* x, Node* y) {
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

double maximizeIndependently(vector<MulDimArr*> & beliefs)
{
   double sum=0.0;
   int max_at; // not actually needed
   for(int i=0; i < beliefs.size(); i++)
   {
      sum += beliefs[i]->Max(max_at);
   }

   return sum;
}

double getValCycle(vector<MulDimArr*> & beliefs, vector<bool> & b_transpose, vector<int> & assignments)
{
   double sum=0.0;
   vector<int> inds; inds.push_back(-1); inds.push_back(-1); // temp

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

double maximizeCycle(vector<MulDimArr*> & beliefs, vector<bool> & b_transpose)
{
   double max_val = -Inf;

   // Fix value of the first variable
   int first_var_size = beliefs[0]->m_base_sizes[b_transpose[0]?1:0];
   int second_var_size = beliefs[0]->m_base_sizes[b_transpose[0]?0:1];
   for(int vo=0; vo < first_var_size; vo++)
   {
      vector<int> inds; inds.push_back(-1); inds.push_back(-1); // temp
      inds[b_transpose[0]?1:0] = vo;

      // Do first edge (construct initial field)
      vector<double> field;
      for(int v2=0; v2 < second_var_size; v2++)
      {
         inds[b_transpose[0]?0:1] = v2;
         field.push_back(beliefs[0]->GetVal(inds));
      }

      // Go over rest of edges, except last (which has to be treated specially)
      for(int i=1; i < beliefs.size()-1; i++)
      {
         vector<double> new_field;
         for(int v2=0; v2 < beliefs[i]->m_base_sizes[b_transpose[i]?0:1]; v2++)
         {
            inds.clear(); inds.push_back(-1); inds.push_back(-1); // temp
            inds[b_transpose[i]?0:1] = v2;

            // Take max
            double tmp_max_val = -Inf;
            for(int v1=0; v1 < field.size(); v1++)
            {
               inds[b_transpose[i]?1:0] = v1;
               tmp_max_val = max(tmp_max_val, field[v1]+beliefs[i]->GetVal(inds));
            }
            new_field.push_back(tmp_max_val);
         }
         field.clear(); // necessary?
         field = new_field;
      }

      // Do last edge (fix endpoint value to vo)
      inds.clear(); inds.push_back(-1); inds.push_back(-1); // temp
      inds[b_transpose[b_transpose.size()-1]?0:1] = vo;

      // Take max
      double tmp_max_val = -Inf;
      for(int v1=0; v1 < field.size(); v1++)
      {
         inds[b_transpose[b_transpose.size()-1]?1:0] = v1;
         tmp_max_val = max(tmp_max_val, field[v1]+beliefs[beliefs.size()-1]->GetVal(inds));
      }

      max_val = max(max_val, tmp_max_val);
   }
   return max_val;
}

/////////////////////////////////////////////////////////////////////////////////
// Implementation of UAI 2008 algorithm (just for triplets; square functionality removed)

int TightenTriplet(MPLPAlg & mplp, int nclus_to_add_min, int nclus_to_add_max, map<vector<int>, bool >& triplet_set, double & promised_bound) {

   int nClustersAdded = 0;
   int nNewClusters = 0;

   if(DEBUG_MODE)
      cout << "Doing pre-processing for adding triplet clusters." << endl;

   // Initialize adjacency list (filled in later) TODO: only do this when needed
   vector<int> adjacency_list[mplp.m_var_sizes.size()];

   // Construct adjacency list for the graph
   // Iterate over all of the edges (we do this by looking at the edge intersection sets)
   for(mapType::const_iterator it = mplp.m_intersect_map.begin(); it != mplp.m_intersect_map.end(); ++it)
   {
      // Get the two nodes i & j
      int i=it->first.first; int j=it->first.second;
      adjacency_list[i].push_back(j);
      adjacency_list[j].push_back(i);
   }

   // Sort the adjacency list, for fast intersections later
   for(int i=0; i < sizeof(adjacency_list)/sizeof(vector<int>); i++) 
   {
      sort(adjacency_list[i].begin(), adjacency_list[i].end());
   }

   // Count the number of triangles
   vector<int>::iterator intersects_iter_end;
   vector<int> commonNodes(mplp.m_var_sizes.size());
   for(mapType::const_iterator it = mplp.m_intersect_map.begin(); it != mplp.m_intersect_map.end(); ++it)
   {

      // Get the two nodes i & j
      int i=it->first.first; int j=it->first.second;

      // Now find all neighbors of both i and j to see where the triangles are
      intersects_iter_end = set_intersection(adjacency_list[i].begin(), adjacency_list[i].end(), adjacency_list[j].begin(), adjacency_list[j].end(), commonNodes.begin());

      for(vector<int>::const_iterator n=commonNodes.begin(); n != intersects_iter_end; ++n)
      {
         // Since a triplet shows up three times as an edge plus
         // a node, we only consider it for the case when n<i and n<j
         if(*n < i && *n < j)
            nNewClusters++;
      }
   }

   if(nNewClusters == 0) {
      if(DEBUG_MODE)
         cout << "nNewClusters = 0. Returning." << endl;

      return 0;
   }

   if(DEBUG_MODE)
      cout << "Looking for triangle clusters to add (" << nNewClusters << " triplets) " << endl;

   // TODO: put this elsewhere so that the space isn't re-allocated continuously?
   // Enumerate over all of the edges
   TripletCluster* newCluster = new TripletCluster[nNewClusters];

   int index=0;

   // Iterate over all of the edge intersection sets
   vector<int> tripAssignment; tripAssignment.push_back(-1); tripAssignment.push_back(-1); tripAssignment.push_back(-1);
   for(mapType::const_iterator it = mplp.m_intersect_map.begin(); it != mplp.m_intersect_map.end(); ++it)
   {

      // Get the two nodes i & j, and the edge intersection index
      int i=it->first.first; int j=it->first.second;
      int ij_intersect_loc = it->second;

      // Now find all neighbors of both i and j to see where the triangles are
      // TEMP TEMP -- fails at i=0, j=1, on i==3.
      intersects_iter_end = set_intersection(adjacency_list[i].begin(), adjacency_list[i].end(), adjacency_list[j].begin(), adjacency_list[j].end(), commonNodes.begin());

      for(vector<int>::const_iterator n=commonNodes.begin(); n != intersects_iter_end; ++n)
      {
         int k = *n;

         // Since a triplet shows up three times as an edge plus
         // a node, we only consider it for the case when k<i and k<j
         if(!(k < i && k < j))
            continue;

         newCluster[index].i = i;
         newCluster[index].j = j;
         newCluster[index].k = k;

         // Find the intersection sets for this triangle
         newCluster[index].ij_intersect_loc = ij_intersect_loc;
         vector<int> jk_edge; jk_edge.push_back(newCluster[index].j); jk_edge.push_back(newCluster[index].k);
         newCluster[index].jk_intersect_loc = mplp.FindIntersectionSet(jk_edge);
         vector<int> ki_edge; ki_edge.push_back(newCluster[index].k); ki_edge.push_back(newCluster[index].i);
         newCluster[index].ki_intersect_loc = mplp.FindIntersectionSet(ki_edge);					

         // Construct the beliefs for each edge, which will be maximized below
         vector<MulDimArr*> beliefs;
         vector<bool> b_transpose;
         beliefs.push_back(&mplp.m_sum_into_intersects[newCluster[index].ij_intersect_loc]);
         b_transpose.push_back(mplp.m_all_intersects[newCluster[index].ij_intersect_loc][0] != newCluster[index].i); // i first
         beliefs.push_back(&mplp.m_sum_into_intersects[newCluster[index].jk_intersect_loc]);
         b_transpose.push_back(mplp.m_all_intersects[newCluster[index].jk_intersect_loc][0] != newCluster[index].j); // then 'j'
         beliefs.push_back(&mplp.m_sum_into_intersects[newCluster[index].ki_intersect_loc]);
         b_transpose.push_back(mplp.m_all_intersects[newCluster[index].ki_intersect_loc][0] != newCluster[index].k); // then 'k'

         double bound_indep = maximizeIndependently(beliefs);

         // Before doing expensive joint maximization, see if we can quickly find optimal assignment
         tripAssignment[0] = mplp.m_decoded_res[i]; tripAssignment[1] = mplp.m_decoded_res[j]; tripAssignment[2] = mplp.m_decoded_res[k];
         double bound_quick = getValCycle(beliefs, b_transpose, tripAssignment);
         if(bound_indep == bound_quick)
         {
            newCluster[index].bound = 0;
         } else {
            // Do expensive joint maximization
            newCluster[index].bound = bound_indep - maximizeCycle(beliefs, b_transpose);
         }

         index++;
      }
   }

   // TODO opt: have a class for a cluster, so we can have different types and sort by bound,
   //       choosing best one by bound.
   //       Make the sorting and adding independent of the type of graph...

   // Sort the clusters by the bound
   sort(newCluster, newCluster+nNewClusters);

   if(DEBUG_MODE)
      printf(" -- Considered %d clusters, smallest bound %g, largest bound %g\n", nNewClusters, newCluster[max(nNewClusters-nclus_to_add_max, 0)].bound, newCluster[nNewClusters-1].bound);

   promised_bound = newCluster[nNewClusters-1].bound;

   // Add the top nclus_to_add clusters to the relaxation
   for(int clusterId = nNewClusters-1; clusterId >= 0 && nClustersAdded < nclus_to_add_max && (nClustersAdded < nclus_to_add_min || ((newCluster[clusterId].bound >= newCluster[nNewClusters-1].bound/5) && newCluster[clusterId].bound >= CLUSTER_THR)) ; clusterId--)
   {
      // Check to see if this triplet is already being used
      vector<int> temp;
      temp.push_back(newCluster[clusterId].i);
      temp.push_back(newCluster[clusterId].j);
      temp.push_back(newCluster[clusterId].k);
      sort(temp.begin(), temp.end());

      map<vector<int>, bool >::iterator t_itr = triplet_set.find(temp);
      if (t_itr == triplet_set.end())
         triplet_set.insert(pair<vector<int>, bool >(temp, true));
      else {
         //	if(DEBUG_MODE)
         //	  cout << "   Triplet was already present. Skipping..." << endl;
         continue; 
      }

      // Now add cluster ijk
      vector<int> ijk_inds;
      ijk_inds.push_back(newCluster[clusterId].i); ijk_inds.push_back(newCluster[clusterId].j); ijk_inds.push_back(newCluster[clusterId].k);

      vector<int> ijk_intersect_inds;
      ijk_intersect_inds.push_back(newCluster[clusterId].ij_intersect_loc);
      ijk_intersect_inds.push_back(newCluster[clusterId].jk_intersect_loc);
      ijk_intersect_inds.push_back(newCluster[clusterId].ki_intersect_loc);

      mplp.AddRegion(ijk_inds, ijk_intersect_inds);

      if(DEBUG_MODE)
         printf("Cluster added on nodes %d, %d, %d\n", newCluster[clusterId].i, newCluster[clusterId].j, newCluster[clusterId].k);
      // TODO: log which clusters are chosen...

      nClustersAdded++;
   }

   delete []newCluster;

   return nClustersAdded;
}


/////////////////////////////////////////////////////////////////////////////////
// Everything that follows is for the UAI 2012 cycle finding algorithm which
// can find arbitrary-length cycles to use in tightening the relaxation.

bool edge_sorting(list<int> i, list<int> j) {
   return i.back() > j.back();
}

bool edge_sorting3(pair<pair<int, int>, double> i, pair<pair<int, int>, double> j) {
   return i.second > j.second;
}

// De-allocate memory relating to the projection graph (TODO: finish writing this)
void delete_projection_graph(int num_vars, vector<vector<int> > &projection_map, vector<int> &projection_imap_var, adj_type &projection_adjacency_list, double* &array_of_sij) {

   //for(int node=0; node < num_vars; node++)
   //  delete []projection_map[node];
   //delete []projection_map;

   //delete projection_imap_var;

   //delete []projection_adjacency_list;

   delete []array_of_sij;
}

// Compute a single edge weight in the projection graph
double find_smn(const vector<int>& partition_i, int var_i_size, const vector<int>& partition_j, int var_j_size, MulDimArr* edge_belief) {

   int* whole_i = new int[var_i_size];
   int* whole_j = new int[var_j_size];  
   for (int i = 0; i < var_i_size; i++) {
      whole_i[i] = 0;
   }
   for (int i = 0; i < var_j_size; i++) {
      whole_j[i] = 0;
   }

   double smn = -Inf;
   for (int i = 0; i < partition_i.size(); i++) {
      whole_i[partition_i[i]] = 1;
      for (int j = 0; j < partition_j.size(); j++) {
         whole_j[partition_j[j]] = 1;
      }
   }

   double sec_max = -Inf;
   for (int i = 0; i < var_i_size; i++) {    
      for (int j = 0; j < var_j_size; j++) {
         if (whole_i[i] == whole_j[j]) {
            vector<int> inds; inds.push_back(i); inds.push_back(j);
            double temp_val = edge_belief->GetVal(inds);
            if (smn < temp_val) smn = temp_val;
         }
         else {
            vector<int> inds; inds.push_back(i); inds.push_back(j);
            double temp_val = edge_belief->GetVal(inds);
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
double find_smn_state_i(int single_i, int var_i_size, const vector<int>& partition_j, int var_j_size, MulDimArr* edge_belief, vector<vector<double> >& max_i_bij_not_xi) {
   double max, sec_max;
   max = sec_max = -Inf;
   double whole_i[var_i_size];
   double whole_j[var_j_size];
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
      vector<int> inds; inds.push_back(single_i); inds.push_back(partition_j[j]);
      double tmp = edge_belief->GetVal(inds);
      if (max < tmp) max = tmp;
   }

   for (int j = 0; j < var_j_size; j++) {
      if (whole_j[j] == 0) {
         //bruteforce max_{pi(x_i)=pi(x_j)=0}bij(x_i,x_j)
         vector<int> inds; inds.push_back(single_i); inds.push_back(j);
         double tmp = edge_belief->GetVal(inds);
         if (sec_max < tmp) sec_max = tmp; 

         if (max < max_i_bij_not_xi[single_i][j]) max = max_i_bij_not_xi[single_i][j];
      }
      else {
         if (sec_max < max_i_bij_not_xi[single_i][j]) sec_max = max_i_bij_not_xi[single_i][j];			
      }
   }	

   return max - sec_max;
}

// Compute a single edge weight in the projection graph (more efficiently)
double find_smn_state_j(const vector<int>& partition_i, int var_i_size, int single_j, int var_j_size, MulDimArr* edge_belief, vector<vector<double> >& max_j_bij_not_xj) {
   double max, sec_max;
   max = sec_max = -Inf;
   double whole_i[var_i_size];
   double whole_j[var_j_size];
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
      vector<int> inds; inds.push_back(partition_i[i]); inds.push_back(single_j);
      double tmp = edge_belief->GetVal(inds);
      if (max < tmp) max = tmp;
   }

   for (int i = 0; i < var_i_size; i++) {
      if (whole_i[i] == 0) {
         //bruteforce max_{pi(x_i)=pi(x_j)=0}bij(x_i,x_j)
         vector<int> inds; inds.push_back(i); inds.push_back(single_j);
         double tmp = edge_belief->GetVal(inds);
         if (sec_max < tmp) sec_max = tmp; 

         if (max < max_j_bij_not_xj[i][single_j]) max = max_j_bij_not_xj[i][single_j];
      }
      else {
         if (sec_max < max_j_bij_not_xj[i][single_j]) sec_max = max_j_bij_not_xj[i][single_j];			
      }
   }	

   return max - sec_max;
}


// This will find the best partitioning of the states of a single edge
// Implements the FindPartition algorithm described in supplementary material of UAI 2012 paper.
void find_partition(vector<map<vector<int>, int> >& partition_set, MPLPAlg& mplp, MulDimArr* edge_belief, int index_i, int index_j) {

   int size_i = mplp.m_var_sizes[index_i];
   int size_j = mplp.m_var_sizes[index_j];

   vector<Node*> x_i; //partition for node i
   vector<Node*> x_j; //partition for node j

   for (int i = 0; i < size_i; i++) {
      Node* t = new Node(2);
      x_i.push_back(t);
   }
   for (int j = 0; j < size_j; j++) {
      Node* t = new Node(1);
      x_j.push_back(t);
   }

   // sort edges by their weights in decreasing order
   vector<pair<pair<int, int>, double> > sorted_edge_set;
   for (int i = 0; i < size_i; i++) {
      for (int j = 0; j < size_j; j++) {
         vector<int> inds; inds.push_back(i); inds.push_back(j);
         double temp_val = edge_belief->GetVal(inds);
         sorted_edge_set.push_back(pair<pair<int, int>, double>(make_pair(make_pair(i, j), temp_val)));
      }
   }
   sort(sorted_edge_set.begin(), sorted_edge_set.end(), edge_sorting3);

   vector<pair<pair<int, int>, double> >::iterator it = sorted_edge_set.begin();
   double val_s = it->second; //val_s = bij(x_i^*, x_j^*);

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
      cout << "something is wrong" << endl;
   }

   if (size_i == 1 || size_j == 1) {
      cout << "should not happen" << endl;
      return;
   }
   vector<int> part_i, part_j; 

   //find the smallest partition and use it as an index
   // If there's a tie, use the partition with the smallest index. Keep it sorted.
   int min = mplp.m_var_sizes[index_i] + 1; 
   Node* head;
   for (int i = 0; i < mplp.m_var_sizes[index_i]; i++) {
      Node* t = find(x_i[i]);
      if ((size_i == 2 || t->bit == 3) && min > t->i_size) {
         min = t->i_size;
         head = t;
      }
   }
   for (int i = 0; i < mplp.m_var_sizes[index_i]; i++) {
      Node* t = find(x_i[i]);
      if (head == t)
         part_i.push_back(i);
   }

   min = mplp.m_var_sizes[index_j] + 1;
   for (int j = 0; j < mplp.m_var_sizes[index_j]; j++) {
      Node* t = find(x_j[j]);
      if ((size_j == 2 || t->bit == 3) && min > t->j_size) {
         min = t->j_size;
         head = t;
      }
   }
   for (int j = 0; j < mplp.m_var_sizes[index_j]; j++) {
      Node* t = find(x_j[j]);
      if (head == t)
         part_j.push_back(j);
   } 

   // We have one map per variable.

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

   for (int i = 0; i < mplp.m_var_sizes[index_i]; i++) {
      delete [] x_i[i];
   }
   for (int j = 0; j < mplp.m_var_sizes[index_j]; j++) {
      delete [] x_j[j];
   }
}


// Create the expanded projection graph by including all singleton partitions and also
// all partitions found by calling FindPartition on all edges.
int create_expanded_projection_graph(MPLPAlg& mplp, vector<int>& projection_imap_var, vector<vector<pair<int, double> > >& projection_adjacency_list, map<pair<int, int>, double>&projection_edge_weights, double* &array_of_sij, int& array_of_sij_size, vector<vector<int> >& partition_imap) {
   // projection_imap_var maps from projection node to variable
   // partition_imap maps from projection node to vector of states

   int num_of_vars = mplp.m_var_sizes.size();
   vector<map<vector<int>,int> > partition_set;
   set<double> set_of_sij;
   int num_projection_nodes = 0;

   for (int i = 0; i < num_of_vars; i++) {
      map<vector<int>,int> p;
      partition_set.push_back(p);

      // push all singleton partitions
      for (int i_state = 0; i_state < mplp.m_var_sizes[i]; i_state++) {
         vector<int> single;
         single.push_back(i_state);
         partition_set[i][single] = 1;
      }
   }

   bool partition = true;
   clock_t find_partition_start_time;
   if (partition) {
      find_partition_start_time = clock();
      //THIS IS THE NEW PARTITIONING ALGORITHM

      for (mapType::const_iterator it = mplp.m_intersect_map.begin(); it != mplp.m_intersect_map.end(); ++it) {    
         int i=it->first.first; int j=it->first.second;
         int ij_intersect_loc = it->second;

         MulDimArr* edge_belief = &mplp.m_sum_into_intersects[ij_intersect_loc];
         if(mplp.m_all_intersects[ij_intersect_loc][0] != i)  { // swap i and j
            int tmp_i = i;
            i = j;
            j = tmp_i;
         }

         // Check to see if i and j have at least two states each -- otherwise, cannot be part of any frustrated edge
         if(mplp.m_var_sizes[i] <= 1 || mplp.m_var_sizes[j] <= 1)
            continue;

         find_partition(partition_set, mplp, edge_belief, i, j);
      }
      clock_t find_partition_end_time = clock();
      double find_partition_total_time = (double)(find_partition_end_time - find_partition_start_time)/CLOCKS_PER_SEC;
      if (DEBUG_MODE) {
         printf(" -- find_partition. Took %lg seconds\n", find_partition_total_time);
      }
   }

   // Create one projection node for each remaining partition. Overload value of partition_set
   // to now refer to the node number in the projection graph.
   for(int i=0; i < partition_set.size(); i++) {
      for(map<vector<int>, int>::iterator it = partition_set[i].begin(); it != partition_set[i].end(); it++) {
         it->second = num_projection_nodes++;
         projection_imap_var.push_back(i);
         partition_imap.push_back(it->first);

         vector<pair<int, double> > temp;
         projection_adjacency_list.push_back(temp);
      }
   }

   clock_t find_smn_start_time = clock();

   // Create projection graph edges for each edge of original graph and each pair of partitions
   for (mapType::const_iterator it = mplp.m_intersect_map.begin(); it != mplp.m_intersect_map.end(); ++it) {    
      int i = it->first.first; int j = it->first.second;
      int ij_intersect_loc = it->second;

      // Do some pre-processing to make the case of single state partitions very fast

      MulDimArr* edge_belief = &mplp.m_sum_into_intersects[ij_intersect_loc];
      if (mplp.m_all_intersects[ij_intersect_loc][0] != i) {
         int tmp_i = i;
         i = j;
         j = tmp_i;
      }

      if (mplp.m_var_sizes[i] <= 1 || mplp.m_var_sizes[j] <= 1) continue;

      vector<vector<double> > max_i_bij_not_xi;
      vector<vector<double> > max_j_bij_not_xj;

      for(int state1=0; state1 < mplp.m_var_sizes[i]; state1++) {

         // Find max over state2
         double largest_val = -Inf;
         int largest_ind = -1;
         for(int state2=0; state2 < mplp.m_var_sizes[j]; state2++) {

            vector<int> inds; inds.push_back(state1); inds.push_back(state2);
            double tmp_val = edge_belief->GetVal(inds);

            if(tmp_val > largest_val) {
               largest_val = tmp_val;
               largest_ind = state2;
            }
         }

         // Find second largest val over state2
         double sec_largest_val = -Inf;
         for(int state2=0; state2 < mplp.m_var_sizes[j]; state2++) {

            vector<int> inds; inds.push_back(state1); inds.push_back(state2);
            double tmp_val = edge_belief->GetVal(inds);

            if(tmp_val > sec_largest_val && state2 != largest_ind) {
               sec_largest_val = tmp_val;
            }
         }

         // Assign values
         for(int state2=0; state2 < mplp.m_var_sizes[j]; state2++) {
            vector<double> state_j; max_j_bij_not_xj.push_back(state_j);
            max_j_bij_not_xj[state1].push_back(largest_val);
         }
         max_j_bij_not_xj[state1][largest_ind] = sec_largest_val;
      }


      for(int state1=0; state1 < mplp.m_var_sizes[j]; state1++) {

         // Find max over state2
         double largest_val = -Inf;
         int largest_ind = -1;
         for(int state2=0; state2 < mplp.m_var_sizes[i]; state2++) {

            vector<int> inds; inds.push_back(state2); inds.push_back(state1);
            double tmp_val = edge_belief->GetVal(inds);

            if(tmp_val > largest_val) {
               largest_val = tmp_val;
               largest_ind = state2;
            }
         }

         // Find second largest val over state2
         double sec_largest_val = -Inf;
         for(int state2=0; state2 < mplp.m_var_sizes[i]; state2++) {

            vector<int> inds; inds.push_back(state2); inds.push_back(state1);
            double tmp_val = edge_belief->GetVal(inds);

            if(tmp_val > sec_largest_val && state2 != largest_ind) {
               sec_largest_val = tmp_val;
            }
         }

         // Assign values
         for(int state2=0; state2 < mplp.m_var_sizes[i]; state2++) {
            vector<double> state_i; max_i_bij_not_xi.push_back(state_i);
            max_i_bij_not_xi[state2].push_back(largest_val);
         }
         max_i_bij_not_xi[largest_ind][state1] = sec_largest_val;
      }


      // decompose: max_{x_j!=x_j}max_{x_i != x_i'}. Then, use the above computations.
      double max_ij_bij_not_xi_xj[mplp.m_var_sizes[i]][mplp.m_var_sizes[j]];
      for(int state1=0; state1 < mplp.m_var_sizes[i]; state1++) {

         // Find max over state2
         double largest_val = -Inf;
         int largest_ind = -1;
         for(int state2=0; state2 < mplp.m_var_sizes[j]; state2++) {
            double tmp_val = max_i_bij_not_xi[state1][state2];

            if(tmp_val > largest_val) {
               largest_val = tmp_val;
               largest_ind = state2;
            }
         }

         // Find second largest val over state2
         double sec_largest_val = -Inf;
         for(int state2=0; state2 < mplp.m_var_sizes[j]; state2++) {
            double tmp_val = max_i_bij_not_xi[state1][state2];

            if(tmp_val > sec_largest_val && state2 != largest_ind) {
               sec_largest_val = tmp_val;
            }
         }

         // Assign values
         for(int state2=0; state2 < mplp.m_var_sizes[j]; state2++)
            max_ij_bij_not_xi_xj[state1][state2] = largest_val;
         max_ij_bij_not_xi_xj[state1][largest_ind] = sec_largest_val;
      }

      // Now, for each partition of node i and each partition of node j, compute edge weights
      // If the edge weight is non-zero, insert edge into adjacency list

      for(map<vector<int>, int>::iterator it_i = partition_set[i].begin(); it_i != partition_set[i].end(); it_i++) {
         int n = it_i->second;

         for(map<vector<int>, int>::iterator it_j = partition_set[j].begin(); it_j != partition_set[j].end(); it_j++) {
            int m = it_j->second;

            double smn = 0;
            if (it_i->first.size() == 1 && it_j->first.size() == 1) {
               int xi = it_i->first[0]; int xj = it_j->first[0];
               vector<int> inds; inds.push_back(xi); inds.push_back(xj);
               double tmp_val = edge_belief->GetVal(inds);
               smn = max(tmp_val, max_ij_bij_not_xi_xj[xi][xj]) - max(max_i_bij_not_xi[xi][xj], max_j_bij_not_xj[xi][xj]);
            }
            else if (it_i->first.size() == 1) {
               int xi = it_i->first[0];
               smn = find_smn_state_i(xi, mplp.m_var_sizes[i], it_j->first, mplp.m_var_sizes[j], edge_belief, max_i_bij_not_xi);
            }
            else if (it_j->first.size() == 1) {
               int xj = it_j->first[0];
               smn = find_smn_state_j(it_i->first, mplp.m_var_sizes[i], xj, mplp.m_var_sizes[j], edge_belief, max_j_bij_not_xj);
            }
            else {
               // This computes smn
               smn = find_smn(it_i->first, mplp.m_var_sizes[i], it_j->first, mplp.m_var_sizes[j], edge_belief);
            }

            if (smn != 0) {
               set_of_sij.insert(fabs(smn));

               projection_adjacency_list[m].push_back(make_pair(n, smn));
               projection_adjacency_list[n].push_back(make_pair(m, smn));

               projection_edge_weights.insert(pair<pair<int, int>, double>(pair<int, int>(n, m), smn));
               projection_edge_weights.insert(pair<pair<int, int>, double>(pair<int, int>(m, n), smn));
            }

         }
      }
   }
   clock_t find_smn_end_time = clock();
   double find_smn_total_time = (double)(find_smn_end_time - find_smn_start_time)/CLOCKS_PER_SEC;
   if (DEBUG_MODE) {
      printf(" -- find_smn. Took %lg seconds\n", find_smn_total_time);
   }

   array_of_sij = new double[set_of_sij.size()];
   array_of_sij_size = 0;

   for (set<double>::iterator set_iter = set_of_sij.begin(); set_iter != set_of_sij.end(); set_iter++) {  //NOTE: PASSED AS FUNCTION ARG, DO NOT DELETE HERE!!
      array_of_sij[array_of_sij_size++] = *set_iter;
   }

   return num_projection_nodes;
}


// Creates the k-projection graph (just a single partition per variable)
void create_k_projection_graph(MPLPAlg& mplp, vector<vector<int> > &projection_map, int& num_projection_nodes, vector<int> &projection_imap_var, vector<vector<int> > &partition_imap, map<pair<int, int>, double>& projection_edge_weights, adj_type &projection_adjacency_list, double* &array_of_sij, int& array_of_sij_size) {

   // TODO: make sure for binary variables that there is only one node per variable (rather than 2).
   // TODO: projection_edge_weights can likely be removed from this function and elsewhere.
   // TODO: most of these ints can be changed to be unsigned and/or fewer bits. Look into memory allocation.

   // Initialize the projection graph
   int projection_node_iter = 0;
   for(int node=0; node < mplp.m_var_sizes.size(); node++) {
      vector<int> i;

      for(int state=0; state < mplp.m_var_sizes[node]; state++) {
         vector<pair<int, double> > temp;
         projection_adjacency_list.push_back(temp);
         i.push_back(projection_node_iter++);
      }
      projection_map.push_back(i);
   }
   // Construct inverse map
   num_projection_nodes = projection_node_iter;

   projection_node_iter = 0;
   for(int node=0; node < mplp.m_var_sizes.size(); node++) {
      for(int state=0; state < mplp.m_var_sizes[node]; state++) {
         projection_imap_var.push_back(node);

         vector<int> tmp_state_vector;
         tmp_state_vector.push_back(state);
         partition_imap.push_back(tmp_state_vector);

         projection_node_iter++;
      }
   }

   // Iterate over all of the edges (we do this by looking at the edge intersection sets)
   set<double> set_of_sij;
   for(mapType::const_iterator it = mplp.m_intersect_map.begin(); it != mplp.m_intersect_map.end(); ++it) {
      // Get the two nodes i & j and the edge intersection set. Put in right order.
      int i=it->first.first; int j=it->first.second;
      int ij_intersect_loc = it->second;

      MulDimArr* edge_belief = &mplp.m_sum_into_intersects[ij_intersect_loc];
      if(mplp.m_all_intersects[ij_intersect_loc][0] != i)  { // swap i and j
         int tmp_i = i;
         i = j;
         j = tmp_i;
      }

      // Check to see if i and j have at least two states each -- otherwise, cannot be part of any frustrated edge
      if(mplp.m_var_sizes[i] <= 1 || mplp.m_var_sizes[j] <= 1)
         continue;

      // Do some pre-processing for speed.
      double max_j_bij_not_xj[mplp.m_var_sizes[i]][mplp.m_var_sizes[j]];
      double max_i_bij_not_xi[mplp.m_var_sizes[i]][mplp.m_var_sizes[j]];

      for(int state1=0; state1 < mplp.m_var_sizes[i]; state1++) {

         // Find max over state2
         double largest_val = -Inf;
         int largest_ind = -1;
         for(int state2=0; state2 < mplp.m_var_sizes[j]; state2++) {

            vector<int> inds; inds.push_back(state1); inds.push_back(state2);
            double tmp_val = edge_belief->GetVal(inds);

            if(tmp_val > largest_val) {
               largest_val = tmp_val;
               largest_ind = state2;
            }
         }

         // Find second largest val over state2
         double sec_largest_val = -Inf;
         for(int state2=0; state2 < mplp.m_var_sizes[j]; state2++) {

            vector<int> inds; inds.push_back(state1); inds.push_back(state2);
            double tmp_val = edge_belief->GetVal(inds);

            if(tmp_val > sec_largest_val && state2 != largest_ind) {
               sec_largest_val = tmp_val;
            }
         }

         // Assign values
         for(int state2=0; state2 < mplp.m_var_sizes[j]; state2++)
            max_j_bij_not_xj[state1][state2] = largest_val;
         max_j_bij_not_xj[state1][largest_ind] = sec_largest_val;
      }


      for(int state1=0; state1 < mplp.m_var_sizes[j]; state1++) {

         // Find max over state2
         double largest_val = -Inf;
         int largest_ind = -1;
         for(int state2=0; state2 < mplp.m_var_sizes[i]; state2++) {

            vector<int> inds; inds.push_back(state2); inds.push_back(state1);
            double tmp_val = edge_belief->GetVal(inds);

            if(tmp_val > largest_val) {
               largest_val = tmp_val;
               largest_ind = state2;
            }
         }

         // Find second largest val over state2
         double sec_largest_val = -Inf;
         for(int state2=0; state2 < mplp.m_var_sizes[i]; state2++) {

            vector<int> inds; inds.push_back(state2); inds.push_back(state1);
            double tmp_val = edge_belief->GetVal(inds);

            if(tmp_val > sec_largest_val && state2 != largest_ind) {
               sec_largest_val = tmp_val;
            }
         }

         // Assign values
         for(int state2=0; state2 < mplp.m_var_sizes[i]; state2++)
            max_i_bij_not_xi[state2][state1] = largest_val;
         max_i_bij_not_xi[largest_ind][state1] = sec_largest_val;
      }

      // decompose: max_{x_j!=x_j}max_{x_i != x_i'}. Then, use the above computations.
      double max_ij_bij_not_xi_xj[mplp.m_var_sizes[i]][mplp.m_var_sizes[j]];
      for(int state1=0; state1 < mplp.m_var_sizes[i]; state1++) {

         // Find max over state2
         double largest_val = -Inf;
         int largest_ind = -1;
         for(int state2=0; state2 < mplp.m_var_sizes[j]; state2++) {
            double tmp_val = max_i_bij_not_xi[state1][state2];

            if(tmp_val > largest_val) {
               largest_val = tmp_val;
               largest_ind = state2;
            }
         }

         // Find second largest val over state2
         double sec_largest_val = -Inf;
         for(int state2=0; state2 < mplp.m_var_sizes[j]; state2++) {
            double tmp_val = max_i_bij_not_xi[state1][state2];

            if(tmp_val > sec_largest_val && state2 != largest_ind) {
               sec_largest_val = tmp_val;
            }
         }

         // Assign values
         for(int state2=0; state2 < mplp.m_var_sizes[j]; state2++)
            max_ij_bij_not_xi_xj[state1][state2] = largest_val;
         max_ij_bij_not_xi_xj[state1][largest_ind] = sec_largest_val;
      }

      // For each of their states
      for(int xi=0; xi < mplp.m_var_sizes[i]; xi++) {
         int m = projection_map[i][xi];

         for(int xj=0; xj < mplp.m_var_sizes[j]; xj++) {
            int n = projection_map[j][xj];

            vector<int> inds; inds.push_back(xi); inds.push_back(xj);
            double tmp_val = edge_belief->GetVal(inds);

            // Compute s_mn for this edge
            double val_s = max(tmp_val, max_ij_bij_not_xi_xj[xi][xj]) - max(max_i_bij_not_xi[xi][xj], max_j_bij_not_xj[xi][xj]);

            // TODO: use threshold here, to make next stage faster
            if(val_s != 0) {          
               projection_adjacency_list[m].push_back(make_pair(n, val_s));
               projection_adjacency_list[n].push_back(make_pair(m, val_s));
               set_of_sij.insert(fabs(val_s));
            }

            // Insert into edge weight map
            projection_edge_weights.insert(pair<pair<int,int>,double>(pair<int,int>(n, m), val_s));
            projection_edge_weights.insert(pair<pair<int,int>,double>(pair<int,int>(m, n), val_s));
         }
      }

   }

   // Sort list_of_sij and remove duplicates
   array_of_sij = new double[set_of_sij.size()];
   array_of_sij_size = 0;

   for(set<double>::iterator set_iter = set_of_sij.begin(); set_iter != set_of_sij.end(); ++set_iter) {
      array_of_sij[array_of_sij_size++] = *set_iter;
   }
}


// Does binary search over the edge weights to find largest edge weight
// such that there is an odd-signed cycle.
double find_optimal_R(adj_type &projection_adjacency_list, double* &array_of_sij, int array_of_sij_size) {

   // Do binary search over sij
   int bin_search_lower_bound = 0;
   int bin_search_upper_bound = array_of_sij_size;
   double sij_min = -1;
   int num_projection_nodes = projection_adjacency_list.size();

   while(bin_search_lower_bound <= bin_search_upper_bound) {

      // Compute mid-point
      int R_pos = floor((bin_search_lower_bound + bin_search_upper_bound)/2);
      double R = array_of_sij[R_pos];
      std::cout << "Test: " << bin_search_lower_bound << " < " << R_pos << " < " << bin_search_upper_bound << ", R = " << R << std::endl;


      // Does there exist an odd signed cycle using just edges with sij >= R? If so, go up. If not, go down.
      bool found_odd_signed_cycle = false;

      // Initialize
      // TODO: do not allocate memory for this every time! Just re-initialize.

      int node_sign[num_projection_nodes];
      for(int m=0; m<num_projection_nodes; m++) {
         node_sign[m] = 0; // denotes "not yet seen"      
      }

      // Graph may be disconnected, so check from all nodes    
      for(int i = 0; i < num_projection_nodes && !found_odd_signed_cycle; i++) {
         if(node_sign[i] == 0) {
            node_sign[i] = 1; // root node
            queue<int> q;
            q.push(i);

            while (!q.empty() && !found_odd_signed_cycle) {
               int current = q.front();
               q.pop();

               for (int j = 0; j < projection_adjacency_list[current].size(); j++) {
                  double smn = projection_adjacency_list[current][j].second;
                  // Ignore edges with weight less than R
                  if (fabs(smn) < R)
                     continue;

                  int next = projection_adjacency_list[current][j].first;
                  int sign_of_smn = (smn > 0) - (smn < 0);

                  if (node_sign[next] == 0) {
                     node_sign[next] = node_sign[current] * sign_of_smn;
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

   return sij_min;
}


// Returns an array which is a random permutation of the numbers 0 through n-1.
int* random_permutation(int n) {
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
void FindCycles(vector<list<int> > &cycle_set, double optimal_R, int ncycles_to_add, adj_type &projection_adjacency_list) {

   double R = optimal_R;
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
   //    cout << "First random node is: " << random_projection_node[0] << endl;
   //  }

   for(int ri = 0; ri < num_projection_nodes; ri++) {
      int i = random_projection_node[ri];    
      if(node_sign[i] == 0) {
         node_sign[i] = 1; // root node
         node_parent[i] = i;
         node_depth[i] = 0;

         queue<int> q;
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
                  double smn = projection_adjacency_list[current][j].second;
                  if (fabs(smn) < R) {
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
   vector<list<int> > edge_set;
   map<pair<int, int>, int> edge_map;

   for (int i = 0; i < num_projection_nodes; i++) {
      for (int j = 0; j < projection_adjacency_list[i].size(); j++) {
         if (node_parent[i] == j || node_parent[j] == i) {							
            continue;
         }
         double smn = projection_adjacency_list[i][j].second;

         if (fabs(smn) < R) {
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
                  list<int> temp;
                  temp.push_back(i);
                  temp.push_back(jj);
                  temp.push_back(depth_of_i - temp_di + depth_of_j - temp_dj);
                  map<pair<int, int>, int>::iterator m_it = edge_map.find(make_pair(i, jj));
                  if (m_it == edge_map.end()) {
                     edge_map[pair<int, int>(i, jj)] = depth_of_i - temp_di + depth_of_j - temp_dj;
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
   sort(edge_set.begin(), edge_set.end(), edge_sorting);	

   // Note: here we do not check for duplicate variables 
   for (int i = 0; i < edge_set.size() && cycle_set.size() < ncycles_to_add; i++) {

      // Get next edge (and depth) from sorted list
      list<int> temp = edge_set.back();
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

      list<int> cycle;

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
// NOTE: this is not currently being used.
void shortcut(list<int> &cycle, vector<int> &projection_imap_var, map<pair<int, int>, double>& projection_edge_weights, int& num_projection_nodes) {

   bool exist_duplicates = true;
   int cycle_array[cycle.size()]; int tmp_ind = 0;
   int cycle_sign[cycle.size()];

   for (list<int>::iterator it=cycle.begin(); it!=cycle.end(); ++it) {
      cycle_array[tmp_ind++] = *it;
   }

   int cycle_start = 0;
   int cycle_end = cycle.size()-1;

   // first, shorten the cycle
   while(exist_duplicates) {
      exist_duplicates = false;

      // This map allows us to quickly find duplicates and their locations within the cycle
      map<int, int> duplicates_map;

      // Must keep track of sign.
      cycle_sign[cycle_start] = 1;

      duplicates_map.insert(pair<int,int>(projection_imap_var[cycle_array[cycle_start]], cycle_start));
      for(int i=cycle_start+1; i <= cycle_end; i++) {

         // Get edge weight
         map< pair<int,int>, double>::iterator weights_iter = projection_edge_weights.find(make_pair(cycle_array[i-1], cycle_array[i]));

         double smn = weights_iter->second;
         int sign_of_smn = (smn > 0) - (smn < 0);
         cycle_sign[i] = cycle_sign[i-1]*sign_of_smn;

         // Does node i already exist in the map?
         map<int, int>::iterator dups_iter = duplicates_map.find(projection_imap_var[cycle_array[i]]);
         if( dups_iter != duplicates_map.end() ) { // Duplicate found
            int first_occurence_index = dups_iter->second;

            // Look to see if the first half of the cycle is violated.
            weights_iter = projection_edge_weights.find(make_pair(cycle_array[first_occurence_index], cycle_array[i-1]));        
            double edge_weight = weights_iter->second;
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
            // Insert into the map
            duplicates_map.insert(pair<int,int>(projection_imap_var[cycle_array[i]], i));
         }
      }
   }

   // Modify the cycle
   cycle.clear();
   for(int i=cycle_start; i <= cycle_end; i++)
      cycle.push_back(cycle_array[i]);
}


int add_cycle(MPLPAlg& mplp, list<int> &cycle, vector<int> &projection_imap_var, map<vector<int>, bool >& triplet_set, int& num_projection_nodes) {

   int nClustersAdded = 0;

   // Number of clusters we're adding is length_cycle - 2
   int nNewClusters = cycle.size() - 2;
   TripletCluster newCluster[nNewClusters];
   int cluster_index = 0;

   // Convert cycle to array
   int cycle_array[cycle.size()]; int tmp_ind = 0;
   for (list<int>::iterator it=cycle.begin(); it!=cycle.end(); ++it) { // this is unnecesary
      //    if (DEBUG_MODE) {
      //      if (*it >= num_projection_nodes) {
      //        cout << "Cycle uses non-trivial partitionining." << endl;
      //      }
      //    }
      cycle_array[tmp_ind++] = *it;
   }

   // Found violated cycle, now triangulate and add to the relaxation!
   for(int i=0; i+1 < cycle.size()-2-i; i++) {

      // Add projection_imap_var applied to [i, i+1, cycle.size()-2-i]
      newCluster[cluster_index].i = projection_imap_var[cycle_array[i]];
      newCluster[cluster_index].j = projection_imap_var[cycle_array[i+1]];
      newCluster[cluster_index].k = projection_imap_var[cycle_array[cycle.size()-2-i]];

      cluster_index++;
   }

   for(int i=cycle.size()-1; i-1 > cycle.size()-1-i; i--) {

      // Add projection_imap_var applied to [i, i-1, cycle.size()-1-i]
      newCluster[cluster_index].i = projection_imap_var[cycle_array[i]];
      newCluster[cluster_index].j = projection_imap_var[cycle_array[i-1]];
      newCluster[cluster_index].k = projection_imap_var[cycle_array[cycle.size()-1-i]];

      cluster_index++;
   }

   // Add the top nclus_to_add clusters to the relaxation
   for(int clusterId = 0; clusterId < nNewClusters; clusterId++) {
      // Check that these clusters and intersection sets haven't already been added
      vector<int> temp;
      temp.push_back(newCluster[clusterId].i);
      temp.push_back(newCluster[clusterId].j);
      temp.push_back(newCluster[clusterId].k);
      sort(temp.begin(), temp.end());

      // Check to see if this cluster involves two of the same variables
      // (could happen because we didn't shortcut)
      if(temp[0] == temp[1] || temp[0] == temp[2] || temp[1] == temp[2]) {
         //      if(DEBUG_MODE)
         //	cout << "skipping this triplet because it is an edge." << endl;
         continue;
      }

      map<vector<int>, bool >::iterator t_itr = triplet_set.find(temp);
      if (t_itr == triplet_set.end())
         triplet_set.insert(pair<vector<int>, bool >(temp, true));
      else {
         //	if(DEBUG_MODE)
         //	  cout << "   Triplet was already present. Skipping..." << endl;
         continue; 
      }

      // Find the intersection sets for this triangle
      vector<int> ij_edge; ij_edge.push_back(newCluster[clusterId].i); ij_edge.push_back(newCluster[clusterId].j);
      newCluster[clusterId].ij_intersect_loc = mplp.FindIntersectionSet(ij_edge);
      // This edge intersection set may not already exist, in which case we should add it
      if(newCluster[clusterId].ij_intersect_loc == -1) {
         newCluster[clusterId].ij_intersect_loc = mplp.AddIntersectionSet(ij_edge);
      }

      vector<int> jk_edge; jk_edge.push_back(newCluster[clusterId].j); jk_edge.push_back(newCluster[clusterId].k);
      newCluster[clusterId].jk_intersect_loc = mplp.FindIntersectionSet(jk_edge);
      // This edge intersection set may not already exist, in which case we should add it
      if(newCluster[clusterId].jk_intersect_loc == -1) {
         newCluster[clusterId].jk_intersect_loc = mplp.AddIntersectionSet(jk_edge);
      }

      vector<int> ki_edge; ki_edge.push_back(newCluster[clusterId].k); ki_edge.push_back(newCluster[clusterId].i);
      newCluster[clusterId].ki_intersect_loc = mplp.FindIntersectionSet(ki_edge);
      // This edge intersection set may not already exist, in which case we should add it
      if(newCluster[clusterId].ki_intersect_loc == -1) {
         newCluster[clusterId].ki_intersect_loc = mplp.AddIntersectionSet(ki_edge);
      }			

      // Now add cluster ijk
      vector<int> ijk_inds;
      ijk_inds.push_back(newCluster[clusterId].i); ijk_inds.push_back(newCluster[clusterId].j); ijk_inds.push_back(newCluster[clusterId].k);

      vector<int> ijk_intersect_inds;
      ijk_intersect_inds.push_back(newCluster[clusterId].ij_intersect_loc);
      ijk_intersect_inds.push_back(newCluster[clusterId].jk_intersect_loc);
      ijk_intersect_inds.push_back(newCluster[clusterId].ki_intersect_loc);

      mplp.AddRegion(ijk_inds, ijk_intersect_inds);

      // TODO: log which clusters are chosen...

      nClustersAdded++;
   }

   return nClustersAdded;
}


/**
 * Main function implementing UAI 2012 cycle tightening algorithm.
 * 
 * method=1: use create_k_projection_graph
 * method=2: use create_expanded_projection_graph
 */
int TightenCycle(MPLPAlg & mplp, int nclus_to_add,  map<vector<int>, bool >& triplet_set, double & promised_bound, int method) {

   int nClustersAdded = 0;
   int nNewClusters;

   if (DEBUG_MODE) cout << "Finding the most violated cycle...." << endl;

   // This map allows us to quickly look up the edge weights
   map<pair<int, int>, double> projection_edge_weights;
   int num_projection_nodes;
   vector<vector<int> > projection_map;
   vector<int> projection_imap_var;  
   vector<vector<pair<int, double> > > projection_adjacency_list;
   vector<vector<int> > partition_imap;

   double* array_of_sij; int array_of_sij_size;

   // Define the projection graph and all edge weights
   if(method == 2)
      num_projection_nodes = create_expanded_projection_graph(mplp, projection_imap_var, projection_adjacency_list, projection_edge_weights, array_of_sij, array_of_sij_size, partition_imap);
   else if(method == 1)
      create_k_projection_graph(mplp, projection_map, num_projection_nodes, projection_imap_var, partition_imap, projection_edge_weights, projection_adjacency_list, array_of_sij, array_of_sij_size);
   else {
      cout << "ERROR: method not defined." << endl;
      return 0;
   }

   vector<list<int> > cycle_set;
   double optimal_R = find_optimal_R(projection_adjacency_list, array_of_sij, array_of_sij_size);
   if (DEBUG_MODE) cout << "R_optimal = " << optimal_R << endl;

   promised_bound = optimal_R;

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
   double total_time = (double)(end_time - start_time)/CLOCKS_PER_SEC;
   if (DEBUG_MODE) {
      printf(" -- FindCycles. Took %lg seconds\n", total_time);
   } 


   // Add all cycles that we've found to the relaxation!
   start_time = clock();
   for (int z = 0; z < cycle_set.size() && nClustersAdded < nclus_to_add; z++) {

      // Check to see if there are any duplicate nodes, and, if so, shortcut
      //    shortcut(cycle_set[z], projection_imap_var, projection_edge_weights, num_projection_nodes);

      // Output the cycle
      if (DEBUG_MODE){
         for (list<int>::iterator it=cycle_set[z].begin(); it!=cycle_set[z].end(); ++it) {
            cout << projection_imap_var[*it] << "(";
            vector<int> temp = partition_imap[*it];          
            for(int s=0; s < temp.size()-1; s++)
               cout << temp[s] << ",";
            cout << temp[temp.size()-1] << "), ";          
         }
         cout << endl;
      } 

      // Add cycle to the relaxation
      nClustersAdded += add_cycle(mplp, cycle_set[z], projection_imap_var, triplet_set, num_projection_nodes);
   }

   end_time = clock();
   total_time = (double)(end_time - start_time)/CLOCKS_PER_SEC;
   if (DEBUG_MODE) {
      printf(" -- shortcut + add_cycles. Took %lg seconds\n", total_time);
   } 

   delete_projection_graph(mplp.m_var_sizes.size(), projection_map, projection_imap_var, projection_adjacency_list, array_of_sij);
   return nClustersAdded;
}

#endif
