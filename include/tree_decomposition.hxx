#ifndef LP_MP_TREE_DECOMPOSITION_HXX
#define LP_MP_TREE_DECOMPOSITION_HXX

#include "LP_MP.h"
#include "serialization.hxx"

namespace LP_MP {

// factors are arranged in trees.
class LP_tree
{
   // make LP_MP::vector out of these
   //std::vector<tree_node> tree_nodes_;
   //std::vector<shared_factor> shared_factors_;

   friend class LP_with_trees;
public:
   // add messages from leaves to root upward
   void AddMessage(MessageTypeAdapter* m, Chirality c) // chirality denotes which factor is upper
   {
      tree_messages_.push_back(std::make_tuple(m, c));
   }

   // to do: give possiblity to specifiy root node
   void init()
   {
      static_assert(sizeof(REAL) == sizeof(double), "needed for SVM implementation");

      std::set<FactorTypeAdapter*> factor_set;
      for(auto tree_msg : tree_messages_) {
         auto* msg = std::get<0>(tree_msg);
         factor_set.insert(msg->GetLeftFactor());
         factor_set.insert(msg->GetRightFactor());
      }
      factors_.assign( factor_set.begin(), factor_set.end() );
      std::sort(factors_.begin(), factors_.end());
      assert(factors_.size() == tree_messages_.size() + 1);

      primal_size_in_bytes_ =  compute_primal_size_in_bytes();
      dual_size_in_bytes_ = compute_dual_size_in_bytes();
   }

   REAL compute_subgradient()
   {
      assert(factors_.size() == tree_messages_.size() + 1); // otherwise call init
      // send messages up the tree
      for(auto it = tree_messages_.begin(); it!= tree_messages_.end(); ++it) {
         auto* m = std::get<0>(*it);
         Chirality c = std::get<1>(*it);
         m->send_message_up(c);
      }
      // compute primal for topmost factor
      // also init primal for top factor, all other primals were initialized already by send_message_up
      REAL subgradient_value = 0.0;
      {
         auto* msg = std::get<0>(tree_messages_.back());
         Chirality c = std::get<1>(tree_messages_.back());
         if(c == Chirality::right) {
            // init primal for right factor!
            msg->GetRightFactorTypeAdapter()->init_primal();
            msg->GetRightFactorTypeAdapter()->MaximizePotentialAndComputePrimal();
            subgradient_value = msg->GetRightFactorTypeAdapter()->EvaluatePrimal();
         } else {
            assert(c == Chirality::left);
            msg->GetLeftFactorTypeAdapter()->init_primal(); 
            msg->GetLeftFactorTypeAdapter()->MaximizePotentialAndComputePrimal(); 
            subgradient_value = msg->GetLeftFactorTypeAdapter()->EvaluatePrimal();
         }
      }
      // track down optimal primal solution
      for(auto it = tree_messages_.rbegin(); it!= tree_messages_.rend(); ++it) {
         auto* msg = std::get<0>(*it);
         Chirality c = std::get<1>(*it);
         msg->track_solution_down(c);
         if(c == Chirality::right) {
            subgradient_value += msg->GetLeftFactorTypeAdapter()->EvaluatePrimal();
         } else {
            assert(c == Chirality::left);
            subgradient_value += msg->GetRightFactorTypeAdapter()->EvaluatePrimal();
         }
      } 

      // check if primal cost is equal to lower bound
      assert(primal_consistent());
      assert(std::abs(lower_bound() - primal_cost()) <= eps);
      assert(std::abs(subgradient_value - primal_cost()) <= eps);
      return subgradient_value;
   }

   template<typename VECTOR1>
   REAL compute_mapped_subgradient(VECTOR1& subgradient)
   {
      const REAL subgradient_value = compute_subgradient();
      std::vector<double> local_subgradient(dual_size(),0.0);
      // write primal solution into subgradient
      for(auto L : Lagrangean_factors_) {
         L.copy_fn(&local_subgradient[0]);
      }
      assert(mapping_.size() >= dual_size());
      for(INDEX i=0; i<dual_size(); ++i) {
         assert(mapping_[i] < subgradient.size());
         subgradient[ mapping_[i] ] += local_subgradient[i];
      } 
      return subgradient_value;
   }

   template<typename VECTOR>
   REAL compute_subgradient(VECTOR& subgradient, const REAL step_size)
   {
      const REAL subgradient_value = compute_subgradient();

      std::fill(subgradient.begin(), subgradient.end(), 0.0);
      // write primal solution into subgradient
      for(auto L : Lagrangean_factors_) {
         L.copy_fn(&subgradient[0]);
      }
      return subgradient_value;
   }

   bool primal_consistent() const 
   {
      for(auto it = tree_messages_.begin(); it!= tree_messages_.end(); ++it) {
         auto* m = std::get<0>(*it);
         if(!m->CheckPrimalConsistency()) {
            return false;
         }
      }
      return true; 
   }

   REAL primal_cost() const
   {
      REAL cost = 0.0;
      for(auto* f : factors_) {
         cost += f->EvaluatePrimal();
         if(cost == std::numeric_limits<REAL>::infinity()) {
            assert(false);
         }
      }
      return cost;
      /*
      for(auto it = tree_messages_.begin(); it!= tree_messages_.end(); ++it) {
         auto* m = std::get<0>(*it);
         Chirality c = std::get<1>(*it);
         if(c == Chirality::right) {
            cost += m->GetLeftFactorTypeAdapter()->EvaluatePrimal(); 
         } else {
            cost += m->GetRightFactorTypeAdapter()->EvaluatePrimal();
         }
      }
      if(std::get<1>(tree_messages_.back()) == Chirality::right) {
         cost += std::get<0>(tree_messages_.back())->GetRightFactorTypeAdapter()->EvaluatePrimal();
      } else {
         cost += std::get<0>(tree_messages_.back())->GetLeftFactorTypeAdapter()->EvaluatePrimal(); 
      }
      return cost;
      */
   }

   REAL lower_bound() const 
   {
      REAL lb = 0.0;
      for(auto* f : factors_) {
         lb += f->LowerBound();
      }
      return lb;
      /*
      for(auto it = tree_messages_.begin(); it!= tree_messages_.end(); ++it) {
         auto* m = std::get<0>(*it);
         Chirality c = std::get<1>(*it);
         if(c == Chirality::right) {
            lb += m->GetLeftFactorTypeAdapter()->LowerBound(); 
         } else {
            lb += m->GetRightFactorTypeAdapter()->LowerBound();
         }
      }
      if(std::get<1>(tree_messages_.back()) == Chirality::right) {
         lb += std::get<0>(tree_messages_.back())->GetRightFactorTypeAdapter()->LowerBound();
      } else {
         lb += std::get<0>(tree_messages_.back())->GetLeftFactorTypeAdapter()->LowerBound(); 
      }
      return lb; 
      */
   }

   template<typename FACTOR_TYPE>
   std::vector<FACTOR_TYPE*> get_factors() const
   {
      std::vector<FACTOR_TYPE*> factors;
      std::set<FACTOR_TYPE*> factor_present;
      for(auto& t : tree_messages_) {

         auto* left = std::get<0>(t)->GetLeftFactorTypeAdapter();
         auto* left_cast = dynamic_cast<FACTOR_TYPE*>(left);
         if(left_cast && factor_present.find(left_cast) == factor_present.end()) {
            factors.push_back(left_cast);
            factor_present.insert(left_cast);
         }

         auto* right = std::get<0>(t)->GetRightFactorTypeAdapter();
         auto* right_cast = dynamic_cast<FACTOR_TYPE*>(right);
         if(right_cast && factor_present.find(right_cast) == factor_present.end()) {
            factors.push_back(right_cast);
            factor_present.insert(right_cast);
         } 

      }
      return std::move(factors);
   }

   // find out necessary size for storing primal solution in archive
   INDEX compute_primal_size_in_bytes()
   {  
      // why not allocate_archive?
      INDEX size = 0;
      for(auto* f : Lagrangean_factors_) {
         size += f->primal_size_in_bytes();
      }
      return size;
   }
   INDEX primal_size_in_bytes()
   { 
      assert(primal_size_in_bytes_ == compute_primal_size_in_bytes());
      return primal_size_in_bytes_; 
   }

   void add_weights(const double* wi, const double scaling)
   {
      for(auto& L : Lagrangean_factors_) {
         L.serialize_Lagrangean(wi, scaling);
      }
   }

  // dual size of Lagrangeans connected to current tree
  INDEX compute_dual_size_in_bytes()
  {
     INDEX s = 0;
     for(auto& L : Lagrangean_factors_) {
        s += L.Lagrangean_vars_size_in_bytes();
     }
     return s;
  }

  INDEX dual_size_in_bytes()
  { 
     assert(dual_size_in_bytes_ == compute_dual_size_in_bytes());
     assert(dual_size_in_bytes_ % sizeof(REAL) == 0);
     return dual_size_in_bytes_; 
  }

  INDEX dual_size() { return dual_size_in_bytes() / sizeof(REAL); }
   
//protected:
   std::vector< std::tuple<MessageTypeAdapter*, Chirality>> tree_messages_; // messages forming a tree. Chirality says which side comprises the lower  factor
   std::vector<FactorTypeAdapter*> factors_;
   INDEX primal_size_in_bytes_;
   INDEX dual_size_in_bytes_;
   // subgradient information = primal solution to tree
   // for sending messages down we need to know to how many lower factors an upper factor is connected and then we need to average messages sent down appropriately.

   // Lagrangean variables come from copying factors
   struct Lagrangean_msg {
      FactorTypeAdapter* f;
      const INDEX no_trees; // number of trees in which factor is
      const INDEX pos; // factor is in pos-th tree that contains it
      INDEX Lagrangean_vars_offset; // at which offset are the Lagrangean variables for factor f stored?
      INDEX dual_size_; // size as multiples of REAL
      // replace by function into factor (serialize dual)

      INDEX dual_size() {
         assert(dual_size_ == f->dual_size());
         return dual_size_; 
      }
      INDEX dual_size_in_bytes() 
      { 
         return dual_size()*sizeof(REAL); 
      }

      INDEX global_offset(const INDEX offset, const INDEX i, const INDEX j) // offset for Lagrangean variables (i,j)
      {
         //assert(no_trees <= 2);
         assert(i < j && j < no_trees);
         return offset + (i*no_trees - i*(i+1)/2 + (j-i-1))*dual_size(); 
      }

      INDEX offset(const INDEX i) // offset for Lagrangean variables (i,j)
      {
         assert(i < no_trees);
         const INDEX o = Lagrangean_vars_offset + i*dual_size();
         return o;
      }

      INDEX Lagrangean_vars_size() 
      {
         assert(no_trees > 1);
         return (no_trees-1)*dual_size(); //((no_trees-1)*no_trees)/2*dual_size;
      }

      INDEX Lagrangean_vars_size_in_bytes() 
      {
         return (no_trees-1)*dual_size_in_bytes(); //((no_trees-1)*no_trees)/2*dual_size;
      }

      // +1 is for loading, -1 is for unloading
      void serialize_Lagrangean(const double* wi, const double scaling)
      {
         for(INDEX i=0; i<pos; ++i) {
            const double* w = wi+offset(i);
            serialization_archive ar(w, Lagrangean_vars_size_in_bytes());
            addition_archive l_ar(ar, 1.0*scaling);
            f->serialize_dual(l_ar);
            ar.release_memory();
         }
         for(INDEX i=pos+1; i<no_trees; ++i) {
            //double* w = wi+offset(pos, i);
            const double* w = wi+offset(i-1);
            serialization_archive ar(w, Lagrangean_vars_size_in_bytes());
            addition_archive l_ar(ar, -1.0*scaling);
            f->serialize_dual(l_ar);
            ar.release_memory();
         }
      } 

      void copy_fn(double* wi)
      {
         for(INDEX i=0; i<pos; ++i) {
            double* w = wi+offset(i);
            f->subgradient(w, +1.0);
         } 
         for(INDEX i=pos+1; i<no_trees; ++i) {
            //double* w = wi+offset(pos,i);
            double* w = wi+offset(i-1);
            f->subgradient(w, -1.0);
         }
      }

      REAL dot_product_fn(double* wi)
      {
         REAL d = 0.0;
         for(INDEX i=0; i<pos; ++i) {
            double* w = wi+offset(i);
            d += f->dot_product(w);
         } 
         for(INDEX i=pos+1; i<no_trees; ++i) {
            double* w = wi+offset(i-1);
            d -= f->dot_product(w);
         }
         return d;
      }
   };

   const std::vector<int>& mapping() const { return mapping_; }
   std::vector<int>& mapping() { return mapping_; }

   // we have Lagrangean variables for every pair of factors that are identical but occur in different trees, hence can be indexed by factor and their pos
   // Lagrangean variables are stored contiguosly for each factor in lexicographic order of pos

   std::vector<Lagrangean_msg> Lagrangean_factors_;
   INDEX subgradient_size;
   std::vector<int> mapping_;

   //serialization_archive potentials_archive_; // possibly only Lagrangean factors need to be stored here
};

// do zrobienia: templatize base class
// rename to LP_Frank_Wolfe
class LP_with_trees : public LP
{
public:
   using LP::LP;

   ~LP_with_trees()
   {
      // delete copies of factors and redirect messages back to original factors
      assert(false);
   }

   void add_tree(LP_tree& t)
   { 
      trees_.push_back(t);
   }

   // find out, which factors are shared between trees and add Lagrangean multipliers for them.
   void Begin()
   {
      LP::Begin();

      // first, go over all Lagrangean factors in each tree and count how often factor is shared
      struct Lagrangean_counting {
         INDEX no_trees = 1; // in how many trees is factor,
         INDEX position = 0; // counter for enumerating in which position (i.e. in how many trees was factor already observed)
         INDEX offset = 0; // offset into dual weights
      };
      std::unordered_map<FactorTypeAdapter*, Lagrangean_counting > Lagrangean_factors;
      for(auto& t : trees_) {
         for(auto* f : t.factors_) {
            auto it = Lagrangean_factors.find(f);
            if(it == Lagrangean_factors.end()) {
               Lagrangean_factors.insert(std::make_pair(f,Lagrangean_counting()));
            } else {
               it->second.no_trees++;
            }
         }
      }
      assert(Lagrangean_factors.size() == this->f_.size()); // otherwise not all factors are covered by trees

      // now set Lagrangean_msg in each tree
      for(auto& t : trees_) {
         for(auto* f : t.factors_) {
            auto it = Lagrangean_factors.find(f);
            const INDEX no_occurences = it->second.no_trees;
            const INDEX pos = it->second.position;
            it->second.position++;
            if(no_occurences > 1) {
               t.Lagrangean_factors_.push_back({f, no_occurences, pos, 0, f->dual_size()});
            }
         }
      }

      // count needed size for Lagrangean variables and index into array for each shared factor
      Lagrangean_vars_size_ = 0;
      for(auto& it : Lagrangean_factors) {
         const INDEX no_occurences = it.second.no_trees;
         if(no_occurences > 1) {
            assert(it.second.no_trees == it.second.position); // number of occurences and positions correctly counted
            it.second.offset = Lagrangean_vars_size_;
            const INDEX delta_offset = ((no_occurences-1)*no_occurences)/2 * it.first->dual_size();
            assert(delta_offset > 0);
            Lagrangean_vars_size_ += delta_offset;
         }
      }

      // set offsets into Lagrangean vars for each shared factor in tree
      // we need this ugly construction of first putting global offsets and only later local ones because of factor splitting, otherwise global offset could be looked up in Lagrangean_factors.
      for(auto& t : trees_) {
         for(auto& L : t.Lagrangean_factors_) {
            auto it = Lagrangean_factors.find(L.f);
            assert(it->second.no_trees > 1);
            const INDEX offset = it->second.offset;
            L.Lagrangean_vars_offset = offset;
         }
      }

      // make copies of shared factors, their number being equal to the number of times they are shared. Redirect links to factors in relevant messages
      for(auto& t : trees_) {
         std::unordered_map<FactorTypeAdapter*, FactorTypeAdapter*> factor_mapping; // original to copied factor
         for(auto& L : t.Lagrangean_factors_) {
            const INDEX no_occurences = L.no_trees;
            assert(no_occurences > 1);
            const INDEX pos = L.pos;
            if(pos > 0) {
               auto* f_copy = L.f->clone(); // do zrobienia: possibly not all pointers to messages have to be cloned as well
               f_copy->divide(no_occurences);
               factor_mapping.insert(std::make_pair(L.f, f_copy));
               L.f = f_copy;
            }
         }
         // redirect links from messages in trees that are directed to current factor
         for(auto tree_msg : t.tree_messages_) {
            auto* m = std::get<0>(tree_msg);
            auto* left = m->GetLeftFactor();
            auto* right = m->GetRightFactor();
            if(factor_mapping.find(left) != factor_mapping.end()) {
               auto* left_copy = factor_mapping.find(left)->second;
               m->SetLeftFactor(left_copy);
            }
            if(factor_mapping.find(right) != factor_mapping.end()) {
               auto* right_copy = factor_mapping.find(right)->second;
               m->SetRightFactor(right_copy);
            }
         }
         // search for factor in tree and change it as well
         for(auto& f : t.factors_) {
            if(factor_mapping.find(f) != factor_mapping.end()) {
               auto* f_copy = factor_mapping.find(f)->second;
               f = f_copy;
            } 
         }
      }

      // divide cost of non-copied factor
      for(auto& t : trees_) {
         for(auto& L : t.Lagrangean_factors_) {
            if(L.no_trees > 1 && L.pos == 0) {
               L.f->divide(L.no_trees); 
            }
         }
      }

      // set primal and dual size for tree
      for(auto& t : trees_) {
         t.init();
      }

      // construct mapping from Lagrangean variables of each tree to whole Lagrangean variables
      //mapping_.reserve(trees_.size());
      for(auto& t : trees_) {
         //std::cout << "tree " << mapping_.size() << "\n";
         INDEX local_offset = 0;
         std::vector<int> m;
         m.reserve(t.dual_size()+1); // +1 because of requirement in SVM that last entry is equal to total size of Lagrangean variables
         for(auto& L : t.Lagrangean_factors_) {

            // change global offset currently stored into local offset
            const INDEX global_offset = L.Lagrangean_vars_offset; 
            L.Lagrangean_vars_offset = m.size();

            //std::cout << L.f << " = " << L.f->dual_size() << "\n";
            assert(L.no_trees > 1);
            const INDEX pos = L.pos;
            for(INDEX i=0; i<pos; ++i) {
               //std::cout << "before pos: ";
               for(INDEX j=0; j<L.dual_size(); ++j) {
                  m.push_back(L.global_offset(global_offset, i, pos) + j);
                  //std::cout << m.back() << ",";
                  assert(m.back() < Lagrangean_vars_size_);
               }
            }
            //std::cout << "\n";
            for(INDEX i=pos+1; i<L.no_trees; ++i) {
               //std::cout << "after pos: ";
               for(INDEX j=0; j<L.dual_size(); ++j) {
                  m.push_back(L.global_offset(global_offset, pos, i) + j);
                  //std::cout << m.back() << ",";
                  assert(m.back() < Lagrangean_vars_size_);
               }
            }
            //std::cout << "\n";
            //assert(local_offset == m.size());
            //local_offset += L.Lagrangean_vars_size();
         }
         assert(m.size() == t.dual_size());
         m.push_back(Lagrangean_vars_size_);
         //assert(local_offset == m.size());
         t.mapping_ = m;
      } 

      // check map validity: each entry in m (except last one) must occur exactly twice
      assert(mapping_valid()); 
   }

   bool mapping_valid() const
   {
      std::vector<INDEX> mapping_count(Lagrangean_vars_size_,0);
      for(const auto& t : trees_) {
         for(INDEX i=0; i<t.mapping().size()-1; ++i) {
            mapping_count[t.mapping()[i]]++;
         }
      }
      for(auto i : mapping_count) {
         if(i != 2) {
            return false;
         }
      }
      return true;
   }

   INDEX no_Lagrangean_vars() const { return Lagrangean_vars_size_; }

   
   void ComputeForwardPassAndPrimal(const INDEX iteration) 
   {
      assert(false);
   }

   void ComputeBackwardPassAndPrimal(const INDEX iteration) 
   {
      assert(false); 
   }

   virtual void ComputePass(const INDEX iteration) = 0;

   // lower bound must be computed over all cloned factors as well
   REAL LowerBound() const
   {
      REAL lb = 0.0;
      for(auto& t : trees_) {
         lb += t.lower_bound();
      }
      return lb;
   }

   void add_weights(const double* w, const REAL scaling) 
   {
      for(INDEX i=0; i<trees_.size(); ++i) {
         auto& tree = trees_[i];
         const auto& m = tree.mapping();
         std::vector<double> local_weights;
         local_weights.reserve(m.size()-1);
         for(INDEX idx=0; idx<m.size()-1; ++idx) {
            local_weights.push_back(w[m[idx]]);
         }
         tree.add_weights(&local_weights[0], scaling);
      }
   }


protected:
   std::vector<LP_tree> trees_;
   INDEX Lagrangean_vars_size_;
};

// perform subgradient ascent with Polyak's step size with estimated optimum
class LP_subgradient_ascent : public LP_with_trees {
public:
   using LP_with_trees::LP_with_trees;

   void Begin()
   {
      best_lower_bound = -std::numeric_limits<REAL>::infinity();
      LP_with_trees::Begin(); 
   }

   void ComputePass(const INDEX iteration)
   {
      REAL current_lower_bound = 0.0;
      std::vector<REAL> subgradient(this->no_Lagrangean_vars(), 0.0);
      for(INDEX i=0; i<trees_.size(); ++i) {
         current_lower_bound += trees_[i].compute_mapped_subgradient(subgradient); // note that mapping has one extra component!
         //trees_[i].primal_cost();
      }
      best_lower_bound = std::max(current_lower_bound, best_lower_bound);
      assert(std::find_if(subgradient.begin(), subgradient.end(), [](auto x) { return x != 0.0 && x != 1.0 && x != -1.0; }) == subgradient.end());
      const REAL subgradient_one_norm = std::accumulate(subgradient.begin(), subgradient.end(), 0.0, [=](REAL s, REAL x) { return s + std::abs(x); });

      const REAL step_size = (best_lower_bound - current_lower_bound + subgradient.size())/(10.0 + iteration) / subgradient_one_norm;

      std::cout << "stepsize = " << step_size << ", absolute value of subgradient = " << subgradient_one_norm << "\n";
      add_weights(&subgradient[0], step_size);
   }

private:
   REAL best_lower_bound;
   REAL step_size;
};

} // end namespace LP_MP

#endif // LP_MP_TREE_DECOMPOSITION_HXX
