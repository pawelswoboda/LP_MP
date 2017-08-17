#ifndef LP_MP_MAIN
#define LP_MP_MAIN

#include "config.hxx"
#include <vector>
#include <valarray>
#include <map>
#include <iostream>
#include <numeric>
#include <algorithm>
#include <functional>
#include <utility>
#include <limits>
#include <exception>
#include "template_utilities.hxx"
#include <assert.h>
#include "topological_sort.hxx"
#include <memory>
#include <iterator>
#include "primal_solution_storage.hxx"
#include "lp_interface/lp_interface.h"
#include "two_dimensional_variable_array.hxx"
#include <thread>
#include <future>
#include "memory_allocator.hxx"
#include "cereal/archives/binary.hpp"
#include "serialization.hxx"
#include "sat_interface.hxx"
#include "tclap/CmdLine.h"

#ifdef LP_MP_PARALLEL
#include <omp.h>
#endif

namespace LP_MP {

// forward declaration
class MessageTypeAdapter;
class MessageIterator;

using weight_vector = two_dim_variable_array<REAL>::ArrayAccessObject;

// pure virtual base class for factor container used by LP class
class FactorTypeAdapter
{
public:
   virtual ~FactorTypeAdapter() {}
   virtual FactorTypeAdapter* clone() const = 0;
   virtual void UpdateFactor(const weight_vector& omega) = 0;
   virtual void UpdateFactorPrimal(const weight_vector& omega, const INDEX iteration) = 0;
   virtual void UpdateFactorSAT(const weight_vector& omega, const REAL th, sat_var begin, sat_vec<sat_literal>& assumptions) = 0;
   //virtual void convert_primal(Glucose::SimpSolver&, sat_var) = 0; // this is not nice: the solver should be templatized
   //virtual void convert_primal(CMSat::SATSolver&, sat_var) = 0; // this is not nice: the solver should be templatized
   virtual void construct_sat_clauses(LGL*) = 0;
   virtual void convert_primal(LGL*, sat_var) = 0; // this is not nice: the solver should be templatized
   virtual bool FactorUpdated() const = 0; // does calling UpdateFactor do anything? If no, it need not be called while in ComputePass, saving time.
   // to do: remove both
   virtual INDEX size() const = 0;
   virtual INDEX PrimalSize() const = 0;
   MessageIterator begin(); 
   MessageIterator end();
   virtual INDEX GetNoMessages() const = 0;
   virtual INDEX no_send_messages() const = 0; // to do: use this function instead of ComputeSendFactorConnection in weight computation whenever possible
   virtual MessageTypeAdapter* GetMessage(const INDEX n) const = 0;
   virtual FactorTypeAdapter* GetConnectedFactor(const INDEX i) const = 0;
   virtual bool CanSendMessage(const INDEX i) const = 0;
   virtual REAL LowerBound() const = 0;
   virtual void init_primal() = 0;
   virtual void MaximizePotentialAndComputePrimal() = 0;

   // for use in tree decomposition:
   // for writing primal solution into subgradient
   // return value is size of subgradient
   virtual INDEX subgradient(double* w, const REAL sign) = 0;
   virtual REAL dot_product(double* w) = 0;

   // for reading reparametrization/labeling out of factor
   virtual void serialize_dual(save_archive&) = 0;
   virtual void serialize_primal(save_archive&) = 0;
   // for writing reparametrization/labeling into factor
   virtual void serialize_dual(load_archive&) = 0;
   virtual void serialize_primal(load_archive&) = 0;
   // for determining size of archive
   virtual void serialize_dual(allocate_archive&) = 0;
   virtual void serialize_primal(allocate_archive&) = 0;
   // for adding weights in Frank Wolfe algorithm
   virtual void serialize_dual(addition_archive&) = 0;

   virtual void divide(const REAL val) = 0; // divide potential by value

   virtual INDEX dual_size() = 0;
   virtual INDEX dual_size_in_bytes() = 0;
   virtual INDEX primal_size_in_bytes() = 0;

   // the offset in the primal storage
   // to do: remove
   virtual void SetPrimalOffset(const INDEX) = 0; // do zrobienia: delete
   virtual INDEX GetPrimalOffset() const = 0; // do zrobienia: delete
   
   virtual void SetAuxOffset(const INDEX n) = 0; // do zrobienia: delete
   virtual INDEX GetAuxOffset() const = 0; // do zrobienia: delete
   
   // do zrobienia: this function is not needed. Evaluation can be performed automatically
   virtual REAL EvaluatePrimal() const = 0;

   // for the LP interface
   virtual INDEX GetNumberOfAuxVariables() const = 0;
   virtual void CreateConstraints(LpInterfaceAdapter* lpInterface) const = 0;
   virtual void ReduceLp(LpInterfaceAdapter* lpInterface) const = 0;
};

class MessageTypeAdapter
{
public:
   virtual ~MessageTypeAdapter() {}
   virtual MessageTypeAdapter* clone(FactorTypeAdapter* l, FactorTypeAdapter* r) const = 0;
   virtual FactorTypeAdapter* GetLeftFactor() const = 0;
   virtual FactorTypeAdapter* GetRightFactor() const = 0;
   virtual void SetLeftFactor(FactorTypeAdapter*) = 0;
   virtual void SetRightFactor(FactorTypeAdapter*) = 0;
   //virtual bool CheckPrimalConsistency(typename PrimalSolutionStorage::Element left, typename PrimalSolutionStorage::Element right) const = 0;
   virtual bool SendsMessageToLeft() const = 0;
   virtual bool SendsMessageToRight() const = 0;
   virtual bool ReceivesMessageFromLeft() const = 0;
   virtual bool ReceivesMessageFromRight() const = 0;
   virtual bool CheckPrimalConsistency() const = 0;
   virtual FactorTypeAdapter* GetRightFactorTypeAdapter() const = 0;
   virtual FactorTypeAdapter* GetLeftFactorTypeAdapter() const = 0;

   virtual void send_message_up(Chirality c) = 0;
   virtual void track_solution_down(Chirality c) = 0;
   
   // Also true, if SendMessagesTo{Left|Right} is active. Used for weight computation. Disregard message in weight computation if it does not send messages at all
   // do zrobienia: throw them out again
   //virtual bool CanSendMessageToLeft() const = 0;
   //virtual bool CanSendMessageToRight() const = 0;

   virtual void construct_sat_clauses(LGL*, sat_var, sat_var) = 0;
   // for the LP interface
   virtual void CreateConstraints(LpInterfaceAdapter* lpInterface) = 0;
};

// primitive iterator class. Access may be slow. A more direct implementation would be more complicated, though.
// used in main LP class
class MessageIterator //: public std::iterator<random_access_iterator>
{
public:
   MessageIterator(FactorTypeAdapter* const factor, const INDEX msg_idx) : factor_(factor), msg_idx_(msg_idx) {}
   MessageIterator* operator++() { ++msg_idx_; return this; } // do zrobienia: is * correct in return value?
   bool operator==(const MessageIterator& rhs) const { return (factor_ == rhs.factor_ && msg_idx_ == rhs.msg_idx_); }
   bool operator!=(const MessageIterator& rhs) const { return !operator==(rhs); }
   MessageTypeAdapter& operator*() const { return *(factor_->GetMessage(msg_idx_)); }
   MessageTypeAdapter& operator->() const { return *(factor_->GetMessage(msg_idx_)); }
   FactorTypeAdapter* GetConnectedFactor() const { return factor_->GetConnectedFactor(msg_idx_); }
   bool CanSendMessage() const  { return factor_->CanSendMessage(msg_idx_); }
private:
   FactorTypeAdapter* const factor_;
   INDEX msg_idx_;
};

inline MessageIterator FactorTypeAdapter::begin() { return MessageIterator(this,0); }
inline MessageIterator FactorTypeAdapter::end()  { return MessageIterator(this, GetNoMessages()); }


class LP {
public:
   LP(TCLAP::CmdLine& cmd) 
#ifdef LP_MP_PARALLEL
      : num_lp_threads_arg_("","numLpThreads","number of threads for message passing, default = 1",false,1,&positiveIntegerConstraint,cmd)
#endif
   {}

   ~LP() 
   {
      for(INDEX i=0; i<m_.size(); i++) { delete m_[i]; }
      for(INDEX i=0; i<f_.size(); i++) { delete f_[i]; }
   }

   // make a deep copy of factors and messages. Adjust pointers to messages and factors
   LP(LP& o) // no const because of o.num_lp_threads_arg_.getValue() not being const!
#ifdef LP_MP_PARALLEL
      : num_lp_threads_arg_("","numLpThreads","number of threads for message passing, default = 1",false,o.num_lp_threads_arg_.getValue(),&positiveIntegerConstraint)
#endif
   {
      f_.reserve(o.f_.size());
      assert(false);
      std::map<FactorTypeAdapter*, FactorTypeAdapter*> factor_map; // translate addresses from o's factors to this' factors
      f_.reserve(o.f_.size());
      for(auto* f : o.f_) {
         auto* clone = f->clone();
         this->AddFactor(clone);
         factor_map.insert(std::make_pair(f, clone));
      }
      m_.reserve(o.m_.size());
      for(auto* m : o.m_) {
         auto* left = m->GetLeftFactorTypeAdapter();
         auto* right = m->GetRightFactorTypeAdapter();
         auto* left_clone = factor_map[left];
         auto* right_clone = factor_map[right];
         auto* clone = m->clone(left_clone, right_clone);
         this->AddMessage(m);
      }

      ordering_valid_ = o.ordering_valid_;
      omega_anisotropic_valid_ = o.omega_anisotropic_valid_ ;
      omega_anisotropic2_valid_ = o.omega_anisotropic2_valid_;
      omega_isotropic_valid_ = o.omega_isotropic_valid_ ;
      omega_isotropic_damped_valid_ = o.omega_isotropic_damped_valid_ ;
      omega_mixed_valid_ = o.omega_mixed_valid_;

      omegaForwardAnisotropic_ = o.omegaForwardAnisotropic_; 
      omegaBackwardAnisotropic_ = o.omegaBackwardAnisotropic_;
      omegaForwardAnisotropic2_ = o.omegaForwardAnisotropic_; 
      omegaBackwardAnisotropic2_ = o.omegaBackwardAnisotropic_;
      omegaForwardIsotropic_ = o.omegaForwardIsotropic_; 
      omegaBackwardIsotropic_ = o.omegaBackwardIsotropic_;
      omegaForwardIsotropicDamped_ = o.omegaForwardIsotropicDamped_; 
      omegaBackwardIsotropicDamped_ = o.omegaBackwardIsotropicDamped_;
      omegaForwardMixed_ = o.omegaForwardMixed_; 
      omegaBackwardMixed_ = o.omegaBackwardMixed_;

      forwardOrdering_.reserve(o.forwardOrdering_.size());
      for(auto* f : o.forwardOrdering_) {
        forwardOrdering_.push_back( factor_map[f] );
      }

      backwardOrdering_.reserve(o.backwardOrdering_.size());
      for(auto* f : o.backwardOrdering_) {
        backwardOrdering_.push_back( factor_map[f] );
      }

      forwardUpdateOrdering_.reserve(o.forwardUpdateOrdering_.size());
      for(auto* f : o.forwardUpdateOrdering_) {
        forwardUpdateOrdering_.push_back( factor_map[f] );
      }

      backwardUpdateOrdering_.reserve(o.backwardUpdateOrdering_.size());
      for(auto* f : o.backwardUpdateOrdering_) {
        backwardUpdateOrdering_.push_back( factor_map[f] );
      }

      forward_pass_factor_rel_.reserve(o.forward_pass_factor_rel_.size());
      for(auto f : o.forward_pass_factor_rel_) {
        forward_pass_factor_rel_.push_back( std::make_pair(factor_map[f.first], factor_map[f.second]) );
      }

      backward_pass_factor_rel_.reserve(o.backward_pass_factor_rel_.size());
      for(auto f : o.backward_pass_factor_rel_) {
        backward_pass_factor_rel_.push_back( std::make_pair(factor_map[f.first], factor_map[f.second]) );
      }

   }
   /*
   std::vector<FactorTypeAdapter*> f_; // note that here the factors are stored in the original order they were given. They will be output in this order as well, e.g. by problemDecomposition
   std::vector<MessageTypeAdapter*> m_;
   std::vector<FactorTypeAdapter*> forwardOrdering_, backwardOrdering_; // separate forward and backward ordering are not needed: Just store factorOrdering_ and generate forward order by begin() and backward order by rbegin().
   std::vector<FactorTypeAdapter*> forwardUpdateOrdering_, backwardUpdateOrdering_; // like forwardOrdering_, but includes only those factors where UpdateFactor actually does something
   std::vector<std::vector<REAL> > omegaForward_, omegaBackward_;
   std::vector<std::pair<FactorTypeAdapter*, FactorTypeAdapter*> > forward_pass_factor_rel_, backward_pass_factor_rel_; // factor ordering relations. First factor must come before second factor. factorRel_ must describe a DAG

   
   //REAL bestLowerBound_ = -std::numeric_limits<REAL>::infinity();
   //REAL currentLowerBound_ = -std::numeric_limits<REAL>::infinity();

   LPReparametrizationMode repamMode_ = LPReparametrizationMode::Undefined;
   */

   virtual INDEX AddFactor(FactorTypeAdapter* f)
   {
      set_flags_dirty();
      assert(factor_address_to_index_.size() == f_.size());
      f_.push_back(f);

      assert(factor_address_to_index_.find(f) == factor_address_to_index_.end());
      factor_address_to_index_.insert(std::make_pair(f,f_.size()-1));
      assert(factor_address_to_index_.find(f)->second == f_.size()-1);

      return f_.size() - 1;
   }

   INDEX GetNumberOfFactors() const { return f_.size(); }
   FactorTypeAdapter* GetFactor(const INDEX i) const { return f_[i]; }
   INDEX size() const { INDEX size=0; for(auto* f : f_) { size += f->size(); } return size; }

   virtual INDEX AddMessage(MessageTypeAdapter* m)
   {
      set_flags_dirty();
      m_.push_back(m);
      // do zrobienia: check whether left and right factors are in f_

      //////////////////////////////////////////////////////
      // check whether left and right factor has all different factors connected to it, likewise with the right one
      /*
         std::vector<FactorTypeAdapter*> fc;
         auto f_left = m->GetLeftFactor();
         for(auto mIt=f_left->begin(); mIt!=f_left->end(); ++mIt) {
         fc.push_back( &*mIt ); // dereference iterator to factor and then take address of factor
         }
         assert( HasUniqueValues(fc) );
         fc.clear();
         auto f_right = m->GetRightFactor();
         for(auto mIt=f_right->begin(); mIt!=f_right->end(); ++mIt) {
         fc.push_back( &*mIt ); // dereference iterator to factor and then take address of factor
         }
         assert( HasUniqueValues(fc) );
         */
      ////////////////////////////////////////////////////

      return m_.size() - 1;
   }

   MessageTypeAdapter* GetMessage(const INDEX i) const { return m_[i]; }
   INDEX GetNumberOfMessages() const { return m_.size(); }
   void AddFactorRelation(FactorTypeAdapter* f1, FactorTypeAdapter* f2) // indicate that factor f1 comes before factor f2
   {
      ForwardPassFactorRelation(f1,f2);
      BackwardPassFactorRelation(f2,f1);
   }

   void ForwardPassFactorRelation(FactorTypeAdapter* f1, FactorTypeAdapter* f2) { set_flags_dirty(); assert(f1!=f2); forward_pass_factor_rel_.push_back({f1,f2}); }
   void BackwardPassFactorRelation(FactorTypeAdapter* f1, FactorTypeAdapter* f2) { set_flags_dirty(); assert(f1!=f2); backward_pass_factor_rel_.push_back({f1,f2}); }

   // initially select branching factor from among those that can be branched on at all.
   FactorTypeAdapter* select_branching_factor()
   {
      assert(false);
      return nullptr;
   }
   template<typename ITERATOR>
   FactorTypeAdapter* select_branching_factor(ITERATOR factor_begin, ITERATOR factor_end)
   {
      // go over all given factors and select the one with largest difference between lower bound and primal cost;
      FactorTypeAdapter* f;
      REAL max_diff = -std::numeric_limits<REAL>::infinity();
      for(; factor_begin!=factor_end; ++factor_begin) {
         const REAL diff = (*factor_begin)->EvaluatePrimal() - factor_begin->LowerBound();
         if(diff > max_diff) {
            f = *factor_begin;
            max_diff = diff;
         } 
      }
      return f; 
   }

   void Begin(); // must be called after all messages and factors have been added

   void CalculatePrimalOffsets()
   {
      INDEX primalOffset = 0;
      for(auto* f : f_) {
         f->SetPrimalOffset(primalOffset);
         primalOffset += f->PrimalSize();
      }
   }
   void SortFactors(
         const std::vector<std::pair<FactorTypeAdapter*, FactorTypeAdapter*>>& factor_rel,
         std::vector<FactorTypeAdapter*>& ordering,
         std::vector<FactorTypeAdapter*>& update_ordering,
         std::vector<INDEX>& f_sorted
         ) 
   {
      // assume that factorRel_ describe a DAG. Compute topological sorting
      Topological_Sort::Graph g(f_.size());

      //std::map<FactorTypeAdapter*,INDEX> factorToIndex; // possibly do it with a hash_map for speed
      //std::map<INDEX,FactorTypeAdapter*> indexToFactor; // do zrobienia: need not be map, oculd be vector!
      //BuildIndexMaps(f_.begin(), f_.end(),factorToIndex,indexToFactor);

      for(auto fRelIt=factor_rel.begin(); fRelIt!=factor_rel.end(); fRelIt++) {
         // do zrobienia: why do these asserts fail?
         assert(factor_address_to_index_.find(fRelIt->first ) != factor_address_to_index_.end());
         assert(factor_address_to_index_.find(fRelIt->second) != factor_address_to_index_.end());
         INDEX f1 = factor_address_to_index_[fRelIt->first];
         INDEX f2 = factor_address_to_index_[fRelIt->second];
         g.addEdge(f1,f2);
      }

      f_sorted = g.topologicalSort();
      //std::vector<INDEX> sortedIndices = g.topologicalSort();
      assert(f_sorted.size() == f_.size());

      std::vector<FactorTypeAdapter*> fSorted;
      fSorted.reserve(f_.size());
      for(INDEX i=0; i<f_sorted.size(); i++) {
         fSorted.push_back( f_[ f_sorted[i] ] );//indexToFactor[sortedIndices[i]] );
      }
      assert(fSorted.size() == f_.size());
      assert(HasUniqueValues(fSorted));
      ordering  = fSorted;
      update_ordering.clear();
      for(auto f : ordering) {
         if(f->FactorUpdated()) {
            update_ordering.push_back(f);
         }
      }
      // check whether sorting was successful
      /*
      std::map<FactorTypeAdapter*, INDEX> factorToIndexSorted;
      std::map<INDEX, FactorTypeAdapter*> indexToFactorSorted;
      BuildIndexMaps(ordering.begin(), ordering.end(), factorToIndexSorted, indexToFactorSorted);
      for(auto rel : factor_rel) {
         const INDEX index_left = factorToIndexSorted[ std::get<0>(rel) ];
         const INDEX index_right = factorToIndexSorted[ std::get<1>(rel) ];
         assert(index_left < index_right);
      }
      */
   }
   void SortFactors()
   {
      if(ordering_valid_) { return; }
      ordering_valid_ = true;

#pragma omp parallel sections
      {
#pragma omp section
         SortFactors(forward_pass_factor_rel_, forwardOrdering_, forwardUpdateOrdering_, f_forward_sorted_);
#pragma omp section
         SortFactors(backward_pass_factor_rel_, backwardOrdering_, backwardUpdateOrdering_, f_backward_sorted_);
      }
   }

   template<typename FACTOR_ITERATOR>
      std::vector<std::vector<FactorTypeAdapter*> > ComputeFactorConnection(FACTOR_ITERATOR factorIt, FACTOR_ITERATOR factorEndIt);
   template<typename FACTOR_ITERATOR>
      std::vector<std::vector<FactorTypeAdapter*> > ComputeSendFactorConnection(FACTOR_ITERATOR factorIt, FACTOR_ITERATOR factorEndIt);

   //void ComputeWeights(const LPReparametrizationMode m);
   void set_reparametrization(const LPReparametrizationMode r) 
   {
      repamMode_ = r;
   }
   void ComputeAnisotropicWeights();
   template<typename FACTOR_ITERATOR, typename FACTOR_SORT_ITERATOR>
   void ComputeAnisotropicWeights(FACTOR_ITERATOR factorIt, FACTOR_ITERATOR factorItEnd, FACTOR_SORT_ITERATOR factor_sort_begin, FACTOR_SORT_ITERATOR factor_sort_end, two_dim_variable_array<REAL>& omega); 

   void ComputeAnisotropicWeights2();
   template<typename FACTOR_ITERATOR, typename FACTOR_SORT_ITERATOR>
   void ComputeAnisotropicWeights2(FACTOR_ITERATOR factorIt, FACTOR_ITERATOR factorItEnd, FACTOR_SORT_ITERATOR factor_sort_begin, FACTOR_SORT_ITERATOR factor_sort_end, two_dim_variable_array<REAL>& omega); 

   void ComputeUniformWeights();

   void ComputeDampedUniformWeights();
   template<typename FACTOR_ITERATOR>
   void ComputeUniformWeights(FACTOR_ITERATOR factorIt, FACTOR_ITERATOR factorEndIt, two_dim_variable_array<REAL>& omega, const REAL leave_weight); // do zrobienia: rename to isotropic weights

   void ComputeMixedWeights(const two_dim_variable_array<REAL>& omega_anisotropic, const two_dim_variable_array<REAL>& omega_damped_uniform, two_dim_variable_array<REAL>& omega); 
   void ComputeMixedWeights()
   {
      if(!omega_mixed_valid_) {
#pragma omp parallel sections
         {
#pragma omp section
            ComputeDampedUniformWeights();
#pragma omp section
            ComputeAnisotropicWeights();
         }
         omega_mixed_valid_ = true;
#pragma omp sections
         {
#pragma omp section
            ComputeMixedWeights(omegaForwardAnisotropic_, omegaForwardIsotropicDamped_, omegaForwardMixed_);
#pragma omp section
            ComputeMixedWeights(omegaBackwardAnisotropic_, omegaBackwardIsotropicDamped_, omegaBackwardMixed_);
         }
      }
   }

   double LowerBound() const
   {
      double lb = 0.0;
#pragma omp parallel for reduction(+:lb)
      for(INDEX i=0; i<f_.size(); ++i) {
         lb += f_[i]->LowerBound();
         assert( f_[i]->LowerBound() > -10000000.0);
         assert(std::isfinite(lb));
      }
      return lb;
   }

   bool CheckPrimalConsistency() const;

   template<typename FACTOR_ITERATOR>
   double EvaluatePrimal(cereal::BinaryInputArchive& primal, FACTOR_ITERATOR factor_begin, FACTOR_ITERATOR factor_end) // read in primal solution from primal archive and evaluate cost
   {
      assert(std::distance(factor_begin, factor_end) == f_.size());
      assert(false); 
      return 0.0;
   }
   double EvaluatePrimal() {
      const bool consistent = CheckPrimalConsistency();
      if(consistent == false) return std::numeric_limits<REAL>::infinity();

      double cost = 0.0;
#pragma omp parallel for reduction(+:cost)
      for(INDEX i=0; i<f_.size(); ++i) {
         cost += f_[i]->EvaluatePrimal();
      }

      if(verbosity >= 2) { std::cout << "primal cost = " << cost << "\n"; }

      return cost;
   }

   void UpdateFactor(FactorTypeAdapter* f, const weight_vector& omega) // perform one block coordinate step for factor f
   {
      f->UpdateFactor(omega);
   }
   void ComputePass(const INDEX iteration);

   void ComputeForwardPass()
   {
      const auto omega = get_omega();
      ComputePass(forwardUpdateOrdering_.begin(), forwardUpdateOrdering_.end(), omega.forward.begin()); 
   }
   void ComputeBackwardPass()
   {
      const auto omega = get_omega();
      ComputePass(backwardUpdateOrdering_.begin(), backwardUpdateOrdering_.end(), omega.backward.begin());
   }


   void ComputeForwardPassAndPrimal(const INDEX iteration)
   {
      const auto omega = get_omega();
      ComputePassAndPrimal(forwardUpdateOrdering_.begin(), forwardUpdateOrdering_.end(), omega.forward.begin(), 2*iteration+1); // timestamp must be > 0, otherwise in the first iteration primal does not get initialized
      //const REAL forward_cost = EvaluatePrimal();
      //std::cout << "forward cost = " << forward_cost << "\n";
   }
   void ComputeBackwardPassAndPrimal(const INDEX iteration)
   {
      const auto omega = get_omega();
      ComputePassAndPrimal(backwardUpdateOrdering_.begin(), backwardUpdateOrdering_.end(), omega.backward.begin(), 2*iteration + 2); 
      //const REAL backward_cost = EvaluatePrimal();
      //std::cout << "backward cost = " << backward_cost << "\n";
   }
   void ComputePassAndPrimal(const INDEX iteration)
   {
      ComputeForwardPassAndPrimal(iteration);
      ComputeBackwardPassAndPrimal(iteration);
   }

   template<typename FACTOR_ITERATOR, typename OMEGA_ITERATOR>
   void ComputePassAndPrimal(FACTOR_ITERATOR factorIt, const FACTOR_ITERATOR factorEndIt, OMEGA_ITERATOR omegaIt, const INDEX iteration);

   template<typename FACTOR_ITERATOR, typename OMEGA_ITERATOR>
   void ComputePass(FACTOR_ITERATOR factorIt, const FACTOR_ITERATOR factorItEnd, OMEGA_ITERATOR omegaIt);
   void UpdateFactorPrimal(FactorTypeAdapter* f, const weight_vector& omega, const INDEX iteration)
   {
      f->UpdateFactorPrimal(omega, iteration);
   }

   const PrimalSolutionStorage& GetBestPrimal() const;
   template<typename FACTOR_ITERATOR>
      void WritePrimal(FACTOR_ITERATOR factorIt, const FACTOR_ITERATOR factorEndIt, std::ofstream& fs) const;
      void WritePrimal(const INDEX factorIndexBegin, const INDEX factorIndexEnd, std::ofstream& fs) const;

   LPReparametrizationMode GetRepamMode() const { return repamMode_; }

   void set_flags_dirty()
   {
      ordering_valid_ = false;
      omega_anisotropic_valid_ = false;
      omega_anisotropic2_valid_ = false;
      omega_isotropic_valid_ = false;
      omega_isotropic_damped_valid_ = false;
      omega_mixed_valid_ = false;
   }

   // return type for get_omega
   struct omega_storage {
      two_dim_variable_array<REAL>& forward;
      two_dim_variable_array<REAL>& backward;
   };

   omega_storage get_omega()
   {
      assert(repamMode_ != LPReparametrizationMode::Undefined);
      SortFactors();
      if(repamMode_ == LPReparametrizationMode::Anisotropic) {
         ComputeAnisotropicWeights();
         return omega_storage{omegaForwardAnisotropic_, omegaBackwardAnisotropic_};
      } else if(repamMode_ == LPReparametrizationMode::Anisotropic2) {
         ComputeAnisotropicWeights2();
         return omega_storage{omegaForwardAnisotropic2_, omegaBackwardAnisotropic2_};
      } else if(repamMode_ == LPReparametrizationMode::Uniform) {
         ComputeUniformWeights();
         return omega_storage{omegaForwardIsotropic_, omegaBackwardIsotropic_};
      } else if(repamMode_ == LPReparametrizationMode::DampedUniform) {
         ComputeDampedUniformWeights();
         return omega_storage{omegaForwardIsotropicDamped_, omegaBackwardIsotropicDamped_};
      } else if(repamMode_ == LPReparametrizationMode::Mixed) {
         ComputeMixedWeights();
         return omega_storage{omegaForwardMixed_, omegaBackwardMixed_};
      } else {
         throw std::runtime_error("no reparametrization mode set");
      }
   }
protected:
   // do zrobienia: possibly hold factors and messages in shared_ptr?
   std::vector<FactorTypeAdapter*> f_; // note that here the factors are stored in the original order they were given. They will be output in this order as well, e.g. by problemDecomposition
   std::vector<MessageTypeAdapter*> m_;

   bool ordering_valid_ = false;
   std::vector<FactorTypeAdapter*> forwardOrdering_, backwardOrdering_; // separate forward and backward ordering are not needed: Just store factorOrdering_ and generate forward order by begin() and backward order by rbegin().
   std::vector<FactorTypeAdapter*> forwardUpdateOrdering_, backwardUpdateOrdering_; // like forwardOrdering_, but includes only those factors where UpdateFactor actually does something

   bool omega_anisotropic_valid_ = false;
   two_dim_variable_array<REAL> omegaForwardAnisotropic_, omegaBackwardAnisotropic_;
   bool omega_anisotropic2_valid_ = false;
   two_dim_variable_array<REAL> omegaForwardAnisotropic2_, omegaBackwardAnisotropic2_;
   bool omega_isotropic_valid_ = false;
   two_dim_variable_array<REAL> omegaForwardIsotropic_, omegaBackwardIsotropic_;
   bool omega_isotropic_damped_valid_ = false;
   two_dim_variable_array<REAL> omegaForwardIsotropicDamped_, omegaBackwardIsotropicDamped_;
   bool omega_mixed_valid_ = false;
   two_dim_variable_array<REAL> omegaForwardMixed_, omegaBackwardMixed_;

   std::vector<std::pair<FactorTypeAdapter*, FactorTypeAdapter*> > forward_pass_factor_rel_, backward_pass_factor_rel_; // factor ordering relations. First factor must come before second factor. factorRel_ must describe a DAG

   
   std::unordered_map<FactorTypeAdapter*,INDEX> factor_address_to_index_;
   std::vector<INDEX> f_forward_sorted_, f_backward_sorted_; // sorted indices in factor vector f_ 

   LPReparametrizationMode repamMode_ = LPReparametrizationMode::Undefined;

#ifdef LP_MP_PARALLEL
   TCLAP::ValueArg<INDEX> num_lp_threads_arg_;
#endif
};

template<typename BASE_LP_CLASS>
class LP_sat : public BASE_LP_CLASS
{
private:
  struct sat_th {
    sat_th() 
      : th(5.0)
    {}

    void adjust_th(const bool feasible)
    {
      if(feasible) {
        th *= 0.8;
      } else {
        th *= 2.0;
      } 
    }

    REAL th;
  };


public:
  using LP_type = LP_sat<BASE_LP_CLASS>;

  LP_sat(TCLAP::CmdLine& cmd) : BASE_LP_CLASS(cmd)
  {
     sat_ = lglinit();
     assert(sat_ != nullptr);
     //sat_.set_no_simplify(); // seems to make solver much faster
  }

  ~LP_sat()
  {
     // if sat solver is still running, wait until it ends
     //assert(!sat_computation_running());
     if(sat_computation_running()) {
        sat_handle_.wait(); 
     }
     lglrelease(sat_);
  }

  LP_sat(LP_sat& o)
    : 
    sat_var_(lglclone(o.sat_var_)),
    forward_sat_th_(o.forward_sat_th_), 
    backward_sat_th_(backward_sat_th_),
    cur_sat_reduction_direction_(o.cur_sat_reduction_direction_)
  {
     assert(sat_handle_.valid()); //should not be copied and we assume that currently no sat solver is running
  }

  void collect_sat_result()
  {
     if(verbosity >= 2) { 
       std::cout << "collect sat result with threshold = ";
       if(cur_sat_reduction_direction_ == Direction::forward) { 
         std::cout << forward_sat_th_.th;
       } else {
         std::cout << backward_sat_th_.th;
       } 
       std::cout << "\n"; // = " << th.th << "\n";
     }
     const bool feasible = sat_handle_.get();
     if(feasible && !sat_dirty_) {

        for(sat_var i=0; i<lglmaxvar(sat_); ++i) {
           //for(sat_var i=0; i<sat_.nVars(); ++i) {
           //assert(sat_.get_model()[i] == CMSat::l_True || sat_.get_model()[1] == CMSat::l_False);
           //}
        }
        // convert sat solution to original solution format and compute primal cost
        for(INDEX i=0; i<this->f_.size(); ++i) {
           assert(this->factor_address_to_index_[this->f_[i]] == i);
           this->f_[i]->convert_primal(sat_, sat_var_[i]);
        }
        // to do: remove this
        REAL primal_cost = this->EvaluatePrimal();
        if(verbosity >= 2) { std::cout << "sat solution cost = " << primal_cost << "\n"; }
     } else {
       if(verbosity >= 2) { std::cout << "sat not feasible with current threshold\n"; }
     }
  }


  virtual INDEX AddFactor(FactorTypeAdapter* f) 
  {
     // we must wait for the sat solver to finish
     if(sat_computation_running()) {
        sat_handle_.wait();
        collect_sat_result();
     }

     sat_var_.push_back( lglmaxvar(sat_) );
     //sat_var_.push_back(sat_.nVars());
     //std::cout << "number of variables in sat_ = " << sat_var_[sat_var_.size()-1] << "\n";
     f->construct_sat_clauses(sat_);
     INDEX n = BASE_LP_CLASS::AddFactor(f);

     sat_dirty_ = true;

     return n;
  }

   virtual INDEX AddMessage(MessageTypeAdapter* m)
   {
      // we must wait for the sat solver to finish
      if(sat_computation_running()) {
         sat_handle_.wait();
         collect_sat_result();
      }

      assert(this->factor_address_to_index_.find(m->GetLeftFactor()) != this->factor_address_to_index_.end());
      assert(this->factor_address_to_index_.find(m->GetRightFactor()) != this->factor_address_to_index_.end());
      const INDEX left_factor_number = this->factor_address_to_index_[m->GetLeftFactor()];
      const INDEX right_factor_number = this->factor_address_to_index_[m->GetRightFactor()];
      m->construct_sat_clauses(sat_, sat_var_[left_factor_number], sat_var_[right_factor_number]);

      sat_dirty_ = true;

      return BASE_LP_CLASS::AddMessage(m);
   }

   // do zrobienia: only pass sat solver, not while LP object
   static bool solve_sat_problem(LP_type* c, sat_vec<sat_literal> assumptions, sat_th* th)
   {
     for(INDEX i=0; i<assumptions.size(); ++i) {
       lglassume( c->sat_, assumptions[i] );
     }
     for(INDEX i=0; i<lglmaxvar(c->sat_); ++i) {
       lglfreeze(c->sat_, to_literal(i));
     }
     const int sat_ret = lglsat(c->sat_);
     if(verbosity >= 2) { std::cout << "solved sat " << sat_ret << "\n"; }

     const bool feasible = (sat_ret == LGL_SATISFIABLE);
     //solve(&assumptions) == CMSat::l_True;
     //const bool feasible = c->sat_.solve(&assumptions) == CMSat::l_True;
     th->adjust_th(feasible);
     return feasible;
   }

   void ComputeForwardPassAndPrimal(const INDEX iteration)
   {
      const auto omega = this->get_omega();
      if(cur_sat_reduction_direction_ == Direction::forward && !sat_computation_running()) {
         compute_pass_reduce_sat(this->forwardUpdateOrdering_.begin(), this->forwardUpdateOrdering_.end(), omega.forward.begin(), forward_sat_th_);
         cur_sat_reduction_direction_ = Direction::backward; 
      } else {
         this->ComputePass(this->forwardUpdateOrdering_.begin(), this->forwardUpdateOrdering_.end(), omega.forward.begin());
      }
   }
   void ComputeBackwardPassAndPrimal(const INDEX iteration)
   {
      const auto omega = this->get_omega();
      if(cur_sat_reduction_direction_ == Direction::backward && !sat_computation_running()) {
         compute_pass_reduce_sat(this->backwardUpdateOrdering_.begin(), this->backwardUpdateOrdering_.end(), omega.backward.begin(), backward_sat_th_);
         cur_sat_reduction_direction_ = Direction::forward; 
      } else {
         for(auto it = this->backwardUpdateOrdering_.begin(); it != this->backwardUpdateOrdering_.end(); ++it) {
            assert((*it)->no_send_messages() == omega.backward[ std::distance(this->backwardUpdateOrdering_.begin(), it) ].size());
         }
         this->ComputePass(this->backwardUpdateOrdering_.begin(), this->backwardUpdateOrdering_.end(), omega.backward.begin());
      }
   }
   void ComputePassAndPrimal(const INDEX iteration)
   {
      assert(false); // is currently not used
      ComputeForwardPassAndPrimal(2*iteration+1);
      ComputeBackwardPassAndPrimal(2*iteration+2);
   }

   bool sat_computation_running() const
   {
     if(!sat_handle_.valid()) {
        return false; // sat not yet begun.
     }
     const auto sat_state = sat_handle_.wait_for(std::chrono::seconds(0));
     assert(sat_state != std::future_status::deferred); // this should not happen as we launch primal computation immediately.
     return sat_state != std::future_status::ready; 
   }

   template<typename FACTOR_ITERATOR, typename WEIGHT_ITERATOR>
   void compute_pass_reduce_sat(FACTOR_ITERATOR factor_begin, FACTOR_ITERATOR factor_end, WEIGHT_ITERATOR omega_begin, sat_th& th)
   {
      assert(!sat_computation_running());

      if(sat_handle_.valid()) { // do not collect sat result in first round
        collect_sat_result();
      }

      sat_vec<sat_literal> assumptions;
      for(auto it=factor_begin; it!=factor_end; ++it, ++omega_begin) {
         const INDEX factor_number = this->factor_address_to_index_[*it];
         (*it)->UpdateFactorSAT(*omega_begin, th.th, sat_var_[factor_number], assumptions);
      }

      // run sat solver on reduced problem asynchronuously
      if(!sat_handle_.valid()) { 
         if(verbosity >= 2) { std::cout << "start sat calculation\n"; }
         sat_handle_ = std::async(std::launch::async, solve_sat_problem, this, assumptions, &th);
         sat_dirty_ = false;
      } else { 
        if(verbosity >= 2) { std::cout << "restart sat calculation\n"; }
        sat_handle_ = std::async(std::launch::async, solve_sat_problem, this, assumptions, &th);
        sat_dirty_ = false;
      }
   }
private:

   std::vector<sat_var> sat_var_;
   //Glucose::SimpSolver sat_;
   //CMSat::SATSolver sat_;
   LGL* sat_;

   sat_th forward_sat_th_, backward_sat_th_;

   decltype(std::async(std::launch::async, solve_sat_problem, nullptr, sat_vec<sat_literal>{}, nullptr)) sat_handle_;
   Direction cur_sat_reduction_direction_ = Direction::forward;
   // possibly not needed anymore
   bool sat_dirty_ = false; // sat solver is run asynchronously. When factor graph changes, then sat solution cannot be read in anymore and has to be discarded. This flag signifies this case
};


inline void LP::Begin()
{
   CalculatePrimalOffsets();

   repamMode_ = LPReparametrizationMode::Undefined;
   assert(f_.size() > 1); // otherwise we need not perform optimization: Just MaximizePotential f_[0]

#ifdef LP_MP_PARALLEL
   omp_set_num_threads(num_lp_threads_arg_.getValue());
   if(verbosity >= 2) { std::cout << "number of threads = " << num_lp_threads_arg_.getValue() << "\n"; }
#endif
}

inline void LP::ComputePass(const INDEX iteration)
{
   const auto omega = get_omega();
   assert(forwardUpdateOrdering_.size() == omega.forward.size());
   assert(forwardUpdateOrdering_.size() == omega.backward.size());
   ComputeForwardPass();
   ComputeBackwardPass();
}

template<typename FACTOR_ITERATOR, typename OMEGA_ITERATOR>
void LP::ComputePass(FACTOR_ITERATOR factorIt, const FACTOR_ITERATOR factorItEnd, OMEGA_ITERATOR omegaIt)
{
   //assert(std::distance(factorItEnd, factorIt) == std::distance(omegaIt, omegaItEnd));
#pragma omp parallel for 
   for(INDEX i=0; i<std::distance(factorIt, factorItEnd); ++i) {
      UpdateFactor(*(factorIt + i), *(omegaIt + i));
   }
   //for(; factorIt!=factorItEnd; ++factorIt, ++omegaIt) {
   //   UpdateFactor(*factorIt, *omegaIt);
   //}
}

inline void LP::ComputeAnisotropicWeights()
{
   if(!omega_anisotropic_valid_) {
      omega_anisotropic_valid_ = true;
#pragma omp sections
      {
#pragma omp section
         ComputeAnisotropicWeights(forwardOrdering_.begin(), forwardOrdering_.end(), f_forward_sorted_.begin(), f_forward_sorted_.end(), omegaForwardAnisotropic_);
#pragma omp section
         ComputeAnisotropicWeights(backwardOrdering_.begin(), backwardOrdering_.end(), f_backward_sorted_.begin(), f_backward_sorted_.end(), omegaBackwardAnisotropic_);
      }
   }
}

inline void LP::ComputeAnisotropicWeights2()
{
   if(!omega_anisotropic2_valid_) {
      omega_anisotropic2_valid_ = true;
#pragma omp sections
      {
#pragma omp section
         ComputeAnisotropicWeights2(forwardOrdering_.begin(), forwardOrdering_.end(), f_forward_sorted_.begin(), f_forward_sorted_.end(), omegaForwardAnisotropic2_);
#pragma omp section
         ComputeAnisotropicWeights2(backwardOrdering_.begin(), backwardOrdering_.end(), f_backward_sorted_.begin(), f_backward_sorted_.end(), omegaBackwardAnisotropic2_);
      }
   }
}

inline void LP::ComputeUniformWeights()
{
   if(!omega_isotropic_valid_) {
      omega_isotropic_valid_ = true;
#pragma omp sections
      {
#pragma omp section
         ComputeUniformWeights(forwardOrdering_.begin(), forwardOrdering_.end(), omegaForwardIsotropic_, 0.0);
#pragma omp section
         ComputeUniformWeights(backwardOrdering_.begin(), backwardOrdering_.end(), omegaBackwardIsotropic_, 0.0);
      }

      assert(this->backwardUpdateOrdering_.size() == omegaBackwardIsotropic_.size());
      for(auto it = this->backwardUpdateOrdering_.begin(); it != this->backwardUpdateOrdering_.end(); ++it) {
         assert((*it)->no_send_messages() == omegaBackwardIsotropic_[ std::distance(this->backwardUpdateOrdering_.begin(), it) ].size());
      } 

      for(auto it = this->forwardUpdateOrdering_.begin(); it != this->forwardUpdateOrdering_.end(); ++it) {
         assert((*it)->no_send_messages() == omegaForwardIsotropic_[ std::distance(this->forwardUpdateOrdering_.begin(), it) ].size());
      } 
   }
}

inline void LP::ComputeDampedUniformWeights()
{
   if(!omega_isotropic_damped_valid_) {
      omega_isotropic_damped_valid_ = true;
#pragma omp sections
      {
#pragma omp section
         ComputeUniformWeights(forwardOrdering_.begin(), forwardOrdering_.end(), omegaForwardIsotropicDamped_, 1.0);
#pragma omp section
         ComputeUniformWeights(backwardOrdering_.begin(), backwardOrdering_.end(), omegaBackwardIsotropicDamped_, 1.0);
      }
   }
}

// Here we check whether messages constraints are satisfied
bool LP::CheckPrimalConsistency() const
{
   volatile bool consistent=true; // or use std::atomic<bool>?

#pragma omp parallel for shared(consistent)
   for(INDEX i=0; i<m_.size(); ++i) {
      if(!consistent) continue;
      if(!m_[i]->CheckPrimalConsistency()) {
         consistent = false;
      }
   }
   if(verbosity >= 2) { std::cout << "primal solution consistent: " << consistent << "\n"; }
   return consistent;
}

// write primal solutions in bounds [factorIndexBegin,factorIndexEnd) to filestream
template<typename FACTOR_ITERATOR>
void LP::WritePrimal(FACTOR_ITERATOR factorIt, const FACTOR_ITERATOR factorEndIt, std::ofstream& fs) const
{
  assert(false); // write primal is implemented in problem constructor
   for(; factorIt!=factorEndIt; ++factorIt) {
      (*factorIt)->WritePrimal(fs);
   }
}

// only compute factors adjacent to which also messages can be send
template<typename FACTOR_ITERATOR>
std::vector<std::vector<FactorTypeAdapter*> > LP::ComputeSendFactorConnection(FACTOR_ITERATOR factorIt, FACTOR_ITERATOR factorEndIt)
{
   std::vector<std::vector<FactorTypeAdapter*> > fc;
   fc.reserve(factorEndIt - factorIt);
   for(;factorIt != factorEndIt; ++factorIt) {
      fc.push_back({});
      for(auto mIt=(*factorIt)->begin(); mIt!=(*factorIt)->end(); ++mIt) {
         if(mIt.CanSendMessage()) {
            fc.back().push_back( mIt.GetConnectedFactor() );
         }
      }
   }

   return fc;
}

template<typename FACTOR_ITERATOR>
std::vector<std::vector<FactorTypeAdapter*> > LP::ComputeFactorConnection(FACTOR_ITERATOR factorIt, FACTOR_ITERATOR factorEndIt)
{
   std::vector<std::vector<FactorTypeAdapter*> > fc;
   fc.reserve(factorEndIt - factorIt);
   for(;factorIt != factorEndIt; ++factorIt) {
      fc.push_back({});
      for(auto mIt=(*factorIt)->begin(); mIt!=(*factorIt)->end(); ++mIt) {
         fc.back().push_back( mIt.GetConnectedFactor() );
      }
   }

   // note: this need not hold for e.g. global assignment factor
   //for(INDEX i=0; i<fc.size(); ++i) {
   //   assert(HasUniqueValues(fc[i]));
   //}

   return fc;
}


template<typename FACTOR_ITERATOR, typename FACTOR_SORT_ITERATOR>
void LP::ComputeAnisotropicWeights2(
      FACTOR_ITERATOR factorIt, FACTOR_ITERATOR factorEndIt, // sorted pointers to factors
      FACTOR_SORT_ITERATOR factor_sort_begin, FACTOR_SORT_ITERATOR factor_sort_end, // sorted factor indices in f_
      two_dim_variable_array<REAL>& omega)
{
   std::vector<INDEX> f_sorted_inverse(std::distance(factor_sort_begin, factor_sort_end)); // factor index in order they were added to sorted order
#pragma omp parallel for
   for(INDEX i=0; i<std::distance(factor_sort_begin, factor_sort_end); ++i) {
      f_sorted_inverse[ factor_sort_begin[i] ] = i;
   }

   assert(std::distance(factorIt,factorEndIt) == f_.size());
   assert(std::distance(factor_sort_begin, factor_sort_end) == f_.size());

   std::vector<INDEX> no_send_messages_later(f_.size(), 0);
   for(INDEX i=0; i<m_.size(); ++i) {
      auto* f_left = m_[i]->GetLeftFactor();
      const INDEX f_index_left = factor_address_to_index_[f_left];
      const INDEX index_left = f_sorted_inverse[f_index_left];
      auto* f_right = m_[i]->GetRightFactor();
      const INDEX f_index_right = factor_address_to_index_[f_right];
      const INDEX index_right = f_sorted_inverse[f_index_right];
      
      if(m_[i]->SendsMessageToRight() && index_left < index_right) {
        no_send_messages_later[index_left]++;
      }
      if(m_[i]->SendsMessageToLeft() && index_right < index_left) {
        no_send_messages_later[index_right]++;
      }
   }

   std::vector<INDEX> omega_size(f_.size());
   {
     INDEX c=0;
     for(auto it=factorIt; it!=factorEndIt; ++it) {
       if((*it)->FactorUpdated()) {
         omega_size[c] = (*it)->no_send_messages();
         ++c;
       }
     }
     omega_size.resize(c);
   }

   omega = two_dim_variable_array<REAL>(omega_size);

   {
      INDEX c=0;
      for(auto it=factorIt; it!=factorEndIt; ++it) {
         const INDEX i = std::distance(factorIt, it);
         assert(i == f_sorted_inverse[ factor_address_to_index_[*it] ]);
         if((*it)->FactorUpdated()) {
            INDEX k=0;
            for(auto mIt=(*it)->begin(); mIt!=(*it)->end(); ++mIt) {
               if(mIt.CanSendMessage()) {
                  auto* f_connected = mIt.GetConnectedFactor();
                  const INDEX j = f_sorted_inverse[ factor_address_to_index_[f_connected] ];
                  assert(i != j);
                  if(i<j) {
                     omega[c][k] = 1.0/REAL(no_send_messages_later[i]);
                  } else {
                     omega[c][k] = 0.0;
                  } 
                  ++k;
               }
            }
            ++c;
         }
      }
   }
}

// do zrobienia: possibly templatize this for use with iterators
// note: this function is not working properly. We should only compute factors for messages which can actually send
template<typename FACTOR_ITERATOR, typename FACTOR_SORT_ITERATOR>
void LP::ComputeAnisotropicWeights(
      FACTOR_ITERATOR factorIt, FACTOR_ITERATOR factorEndIt, // sorted pointers to factors
      FACTOR_SORT_ITERATOR factor_sort_begin, FACTOR_SORT_ITERATOR factor_sort_end, // sorted factor indices in f_
      two_dim_variable_array<REAL>& omega)
{
   std::vector<INDEX> f_sorted_inverse(std::distance(factor_sort_begin, factor_sort_end)); // factor index in order they were added to sorted order
#pragma omp parallel for
   for(INDEX i=0; i<std::distance(factor_sort_begin, factor_sort_end); ++i) {
      f_sorted_inverse[ factor_sort_begin[i] ] = i;
   }

   assert(std::distance(factorIt,factorEndIt) == f_.size());
   assert(std::distance(factor_sort_begin, factor_sort_end) == f_.size());

   // compute the following numbers: 
   // 1) #{factors after current one, to which messages are sent from current factor}
   // 2) #{factors after current one, which receive messages from current one}
#ifdef LP_MP_PARALLEL
   std::atomic<INDEX>* no_send_factors = new std::atomic<INDEX>[f_.size()];
   std::fill(no_send_factors, no_send_factors + f_.size(), 0);
   std::atomic<INDEX>* no_send_factors_later = new std::atomic<INDEX>[f_.size()];
   std::fill(no_send_factors_later, no_send_factors_later + f_.size(), 0);
   std::atomic<INDEX>* no_receiving_factors_later = new std::atomic<INDEX>[f_.size()];
   std::fill(no_receiving_factors_later, no_receiving_factors_later + f_.size(), 0);
   std::atomic<INDEX>* last_receiving_factor = new std::atomic<INDEX>[f_.size()];
   std::fill(last_receiving_factor, last_receiving_factor + f_.size(), 0);
#else
   std::vector<INDEX> no_send_factors(f_.size(),0);
   std::vector<INDEX> no_send_factors_later(f_.size(),0);
   std::vector<INDEX> no_receiving_factors_later(f_.size(),0); // numer of factors later than current one receiving message from current one
   std::vector<INDEX> last_receiving_factor(f_.size(), 0); // what is the last (in the order given by factor iterator) factor that receives a message?
#endif

   // do zrobienia: if factor is not visited at all, then omega is not needed for that entry. We must filter out such entries still
#pragma omp parallel for
   for(INDEX i=0; i<m_.size(); ++i) {
      auto* f_left = m_[i]->GetLeftFactor();
      const INDEX f_index_left = factor_address_to_index_[f_left];
      const INDEX index_left = f_sorted_inverse[f_index_left];
      auto* f_right = m_[i]->GetRightFactor();
      const INDEX f_index_right = factor_address_to_index_[f_right];
      const INDEX index_right = f_sorted_inverse[f_index_right];
      
      if(m_[i]->ReceivesMessageFromLeft()) {
         if(index_left < index_right) {
            no_receiving_factors_later[index_left]++;
         }
#ifdef LP_MP_PARALLEL
         INDEX old_val = last_receiving_factor[index_left];
         const INDEX new_val = std::max(old_val, index_right);
         while(old_val < new_val && !last_receiving_factor[index_left].compare_exchange_weak(old_val, new_val)) ;
#else
         last_receiving_factor[index_left] = std::max(last_receiving_factor[index_left], index_right);
#endif
      }

      if(m_[i]->ReceivesMessageFromRight()) {
         if(index_left > index_right) {
            no_receiving_factors_later[index_right]++;
         }
#ifdef LP_MP_PARALLEL
         INDEX old_val = last_receiving_factor[index_right];
         const INDEX new_val = std::max(old_val, index_left);
         while(old_val < new_val && !last_receiving_factor[index_right].compare_exchange_weak(old_val, new_val)) ;
#else
         last_receiving_factor[index_right] = std::max(last_receiving_factor[index_right], index_left);
#endif
      }
   }

#pragma omp parallel for
   for(INDEX i=0; i<m_.size(); ++i) {
      auto* f_left = m_[i]->GetLeftFactor();
      const INDEX f_index_left = factor_address_to_index_[f_left];
      const INDEX index_left = f_sorted_inverse[f_index_left];
      auto* f_right = m_[i]->GetRightFactor();
      const INDEX f_index_right = factor_address_to_index_[f_right];
      const INDEX index_right = f_sorted_inverse[f_index_right];

      if(m_[i]->SendsMessageToRight()) {
         no_send_factors[index_left]++;
         if(index_left < index_right || last_receiving_factor[index_right] > index_left) {
            no_send_factors_later[index_left]++;
         }
      }
      if(m_[i]->SendsMessageToLeft()) {
         no_send_factors[index_right]++;
         if(index_right < index_left || last_receiving_factor[index_left] > index_right) {
            no_send_factors_later[index_right]++;
         }
      }
   }

   std::vector<INDEX> omega_size(f_.size());
   {
     INDEX c=0;
     for(auto it=factorIt; it!=factorEndIt; ++it) {
       if((*it)->FactorUpdated()) {
         omega_size[c] = (*it)->no_send_messages();
         ++c;
       }
     }
     omega_size.resize(c);
   }

   omega = two_dim_variable_array<REAL>(omega_size);

   {
      INDEX c=0;
      for(auto it=factorIt; it!=factorEndIt; ++it) {
         const INDEX i = std::distance(factorIt, it);
         assert(i == f_sorted_inverse[ factor_address_to_index_[*it] ]);
         if((*it)->FactorUpdated()) {
            INDEX k=0;
            for(auto mIt=(*it)->begin(); mIt!=(*it)->end(); ++mIt) {
               if(mIt.CanSendMessage()) {
                  auto* f_connected = mIt.GetConnectedFactor();
                  const INDEX j = f_sorted_inverse[ factor_address_to_index_[f_connected] ];
                  assert(i != j);
                  if(i<j || last_receiving_factor[j] > i) {
                     omega[c][k] = (1.0/REAL(no_receiving_factors_later[i] + std::max(INDEX(no_send_factors_later[i]), INDEX(no_send_factors[i]) - INDEX(no_send_factors_later[i]))));
                  } else {
                     omega[c][k] = 0.0;
                  } 
                  ++k;
               }
            }
            ++c;
         }
      }
   }


   // check whether all messages were added to m_. Possibly, this can be automated: Traverse all factors, get all messages, add them to m_ and avoid duplicates along the way.
   assert(2*m_.size() == std::accumulate(f_.begin(), f_.end(), 0, [](INDEX sum, auto* f){ return sum + f->GetNoMessages(); }));
   assert(HasUniqueValues(m_));
   for(INDEX i=0; i<omega.size(); ++i) {
      //assert(omega[i].size() <= (*(factorIt+i))->GetNoMessages());
      //const REAL omega_sum = std::accumulate(omega[i].begin(), omega[i].end(), 0.0); 
      assert(std::accumulate(omega[i].begin(), omega[i].end(), 0.0) <= 1.0 + eps);
   }
#ifdef LP_MP_PARALLEL
   delete[] no_send_factors;
   delete[] no_send_factors_later;
   delete[] no_receiving_factors_later;
   delete[] last_receiving_factor; 
#endif
}

// compute uniform weights so as to help decoding for obtaining primal solutions
// leave_weight signals how much weight to leave in sending factor. Important for rounding and tightening
// to do: possibly pass update factor list (i.e. only those factors which are updated)
template<typename FACTOR_ITERATOR>
void LP::ComputeUniformWeights(FACTOR_ITERATOR factorIt, FACTOR_ITERATOR factorEndIt, two_dim_variable_array<REAL>& omega, const REAL leave_weight)
{
   assert(leave_weight >= 0.0 && leave_weight <= 1.0);
   assert(factorEndIt - factorIt == f_.size());

   std::vector<INDEX> omega_size(std::distance(factorIt, factorEndIt),0);
   INDEX c=0;
   for(; factorIt != factorEndIt; ++factorIt) {
      if((*factorIt)->FactorUpdated()) {
        omega_size[c] = (*factorIt)->no_send_messages();
        ++c;
      }
   }

   omega_size.resize(c);
   omega = two_dim_variable_array<REAL>(omega_size);

   assert(omega.size() == omega_size.size());
   for(INDEX i=0; i<omega.size(); ++i) {
      assert(omega[i].size() == omega_size[i]);
   }

#pragma omp parallel for
   for(INDEX i=0; i<omega.size(); ++i) {
     for(INDEX j=0; j<omega_size[i]; ++j) {
       omega[i][j] = 1.0/REAL( omega_size[i] + leave_weight );
     }
   }
}

// compute anisotropic and damped uniform weights, then average them
void LP::ComputeMixedWeights(
      const two_dim_variable_array<REAL>& omega_anisotropic,
      const two_dim_variable_array<REAL>& omega_damped_uniform,
      two_dim_variable_array<REAL>& omega)
{
   omega = omega_anisotropic;
   
   assert(omega_damped_uniform.size() == omega.size());
#pragma omp parallel for
   for(INDEX i=0; i<omega.size(); ++i) {
      assert(omega_damped_uniform[i].size() == omega[i].size());
      for(INDEX j=0; j<omega[i].size(); ++j) {
         omega[i][j] = 0.5*(omega[i][j] + omega_damped_uniform[i][j]);
      }
   }
}

template<typename FACTOR_ITERATOR, typename OMEGA_ITERATOR>
void LP::ComputePassAndPrimal(FACTOR_ITERATOR factorIt, const FACTOR_ITERATOR factorEndIt, OMEGA_ITERATOR omegaIt, INDEX iteration)
{
   //possibly do not use parallelization here
#pragma omp parallel for 
   for(INDEX i=0; i<std::distance(factorIt, factorEndIt); ++i) {
      UpdateFactorPrimal(*(factorIt + i), *(omegaIt + i), iteration);
   }
}


} // end namespace LP_MP

#endif // LP_MP_MAIN


