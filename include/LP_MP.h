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
#include <unordered_map>
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
#include "serialization.hxx"
#include "tclap/CmdLine.h"
#include "DD_ILP.hxx"

#ifdef LP_MP_PARALLEL
#include <omp.h>
#endif

namespace LP_MP {

// forward declaration
class MessageTypeAdapter;
class MessageIterator;

using weight_array = two_dim_variable_array<REAL>;
using weight_slice = two_dim_variable_array<REAL>::ArrayAccessObject;
using receive_array = two_dim_variable_array<unsigned char>;
using receive_slice = two_dim_variable_array<unsigned char>::ArrayAccessObject;

// pure virtual base class for factor container used by LP class
class FactorTypeAdapter
{
public:
   virtual ~FactorTypeAdapter() {}
   virtual FactorTypeAdapter* clone() const = 0;
   virtual void UpdateFactor(const weight_slice& omega, const receive_slice& receive_mask) = 0;
   virtual void update_factor_residual(const weight_slice& omega) = 0;
   virtual void UpdateFactorPrimal(const weight_slice& omega, const receive_slice& receive_mask, const INDEX iteration) = 0;
#ifdef LP_MP_PARALLEL
   virtual void UpdateFactorSynchronized(const weight_slice& omega) = 0;
   virtual void UpdateFactorPrimalSynchronized(const weight_slice& omega, const INDEX iteration) = 0;
#endif
   virtual bool SendsMessage(const INDEX msg_idx) const = 0;
   virtual bool ReceivesMessage(const INDEX msg_idx) const = 0;
   virtual bool FactorUpdated() const = 0; // does calling UpdateFactor do anything? If no, it need not be called while in ComputePass, saving time.
   // to do: remove both
   MessageIterator begin(); 
   MessageIterator end();
   virtual INDEX no_messages() const = 0;
   virtual INDEX no_send_messages() const = 0;
   virtual INDEX no_receive_messages() const = 0;
   virtual const MessageTypeAdapter* GetMessage(const INDEX n) const = 0;
   virtual FactorTypeAdapter* GetConnectedFactor(const INDEX i) const = 0;
   virtual REAL LowerBound() const = 0;
   virtual void init_primal() = 0;
   virtual void MaximizePotentialAndComputePrimal() = 0;
   virtual void propagate_primal_through_messages() = 0;
   virtual bool check_primal_consistency() = 0;

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

   // do zrobienia: this function is not needed. Evaluation can be performed automatically
   virtual REAL EvaluatePrimal() const = 0;

   // external ILP-interface
   virtual void construct_constraints(DD_ILP::external_solver_interface<DD_ILP::sat_solver>& solver) = 0;
   virtual void load_costs(DD_ILP::external_solver_interface<DD_ILP::sat_solver>& solver) = 0;
   virtual void convert_primal(DD_ILP::external_solver_interface<DD_ILP::sat_solver>& solver) = 0; 

   virtual void construct_constraints(DD_ILP::external_solver_interface<DD_ILP::problem_export>& solver) = 0;
   virtual void load_costs(DD_ILP::external_solver_interface<DD_ILP::problem_export>& solver) = 0; 
   virtual void convert_primal(DD_ILP::external_solver_interface<DD_ILP::problem_export>& solver) = 0;

#ifdef DD_ILP_WITH_GUROBI
   virtual void construct_constraints(DD_ILP::external_solver_interface<DD_ILP::gurobi_interface>& solver) = 0;
   virtual void load_costs(DD_ILP::external_solver_interface<DD_ILP::gurobi_interface>& solver) = 0;
   virtual void convert_primal(DD_ILP::external_solver_interface<DD_ILP::gurobi_interface>& solver) = 0;
#endif

   // estimate of how long a factor update will take
   virtual INDEX runtime_estimate() = 0;
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
   
   // external ILP-interface
   virtual void construct_constraints(DD_ILP::external_solver_interface<DD_ILP::sat_solver>&, const DD_ILP::variable_counters&, const DD_ILP::variable_counters&) = 0;
   virtual void construct_constraints(DD_ILP::external_solver_interface<DD_ILP::problem_export>&, const DD_ILP::variable_counters&, const DD_ILP::variable_counters&) = 0;
#ifdef DD_ILP_WITH_GUROBI
   virtual void construct_constraints(DD_ILP::external_solver_interface<DD_ILP::gurobi_interface>&, const DD_ILP::variable_counters&, const DD_ILP::variable_counters&) = 0;
#endif
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
   const MessageTypeAdapter& operator*() const { return *(factor_->GetMessage(msg_idx_)); }
   const MessageTypeAdapter& operator->() const { return *(factor_->GetMessage(msg_idx_)); }
   FactorTypeAdapter* GetConnectedFactor() const { return factor_->GetConnectedFactor(msg_idx_); }
   bool SendsMessage() const  { return factor_->SendsMessage(msg_idx_); }
   bool ReceivesMessage() const  { return factor_->ReceivesMessage(msg_idx_); }
private:
   FactorTypeAdapter* const factor_;
   INDEX msg_idx_;
};

inline MessageIterator FactorTypeAdapter::begin() { return MessageIterator(this,0); }
inline MessageIterator FactorTypeAdapter::end()  { return MessageIterator(this, no_messages()); }


class LP {
   struct message_trait
   {
       FactorTypeAdapter* left;
       FactorTypeAdapter* right;
       const bool sends_message_to_left, sends_message_to_right, receives_message_from_left, receives_message_from_right;
   };

public:
   LP(TCLAP::CmdLine& cmd);
   ~LP();
   LP(LP& o);
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

   virtual INDEX AddFactor(FactorTypeAdapter* f);
   INDEX GetNumberOfFactors() const { return f_.size(); }
   FactorTypeAdapter* GetFactor(const INDEX i) const { return f_[i]; }

   template<typename MESSAGE_CONTAINER_TYPE, typename LEFT_FACTOR, typename RIGHT_FACTOR, typename... ARGS>
   INDEX add_message(LEFT_FACTOR* l, RIGHT_FACTOR* r, ARGS... args)
   {
       set_flags_dirty();

       auto* m_l = l->template add_message<MESSAGE_CONTAINER_TYPE,Chirality::left>(r,args...);
       auto* m_r = r->template add_message<MESSAGE_CONTAINER_TYPE,Chirality::right>(l,args...);

       l->set_left_msg(m_r);
       r->set_right_msg(m_l); 

       auto* m = (m_l != nullptr ? m_l : m_r);
       assert(m != nullptr);
       m_.push_back({l,r, m->SendsMessageToLeft(), m->SendsMessageToRight(), m->ReceivesMessageFromLeft(), m->ReceivesMessageFromRight()});

       return m_.size() + 1;
   }
   //virtual INDEX AddMessage(MessageTypeAdapter* m);
   message_trait GetMessage(const INDEX i) const { return m_[i]; }
   INDEX GetNumberOfMessages() const { return m_.size(); }

   void AddFactorRelation(FactorTypeAdapter* f1, FactorTypeAdapter* f2); // indicate that factor f1 comes before factor f2
   void ForwardPassFactorRelation(FactorTypeAdapter* f1, FactorTypeAdapter* f2);
   void BackwardPassFactorRelation(FactorTypeAdapter* f1, FactorTypeAdapter* f2);

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
   void End() {};

   void SortFactors(
         const std::vector<std::pair<FactorTypeAdapter*, FactorTypeAdapter*>>& factor_rel,
         std::vector<FactorTypeAdapter*>& ordering,
         std::vector<FactorTypeAdapter*>& update_ordering,
         std::vector<INDEX>& f_sorted
         );

   void SortFactors();

   //void ComputeWeights(const LPReparametrizationMode m);
   void set_reparametrization(const LPReparametrizationMode r) { repamMode_ = r; }

   bool omega_valid(const weight_array& omega) const;

   void ComputeAnisotropicWeights();
   template<typename FACTOR_ITERATOR, typename FACTOR_SORT_ITERATOR, typename FACTOR_MASK_ITERATOR>
   void ComputeAnisotropicWeights(FACTOR_ITERATOR factorIt, FACTOR_ITERATOR factorItEnd, FACTOR_SORT_ITERATOR factor_sort_begin, FACTOR_SORT_ITERATOR factor_sort_end, FACTOR_MASK_ITERATOR factor_mask_begin, FACTOR_MASK_ITERATOR factor_mask_end, weight_array& omega, receive_array& receive_mask); 

   void ComputeAnisotropicWeights2();
   template<typename FACTOR_ITERATOR, typename FACTOR_SORT_ITERATOR>
   void ComputeAnisotropicWeights2(FACTOR_ITERATOR factorIt, FACTOR_ITERATOR factorItEnd, FACTOR_SORT_ITERATOR factor_sort_begin, FACTOR_SORT_ITERATOR factor_sort_end, two_dim_variable_array<REAL>& omega); 

   void ComputeUniformWeights();

   void ComputeDampedUniformWeights();
   template<typename FACTOR_ITERATOR, typename FACTOR_MASK_ITERATOR>
   void ComputeUniformWeights(FACTOR_ITERATOR factorIt, FACTOR_ITERATOR factorEndIt, FACTOR_MASK_ITERATOR factor_mask_begin, FACTOR_MASK_ITERATOR factor_mask_end, two_dim_variable_array<REAL>& omega, const REAL leave_weight); // do zrobienia: rename to isotropic weights

   void ComputeMixedWeights(const two_dim_variable_array<REAL>& omega_anisotropic, const two_dim_variable_array<REAL>& omega_damped_uniform, two_dim_variable_array<REAL>& omega); 
   void ComputeMixedWeights();

   void compute_full_receive_mask();
   template<typename FACTOR_ITERATOR>
   void compute_full_receive_mask(FACTOR_ITERATOR factor_begin, FACTOR_ITERATOR factor_end, receive_array& receive_mask);

   double LowerBound() const;
   double EvaluatePrimal();

   bool CheckPrimalConsistency() const;

   void ComputePass(const INDEX iteration);
   template<typename FACTOR_MASK_ITERATOR>
   void ComputePass(const INDEX iteration, FACTOR_MASK_ITERATOR factor_mask_begin, FACTOR_MASK_ITERATOR factor_mask_end);
   void ComputeForwardPass();
   void ComputeBackwardPass();

   void ComputePassAndPrimal(const INDEX iteration);
   void ComputeForwardPassAndPrimal(const INDEX iteration);
   void ComputeBackwardPassAndPrimal(const INDEX iteration);

   template<typename FACTOR_ITERATOR, typename OMEGA_ITERATOR, typename RECEIVE_MASK_ITERATOR>
   void ComputePassAndPrimal(FACTOR_ITERATOR factorIt, const FACTOR_ITERATOR factorEndIt, OMEGA_ITERATOR omegaIt, RECEIVE_MASK_ITERATOR receive_mask_it, const INDEX iteration);

   template<typename FACTOR_ITERATOR, typename OMEGA_ITERATOR, typename RECEIVE_MASK_ITERATOR>
   void ComputePass(FACTOR_ITERATOR factorIt, const FACTOR_ITERATOR factorItEnd, OMEGA_ITERATOR omegaIt, RECEIVE_MASK_ITERATOR receive_it);

#ifdef LP_MP_PARALLEL
   template<typename FACTOR_ITERATOR, typename OMEGA_ITERATOR, typename SYNCHRONIZATION_ITERATOR>
   void ComputePassSynchronized(
       FACTOR_ITERATOR factorIt, const FACTOR_ITERATOR factorItEnd, 
       OMEGA_ITERATOR omega_begin, OMEGA_ITERATOR omega_end,
       SYNCHRONIZATION_ITERATOR synchronization_begin, SYNCHRONIZATION_ITERATOR synchronization_end);

   template<typename FACTOR_ITERATOR, typename OMEGA_ITERATOR, typename SYNCHRONIZATION_ITERATOR>
   void ComputePassAndPrimalSynchronized(FACTOR_ITERATOR factorIt, const FACTOR_ITERATOR factorEndIt, OMEGA_ITERATOR omegaIt, SYNCHRONIZATION_ITERATOR, const INDEX iteration);
#endif

   //const PrimalSolutionStorage& GetBestPrimal() const;

   LPReparametrizationMode GetRepamMode() const { return repamMode_; }

   void set_flags_dirty();

   // return type for get_omega
   struct omega_storage {
      two_dim_variable_array<REAL>& forward;
      two_dim_variable_array<REAL>& backward;
      receive_array& receive_mask_forward;
      receive_array& receive_mask_backward;
   };

   omega_storage get_omega()
   {
      assert(repamMode_ != LPReparametrizationMode::Undefined);
      SortFactors();

#ifdef LP_MP_PARALLEL
      compute_synchronization();
#endif 
      if(repamMode_ != LPReparametrizationMode::Anisotropic) {
          if(!full_receive_mask_valid_) {
              compute_full_receive_mask();
              full_receive_mask_valid_ = true;
          }
      }

      if(repamMode_ == LPReparametrizationMode::Anisotropic) {
        if(!omega_anisotropic_valid_) {
          ComputeAnisotropicWeights();
          omega_anisotropic_valid_ = true;
        }
        return omega_storage{omegaForwardAnisotropic_, omegaBackwardAnisotropic_, anisotropic_receive_mask_forward_, anisotropic_receive_mask_backward_};
      } else if(repamMode_ == LPReparametrizationMode::Anisotropic2) {
        if(!omega_anisotropic2_valid_) {
          ComputeAnisotropicWeights2();
          omega_anisotropic2_valid_ = true;
        }
        return omega_storage{omegaForwardAnisotropic2_, omegaBackwardAnisotropic2_, full_receive_mask_backward_, full_receive_mask_backward_};
      } else if(repamMode_ == LPReparametrizationMode::Uniform) {
        if(!omega_isotropic_valid_) {
          ComputeUniformWeights();
          omega_isotropic_valid_ = true;
        }
        return omega_storage{omegaForwardIsotropic_, omegaBackwardIsotropic_, full_receive_mask_forward_, full_receive_mask_backward_};
      } else if(repamMode_ == LPReparametrizationMode::DampedUniform) {
        if(!omega_isotropic_damped_valid_) {
          ComputeDampedUniformWeights();
          omega_isotropic_damped_valid_ = true;
        }
        return omega_storage{omegaForwardIsotropicDamped_, omegaBackwardIsotropicDamped_, full_receive_mask_backward_, full_receive_mask_backward_};
      } else if(repamMode_ == LPReparametrizationMode::Mixed) {
        if(!omega_mixed_valid_) {
          ComputeMixedWeights();
          omega_mixed_valid_ = true;
        }
        return omega_storage{omegaForwardMixed_, omegaBackwardMixed_, full_receive_mask_backward_, full_receive_mask_backward_};
      } else {
        throw std::runtime_error("no reparametrization mode set");
      }

   }

   void add_to_constant(const REAL x) { constant_ += x; }
protected:

   // do zrobienia: possibly hold factors and messages in shared_ptr?
   std::vector<FactorTypeAdapter*> f_; // note that here the factors are stored in the original order they were given. They will be output in this order as well, e.g. by problemDecomposition
   std::vector<message_trait> m_;

   bool ordering_valid_ = false;
   std::vector<FactorTypeAdapter*> forwardOrdering_, backwardOrdering_; // separate forward and backward ordering are not needed: Just store factorOrdering_ and generate forward order by begin() and backward order by rbegin().
   std::vector<FactorTypeAdapter*> forwardUpdateOrdering_, backwardUpdateOrdering_; // like forwardOrdering_, but includes only those factors where UpdateFactor actually does something

   bool omega_anisotropic_valid_ = false;
   two_dim_variable_array<REAL> omegaForwardAnisotropic_, omegaBackwardAnisotropic_;
   receive_array anisotropic_receive_mask_forward_, anisotropic_receive_mask_backward_;
   bool omega_anisotropic2_valid_ = false;
   two_dim_variable_array<REAL> omegaForwardAnisotropic2_, omegaBackwardAnisotropic2_;
   bool omega_isotropic_valid_ = false;
   two_dim_variable_array<REAL> omegaForwardIsotropic_, omegaBackwardIsotropic_;
   bool omega_isotropic_damped_valid_ = false;
   two_dim_variable_array<REAL> omegaForwardIsotropicDamped_, omegaBackwardIsotropicDamped_;
   bool omega_mixed_valid_ = false;
   two_dim_variable_array<REAL> omegaForwardMixed_, omegaBackwardMixed_;

   bool full_receive_mask_valid_ = false;
   receive_array full_receive_mask_forward_, full_receive_mask_backward_;

   std::vector<std::pair<FactorTypeAdapter*, FactorTypeAdapter*> > forward_pass_factor_rel_, backward_pass_factor_rel_; // factor ordering relations. First factor must come before second factor. factorRel_ must describe a DAG

   
   std::unordered_map<FactorTypeAdapter*,INDEX> factor_address_to_index_;
   std::vector<INDEX> f_forward_sorted_, f_backward_sorted_; // sorted indices in factor vector f_ 

   LPReparametrizationMode repamMode_ = LPReparametrizationMode::Undefined;

   TCLAP::ValueArg<std::string> reparametrization_type_arg_; // shared|residual
   enum class reparametrization_type {shared,residual};
   reparametrization_type reparametrization_type_;
#ifdef LP_MP_PARALLEL
   TCLAP::ValueArg<INDEX> num_lp_threads_arg_;
   bool synchronization_valid_ = false;
   std::vector<bool> synchronize_forward_;
   std::vector<bool> synchronize_backward_;

   template<typename ITERATOR>
   std::vector<bool> compute_synchronization(ITERATOR factor_begin, ITERATOR factor_end);

   // determine for which factor updates synchronization must be enabled
   void compute_synchronization()
   {
     assert(ordering_valid_);
     if(synchronization_valid_) { return; }
     synchronization_valid_ = true;

     synchronize_forward_ = compute_synchronization(forwardUpdateOrdering_.begin(), forwardUpdateOrdering_.end());
     synchronize_backward_ = compute_synchronization(backwardUpdateOrdering_.begin(), backwardUpdateOrdering_.end()); 
   }
#endif

   std::vector<bool> get_inconsistent_mask(const std::size_t no_fatten_rounds = 1);
   template<typename FACTOR_ITERATOR, typename FACTOR_MASK_ITERATOR>
   std::vector<FactorTypeAdapter*> get_masked_factors( FACTOR_ITERATOR factor_begin, FACTOR_ITERATOR factor_end, FACTOR_MASK_ITERATOR factor_mask_begin, FACTOR_MASK_ITERATOR factor_mask_end);
   void reduce_optimization_factors();

   std::vector<bool> factor_mask;

   REAL constant_ = 0;
};

LP::LP(TCLAP::CmdLine& cmd)
: reparametrization_type_arg_("","reparametrizationType","message sending type: ", false, "shared", "{shared|residual}", cmd)
#ifdef LP_MP_PARALLEL
, num_lp_threads_arg_("","numLpThreads","number of threads for message passing, default = 1",false,1,&positiveIntegerConstraint,cmd)
#endif
{}

LP::~LP()
{
  for(INDEX i=0; i<f_.size(); i++) { delete f_[i]; }
}

// make a deep copy of factors and messages. Adjust pointers to messages and factors
LP::LP(LP& o) // no const because of o.num_lp_threads_arg_.getValue() not being const!
  : reparametrization_type_arg_("","reparametrizationType","message sending type: ", false, o.reparametrization_type_arg_.getValue(), "{shared|residual}" )
#ifdef LP_MP_PARALLEL
    , num_lp_threads_arg_("","numLpThreads","number of threads for message passing, default = 1",false,o.num_lp_threads_arg_.getValue(),&positiveIntegerConstraint)
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
  for(auto m : o.m_) {
    auto* left = m.left;
    auto* right = m.right;
    auto* left_clone = factor_map[left];
    auto* right_clone = factor_map[right];
    m_.push_back({left_clone, right_clone, m.sends_message_to_left, m.sends_message_to_right, m.receives_message_from_left, m.receives_message_from_right});
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

  constant_ = o.constant_;
}

INDEX LP::AddFactor(FactorTypeAdapter* f)
{
  set_flags_dirty();
  assert(factor_address_to_index_.size() == f_.size());
  f_.push_back(f);

  assert(factor_address_to_index_.find(f) == factor_address_to_index_.end());
  factor_address_to_index_.insert(std::make_pair(f,f_.size()-1));
  assert(factor_address_to_index_.find(f)->second == f_.size()-1);

  factor_mask.push_back(true);
  assert(factor_mask.size() == f_.size());

  return f_.size() - 1;
}

//INDEX LP::AddMessage(MessageTypeAdapter* m)
//{
//  set_flags_dirty();
//  m_.push_back(m);
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

//  return m_.size() - 1;
//}

void LP::ForwardPassFactorRelation(FactorTypeAdapter* f1, FactorTypeAdapter* f2) 
{ 
  set_flags_dirty();
  assert(f1!=f2);
  forward_pass_factor_rel_.push_back({f1,f2}); 
}

void LP::BackwardPassFactorRelation(FactorTypeAdapter* f1, FactorTypeAdapter* f2) 
{ 
  set_flags_dirty(); 
  assert(f1!=f2); 
  backward_pass_factor_rel_.push_back({f1,f2}); 
}

void LP::AddFactorRelation(FactorTypeAdapter* f1, FactorTypeAdapter* f2)
{
  ForwardPassFactorRelation(f1,f2);
  BackwardPassFactorRelation(f2,f1);
}

inline void LP::Begin()
{
   repamMode_ = LPReparametrizationMode::Undefined;
   assert(f_.size() > 1); // otherwise we need not perform optimization: Just MaximizePotential f_[0]

   if(reparametrization_type_arg_.getValue() == "shared") {
     reparametrization_type_ = reparametrization_type::shared;
   } else if(reparametrization_type_arg_.getValue() == "residual") {
     reparametrization_type_ = reparametrization_type::residual;
   } else {
     assert(false);
   }

#ifdef LP_MP_PARALLEL
   omp_set_num_threads(num_lp_threads_arg_.getValue());
   if(debug()) { std::cout << "number of threads = " << num_lp_threads_arg_.getValue() << "\n"; }
#endif 
}

void LP::SortFactors(
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
    if(factor_mask[f1] && factor_mask[f2]) {
      g.addEdge(f1,f2);
    }
  }

  f_sorted = g.topologicalSort();
  //std::vector<INDEX> sortedIndices = g.topologicalSort();
  assert(f_sorted.size() == f_.size());

  const INDEX no_active_factors = std::count(factor_mask.begin(), factor_mask.end(), true);
  std::vector<FactorTypeAdapter*> fSorted;
  fSorted.reserve(no_active_factors);
  for(INDEX i=0; i<f_sorted.size(); i++) {
    if(factor_mask[ f_sorted[i] ]) {
      fSorted.push_back( f_[ f_sorted[i] ] );
    }
  }
  assert(fSorted.size() == no_active_factors);
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

void LP::SortFactors()
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


#ifdef LP_MP_PARALLEL
// a factor needs to be called with enabled synchronization only if one of its neighbots of distance 2 is updated by another thread
  template<typename ITERATOR>
inline std::vector<bool> LP::compute_synchronization(ITERATOR factor_begin, ITERATOR factor_end)
{
  const INDEX n = std::distance(factor_begin, factor_end);
  assert(n > 0);

  std::vector<INDEX> thread_number(this->f_.size(), std::numeric_limits<INDEX>::max());
  std::cout << "compute " << n << " factors to be synchronized\n";
#pragma omp parallel num_threads(num_lp_threads_arg_.getValue())
  {
    assert(num_lp_threads_arg_.getValue() == omp_get_num_threads());
    const int nthreads = num_lp_threads_arg_.getValue();
    const int ithread = omp_get_thread_num();
    assert(0 <= ithread && ithread < num_lp_threads_arg_.getValue());
    const int start = (ithread*n)/nthreads;
    const int finish = ((ithread+1)*n)/nthreads;

    for(INDEX i=start; i<finish; ++i) {
      const INDEX factor_number = factor_address_to_index_[*(factor_begin+i)];
      thread_number[factor_number] = ithread;
    }
  }

  // check for every factor all its neighbors and see whether more than two possible threads access it.
  std::vector<bool> conflict_factor(this->f_.size(), false);
#pragma omp parallel for
  for(INDEX i=0; i<this->f_.size(); ++i) {
    auto *f = f_[i];
    INDEX prev_adjacent_thread_number = thread_number[i];
    for(auto m_it=f->begin(); m_it!=f->end(); ++m_it) {
      const INDEX adjacent_factor_number = factor_address_to_index_[m_it.GetConnectedFactor()];
      const INDEX adjacent_thread_number = thread_number[adjacent_factor_number];
      if(adjacent_thread_number != std::numeric_limits<INDEX>::max()) {
        if(prev_adjacent_thread_number != std::numeric_limits<INDEX>::max() && adjacent_thread_number != prev_adjacent_thread_number) {
          conflict_factor[i] = true;
        }
        prev_adjacent_thread_number = adjacent_thread_number;
      }
    }
  }
  std::cout << "# conflict factors = " << std::count(conflict_factor.begin(), conflict_factor.end(), true) << "\n";

  // if a factor is adjacent to a conflict factor or is itself one, then it needs to be synchronized
  std::vector<bool> synchronize(n);
#pragma omp parallel for
  for(INDEX i=0; i<n; ++i) {
    auto* f = *(factor_begin+i);
    const INDEX factor_number = factor_address_to_index_[f];
    for(auto m_it=f->begin(); m_it!=f->end(); ++m_it) {
      const INDEX adjacent_factor_number = factor_address_to_index_[m_it.GetConnectedFactor()];
      if(conflict_factor[adjacent_factor_number]) {
        synchronize[i] = true;
      }
    }
    if(conflict_factor[factor_number]) {
      synchronize[i] = true;
    }
  }

  if(debug()) {
    std::cout << std::count(synchronize.begin(), synchronize.end(), true) << ";" << synchronize.size() << "\n";
    std::cout << "\%factors to synchronize = " << REAL(std::count(synchronize.begin(), synchronize.end(), true)) / REAL(synchronize.size()) << "\n";
  }
  return std::move(synchronize);
}
#endif

inline void LP::ComputePass(const INDEX iteration)
{
#ifdef LP_MP_PARALLEL
   compute_synchronization();
#endif
   ComputeForwardPass();
   ComputeBackwardPass();
}

void LP::ComputeForwardPass()
{
  const auto omega = get_omega();
#ifdef LP_MP_PARALLEL
  ComputePassSynchronized(forwardUpdateOrdering_.begin(), forwardUpdateOrdering_.end(), omega.forward.begin(), omega.forward.end(), synchronize_forward_.begin(), synchronize_forward_.end()); 
#else
  ComputePass(forwardUpdateOrdering_.begin(), forwardUpdateOrdering_.end(), omega.forward.begin(), omega.receive_mask_forward.begin()); 
#endif
}
void LP::ComputeBackwardPass()
{
  const auto omega = get_omega();
#ifdef LP_MP_PARALLEL
  ComputePassSynchronized(backwardUpdateOrdering_.begin(), backwardUpdateOrdering_.end(), omega.backward.begin(), omega.backward.end(), synchronize_backward_.begin(), synchronize_backward_.end()); 
#else
  ComputePass(backwardUpdateOrdering_.begin(), backwardUpdateOrdering_.end(), omega.backward.begin(), omega.receive_mask_backward.begin());
#endif
}

void LP::ComputeForwardPassAndPrimal(const INDEX iteration)
{
  const auto omega = get_omega();
#ifdef LP_MP_PARALLEL
  ComputePassAndPrimalSynchronized(forwardUpdateOrdering_.begin(), forwardUpdateOrdering_.end(), omega.forward.begin(), synchronize_forward_.begin(), 2*iteration+1); // timestamp must be > 0, otherwise in the first iteration primal does not get initialized
#else
  ComputePassAndPrimal(forwardUpdateOrdering_.begin(), forwardUpdateOrdering_.end(), omega.forward.begin(), omega.receive_mask_forward.begin(), 2*iteration+1); // timestamp must be > 0, otherwise in the first iteration primal does not get initialized
#endif
}
void LP::ComputeBackwardPassAndPrimal(const INDEX iteration)
{
  const auto omega = get_omega();
#ifdef LP_MP_PARALLEL
  ComputePassAndPrimalSynchronized(backwardUpdateOrdering_.begin(), backwardUpdateOrdering_.end(), omega.backward.begin(), synchronize_backward_.begin(), 2*iteration + 2); 
#else
  ComputePassAndPrimal(backwardUpdateOrdering_.begin(), backwardUpdateOrdering_.end(), omega.backward.begin(), omega.receive_mask_backward.begin(), 2*iteration + 2); 
#endif
}

void LP::ComputePassAndPrimal(const INDEX iteration)
{
  ComputeForwardPassAndPrimal(iteration);
  ComputeBackwardPassAndPrimal(iteration);
}

#ifdef LP_MP_PARALLEL
template<typename FACTOR_ITERATOR, typename OMEGA_ITERATOR, typename SYNCHRONIZATION_ITERATOR>
void LP::ComputePassSynchronized(
       FACTOR_ITERATOR factorIt, const FACTOR_ITERATOR factorItEnd, 
       OMEGA_ITERATOR omega_begin, OMEGA_ITERATOR omega_end,
       SYNCHRONIZATION_ITERATOR synchronization_begin, SYNCHRONIZATION_ITERATOR synchronization_end)

{
  const INDEX n = std::distance(factorIt, factorItEnd);
  assert(std::distance(factorIt, factorItEnd) == std::distance(omega_begin, omega_end));
  assert(std::distance(factorIt, factorItEnd) == std::distance(synchronization_begin, synchronization_end));

  //std::cout << "# synchronization calls = " << std::count(synchronization_begin, synchronization_end, true) << "\n";
  //for(INDEX i=0; i<n; ++i) {
  //  std::cout << i << ": " << (*(synchronization_begin + i) == true ? "true" : "false") << "\n";
  //}
#pragma omp parallel num_threads(num_lp_threads_arg_.getValue())
  {
    assert(num_lp_threads_arg_.getValue() == omp_get_num_threads());
    const int nthreads = num_lp_threads_arg_.getValue();
    const int ithread = omp_get_thread_num();
    assert(0 <= ithread && ithread < num_lp_threads_arg_.getValue());
    const int start = (ithread*n)/nthreads;
    const int finish = ((ithread+1)*n)/nthreads;

    for(INDEX i=start; i<finish; ++i) {
      auto* f = *(factorIt + i); 
      if(*(synchronization_begin+i)) {
        f->UpdateFactorSynchronized(*(omega_begin + i));
        //f->UpdateFactor(*(omega_begin + i));
      } else {
        f->UpdateFactor(*(omega_begin + i));
        //f->UpdateFactorSynchronized(*(omegaIt + i));
      }
    }
  } 
}
#endif

template<typename FACTOR_ITERATOR, typename OMEGA_ITERATOR, typename RECEIVE_MASK_ITERATOR>
void LP::ComputePass(FACTOR_ITERATOR factorIt, const FACTOR_ITERATOR factorItEnd, OMEGA_ITERATOR omegaIt, RECEIVE_MASK_ITERATOR receive_it)
{
  //assert(std::distance(factorItEnd, factorIt) == std::distance(omegaIt, omegaItEnd));
  const INDEX n = std::distance(factorIt, factorItEnd);
//#pragma omp parallel for schedule(static)
  if(reparametrization_type_ == reparametrization_type::shared) {
    for(INDEX i=0; i<n; ++i) {
      auto* f = *(factorIt + i);
      f->UpdateFactor(*(omegaIt + i), *(receive_it + i));
    }
  } else {
    assert(reparametrization_type_ == reparametrization_type::residual);
    for(INDEX i=0; i<n; ++i) {
      auto* f = *(factorIt + i);
      f->update_factor_residual(*(omegaIt + i));
    }
  }
}

bool LP::omega_valid(const weight_array& omega) const
{
    for(std::size_t i=0; i<omega.size(); ++i) {
        assert(*std::min_element(omega[i].begin(), omega[i].end()) >= 0.0);
        assert(std::accumulate(omega[i].begin(), omega[i].end(), 0.0) <= 1.0 + eps);
    }
}

inline void LP::ComputeAnisotropicWeights()
{
  ComputeAnisotropicWeights(forwardOrdering_.begin(), forwardOrdering_.end(), f_forward_sorted_.begin(), f_forward_sorted_.end(), factor_mask.begin(), factor_mask.end(), omegaForwardAnisotropic_, anisotropic_receive_mask_forward_);
  ComputeAnisotropicWeights(backwardOrdering_.begin(), backwardOrdering_.end(), f_backward_sorted_.begin(), f_backward_sorted_.end(), factor_mask.begin(), factor_mask.end(), omegaBackwardAnisotropic_, anisotropic_receive_mask_backward_);

  omega_valid(omegaForwardAnisotropic_);
  omega_valid(omegaBackwardAnisotropic_);
}

inline void LP::ComputeAnisotropicWeights2()
{
  ComputeAnisotropicWeights2(forwardOrdering_.begin(), forwardOrdering_.end(), f_forward_sorted_.begin(), f_forward_sorted_.end(), omegaForwardAnisotropic2_);
  ComputeAnisotropicWeights2(backwardOrdering_.begin(), backwardOrdering_.end(), f_backward_sorted_.begin(), f_backward_sorted_.end(), omegaBackwardAnisotropic2_);

  omega_valid(omegaForwardAnisotropic2_);
  omega_valid(omegaBackwardAnisotropic2_);
}

inline void LP::ComputeUniformWeights()
{
  ComputeUniformWeights(forwardOrdering_.begin(), forwardOrdering_.end(), factor_mask.begin(), factor_mask.end(), omegaForwardIsotropic_, 0.0);
  ComputeUniformWeights(backwardOrdering_.begin(), backwardOrdering_.end(), factor_mask.begin(), factor_mask.end(), omegaBackwardIsotropic_, 0.0);

  omega_valid(omegaForwardIsotropic_);
  omega_valid(omegaBackwardIsotropic_);

  assert(this->backwardUpdateOrdering_.size() == omegaBackwardIsotropic_.size());
  for(auto it = this->backwardUpdateOrdering_.begin(); it != this->backwardUpdateOrdering_.end(); ++it) {
    assert((*it)->no_send_messages() == omegaBackwardIsotropic_[ std::distance(this->backwardUpdateOrdering_.begin(), it) ].size());
  } 

  for(auto it = this->forwardUpdateOrdering_.begin(); it != this->forwardUpdateOrdering_.end(); ++it) {
    assert((*it)->no_send_messages() == omegaForwardIsotropic_[ std::distance(this->forwardUpdateOrdering_.begin(), it) ].size());
  } 
}

inline void LP::ComputeDampedUniformWeights()
{
  ComputeUniformWeights(forwardOrdering_.begin(), forwardOrdering_.end(), factor_mask.begin(), factor_mask.end(), omegaForwardIsotropicDamped_, 1.0);
  ComputeUniformWeights(backwardOrdering_.begin(), backwardOrdering_.end(), factor_mask.begin(), factor_mask.end(), omegaBackwardIsotropicDamped_, 1.0);

  omega_valid(omegaForwardIsotropicDamped_);
  omega_valid(omegaBackwardIsotropicDamped_);
}

// Here we check whether messages constraints are satisfied
inline bool LP::CheckPrimalConsistency() const
{
   volatile bool consistent=true; // or use std::atomic<bool>?

#pragma omp parallel for shared(consistent)
   for(INDEX i=0; i<f_.size(); ++i) {
      if(!consistent) continue;
      if(!f_[i]->check_primal_consistency()) {
          consistent = false;
      }
   }
   if(debug()) { std::cout << "primal solution consistent: " << (consistent ? "true" : "false") << "\n"; }
   return consistent;
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
      auto* f_left = m_[i].left;
      const INDEX f_index_left = factor_address_to_index_[f_left];
      const INDEX index_left = f_sorted_inverse[f_index_left];
      auto* f_right = m_[i].right;
      const INDEX f_index_right = factor_address_to_index_[f_right];
      const INDEX index_right = f_sorted_inverse[f_index_right];
      
      if(m_[i].sends_message_to_right && index_left < index_right) {
        no_send_messages_later[index_left]++;
      }
      if(m_[i].sends_message_to_left && index_right < index_left) {
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
               if(mIt.SendsMessage()) {
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
template<typename FACTOR_ITERATOR, typename FACTOR_SORT_ITERATOR, typename FACTOR_MASK_ITERATOR>
void LP::ComputeAnisotropicWeights(
      FACTOR_ITERATOR factorIt, FACTOR_ITERATOR factorEndIt, // sorted pointers to factors
      FACTOR_SORT_ITERATOR factor_sort_begin, FACTOR_SORT_ITERATOR factor_sort_end, // sorted factor indices in f_
      FACTOR_MASK_ITERATOR factor_mask_begin, FACTOR_MASK_ITERATOR factor_mask_end,
      weight_array& omega, receive_array& receive_mask)
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
   std::vector<INDEX> no_send_factors(f_.size(),0); // no of messages over which messages can be sent from given factor
   std::vector<INDEX> no_send_factors_later(f_.size(),0); // no of messages over which messages are sent from given factor in given direction
   std::vector<INDEX> no_receiving_factors_later(f_.size(),0); // number of factors later than current one receiving message from current one
   std::vector<INDEX> last_receiving_factor(f_.size(), 0); // what is the last (in the order given by factor iterator) factor that receives a message?
   std::vector<INDEX> first_receiving_factor(f_.size(), 0); // what is the last (in the order given by factor iterator) factor that receives a message?
#endif

   // do zrobienia: if factor is not visited at all, then omega is not needed for that entry. We must filter out such entries still
#pragma omp parallel for
   for(INDEX i=0; i<m_.size(); ++i) {
      auto* f_left = m_[i].left;
      const INDEX f_index_left = factor_address_to_index_[f_left];
      const INDEX index_left = f_sorted_inverse[f_index_left];
      auto* f_right = m_[i].right;
      const INDEX f_index_right = factor_address_to_index_[f_right];
      const INDEX index_right = f_sorted_inverse[f_index_right];
      
      if(!(factor_mask_begin[f_index_left] && factor_mask_begin[f_index_right])) {
        continue; 
      }
      
      if(m_[i].receives_message_from_left) {
         if(index_left < index_right) {
            no_receiving_factors_later[index_left]++;
         }
#ifdef LP_MP_PARALLEL
         INDEX old_val = last_receiving_factor[index_left];
         const INDEX new_val = std::max(old_val, index_right);
         while(old_val < new_val && !last_receiving_factor[index_left].compare_exchange_weak(old_val, new_val)) ;
#else
         last_receiving_factor[index_left] = std::max(last_receiving_factor[index_left], index_right);
         first_receiving_factor[index_left] = std::min(first_receiving_factor[index_left], index_right);
#endif
      }

      if(m_[i].receives_message_from_right) {
         if(index_left > index_right) {
            no_receiving_factors_later[index_right]++;
         }
#ifdef LP_MP_PARALLEL
         INDEX old_val = last_receiving_factor[index_right];
         const INDEX new_val = std::max(old_val, index_left);
         while(old_val < new_val && !last_receiving_factor[index_right].compare_exchange_weak(old_val, new_val)) ;
#else
         last_receiving_factor[index_right] = std::max(last_receiving_factor[index_right], index_left);
         first_receiving_factor[index_right] = std::min(first_receiving_factor[index_right], index_left);
#endif
      }
   }

#pragma omp parallel for
   for(INDEX i=0; i<m_.size(); ++i) {
      auto* f_left = m_[i].left;
      const INDEX f_index_left = factor_address_to_index_[f_left];
      const INDEX index_left = f_sorted_inverse[f_index_left];
      auto* f_right = m_[i].right;
      const INDEX f_index_right = factor_address_to_index_[f_right];
      const INDEX index_right = f_sorted_inverse[f_index_right];

      if(!(factor_mask_begin[f_index_left] && factor_mask_begin[f_index_right])) {
        continue; 
      }

      if(m_[i].sends_message_to_right) {
         no_send_factors[index_left]++;
         if(index_left < index_right || last_receiving_factor[index_right] > index_left) {
            no_send_factors_later[index_left]++;
         }
      }
      if(m_[i].sends_message_to_left) {
         no_send_factors[index_right]++;
         if(index_right < index_left || last_receiving_factor[index_left] > index_right) {
            no_send_factors_later[index_right]++;
         }
      }
   }

   std::vector<INDEX> omega_size;
   omega_size.reserve(f_.size());
   std::vector<INDEX> receive_size;
   receive_size.reserve(f_.size());
   for(auto it=factorIt; it!=factorEndIt; ++it) {
       const INDEX f_index = factor_address_to_index_[*it];
       if(!factor_mask_begin[f_index]) {
           continue;
       }
       if((*it)->FactorUpdated()) {
           omega_size.push_back( (*it)->no_send_messages() );
           receive_size.push_back( (*it)->no_receive_messages() );
       }
   }

   omega = weight_array(omega_size);
   receive_mask = receive_array(receive_size);

   {
      INDEX c=0;
      for(auto it=factorIt; it!=factorEndIt; ++it) {
         const INDEX f_index = factor_address_to_index_[*it];
         if(!factor_mask_begin[f_index]) {
           continue;
         }
         const INDEX i = std::distance(factorIt, it);
         assert(i == f_sorted_inverse[ factor_address_to_index_[*it] ]);
         if((*it)->FactorUpdated()) {
            INDEX k_send = 0;
            INDEX k_receive = 0;
            for(auto mIt=(*it)->begin(); mIt!=(*it)->end(); ++mIt) {
                auto* f_connected = mIt.GetConnectedFactor();
                const INDEX j_index = factor_address_to_index_[f_connected];
                const INDEX j = f_sorted_inverse[ j_index ];

               if(mIt.SendsMessage()) {
                  assert(i != j);
                  if(i<j || last_receiving_factor[j] > i && factor_mask_begin[j_index]) {
                    omega[c][k_send] = (1.0/REAL(no_receiving_factors_later[i] + std::max(INDEX(no_send_factors_later[i]), INDEX(no_send_factors[i]) - INDEX(no_send_factors_later[i]))));
                    //if(no_receiving_factors_later[i] > 0) {
                    //  omega[c][k] = (1.0/REAL(1 + std::max(INDEX(no_send_factors_later[i]), INDEX(no_send_factors[i]) - INDEX(no_send_factors_later[i]))));
                    //} else {
                    //  omega[c][k] = (1.0/REAL(no_send_factors_later[i]));
                    //}
                  } else {
                     omega[c][k_send] = 0.0;
                  } 
                  ++k_send;
               }
               if(mIt.ReceivesMessage()) {
                   if(j<i && first_receiving_factor[i] < j) {
                       receive_mask[c][k_receive] = true; 
                   } else {
                       receive_mask[c][k_receive] = false; 
                   }
                   k_receive++; 
               }
            }
            assert(std::accumulate(omega[c].begin(), omega[c].end(), 0.0) <= 1.0 + eps);
            ++c;
         }
      }
   }

   // check whether all messages were added to m_. Possibly, this can be automated: Traverse all factors, get all messages, add them to m_ and avoid duplicates along the way.
   //std::cout << 2*m_.size() << " == "<< std::accumulate(f_.begin(), f_.end(), 0, [](INDEX sum, auto* f){ return sum + f->no_messages(); }) << "\n";
   assert(2*m_.size() == std::accumulate(f_.begin(), f_.end(), 0, [](INDEX sum, auto* f){ return sum + f->no_messages(); }));
   for(INDEX i=0; i<omega.size(); ++i) {
      //assert(omega[i].size() <= (*(factorIt+i))->no_messages());
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
template<typename FACTOR_ITERATOR, typename FACTOR_MASK_ITERATOR>
void LP::ComputeUniformWeights(
    FACTOR_ITERATOR factorIt, FACTOR_ITERATOR factorEndIt,
    FACTOR_MASK_ITERATOR factor_mask_begin, FACTOR_MASK_ITERATOR factor_mask_end,
    weight_array& omega, const REAL leave_weight)
{
   assert(leave_weight >= 0.0 && leave_weight <= 1.0);
   assert(factorEndIt - factorIt == f_.size());

   std::vector<INDEX> omega_size(std::distance(factorIt, factorEndIt),0);
   INDEX c=0;
   for(auto it=factorIt; it != factorEndIt; ++it) {
     const INDEX f_index = factor_address_to_index_[ *it ];
      if((*it)->FactorUpdated() && factor_mask_begin[f_index]) {
        omega_size[c] = (*it)->no_send_messages();
        ++c;
      }
   }

   omega_size.resize(c);
   omega = two_dim_variable_array<REAL>(omega_size);

   assert(omega.size() == omega_size.size());
   for(INDEX i=0; i<omega.size(); ++i) {
      assert(omega[i].size() == omega_size[i]);
   }

   c=0;
   for(auto it=factorIt; it != factorEndIt; ++it) {
     const INDEX f_index = factor_address_to_index_[ *it ];
     if((*it)->FactorUpdated() && factor_mask_begin[f_index]) {
       INDEX k=0;
       for(auto mIt=(*it)->begin(); mIt!=(*it)->end(); ++mIt) {
         if(mIt.SendsMessage()) {
           auto* f_connected = mIt.GetConnectedFactor();
           const INDEX connected_index = factor_address_to_index_[f_connected];
           if(factor_mask_begin[connected_index]) {
             omega[c][k] = (1.0/REAL(omega[c].size() + leave_weight));
           } else {
             omega[c][k] = 0.0;
           }
           ++k;
         }
       }
       assert(k == omega[c].size());
       c++;
     }
   }
   assert(c == omega.size());
}

// compute anisotropic and damped uniform weights, then average them
void LP::ComputeMixedWeights()
{
  assert(false); // check for valid flags!
  ComputeDampedUniformWeights();
  ComputeAnisotropicWeights();
  ComputeMixedWeights(omegaForwardAnisotropic_, omegaForwardIsotropicDamped_, omegaForwardMixed_);
  ComputeMixedWeights(omegaBackwardAnisotropic_, omegaBackwardIsotropicDamped_, omegaBackwardMixed_);
} 

inline void LP::ComputeMixedWeights(
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

void LP::compute_full_receive_mask()
{
  compute_full_receive_mask(forwardOrdering_.begin(), forwardOrdering_.end(), full_receive_mask_forward_);
  compute_full_receive_mask(backwardOrdering_.begin(), backwardOrdering_.end(), full_receive_mask_backward_); 
}

template<typename FACTOR_ITERATOR>
void LP::compute_full_receive_mask(FACTOR_ITERATOR factor_begin, FACTOR_ITERATOR factor_end, receive_array& receive_mask)
{
    std::vector<INDEX> mask_size;
    mask_size.reserve(std::distance(factor_begin, factor_end));
    for(auto it=factor_begin; it!=factor_end; ++it) {
        if((*it)->FactorUpdated()) {
            mask_size.push_back( (*it)->no_receive_messages() );
        } 
    }
    receive_mask = receive_array(mask_size); 
    for(std::size_t i=0; i<receive_mask.size(); ++i) {
        for(std::size_t j=0; j<receive_mask[i].size(); ++j) {
            receive_mask(i,j) = true;
        }
    }

}

double LP::LowerBound() const
{
  double lb = constant_;
#pragma omp parallel for reduction(+:lb)
  for(INDEX i=0; i<f_.size(); ++i) {
    lb += f_[i]->LowerBound();
    assert( f_[i]->LowerBound() > -10000000.0);
    assert(std::isfinite(lb));
  }
  return lb;
}

double LP::EvaluatePrimal() {
  const bool consistent = CheckPrimalConsistency();
  if(consistent == false) return std::numeric_limits<REAL>::infinity();

  double cost = constant_;
#pragma omp parallel for reduction(+:cost)
  for(INDEX i=0; i<f_.size(); ++i) {
    assert(f_[i]->LowerBound() <= f_[i]->EvaluatePrimal() + eps);
    cost += f_[i]->EvaluatePrimal();
  }

  if(debug()) { std::cout << "primal cost = " << cost << "\n"; }

  return cost;
}


template<typename FACTOR_ITERATOR, typename OMEGA_ITERATOR, typename RECEIVE_MASK_ITERATOR>
void LP::ComputePassAndPrimal(FACTOR_ITERATOR factorIt, const FACTOR_ITERATOR factorEndIt, OMEGA_ITERATOR omegaIt, RECEIVE_MASK_ITERATOR receive_mask_it, INDEX iteration)
{
   //possibly do not use parallelization here
//#pragma omp parallel for schedule(static)
   for(INDEX i=0; i<std::distance(factorIt, factorEndIt); ++i) {
      auto* f = *(factorIt+i);
      f->UpdateFactorPrimal(*(omegaIt + i), *(receive_mask_it + i), iteration);
   }
}

#ifdef LP_MP_PARALLEL
template<typename FACTOR_ITERATOR, typename OMEGA_ITERATOR, typename SYNCHRONIZATION_ITERATOR>
void LP::ComputePassAndPrimalSynchronized(FACTOR_ITERATOR factorIt, const FACTOR_ITERATOR factorEndIt, OMEGA_ITERATOR omegaIt, SYNCHRONIZATION_ITERATOR synchronization_begin, INDEX iteration)
{
   //possibly do not use parallelization here
#pragma omp parallel for schedule(static)
  for(INDEX i=0; i<std::distance(factorIt, factorEndIt); ++i) {
    auto* f = *(factorIt + i);
    if(*(synchronization_begin+i)) {
      f->UpdateFactorPrimalSynchronized(*(omegaIt + i), iteration);
    } else {
      f->UpdateFactorPrimal(*(omegaIt + i), iteration);
    }
  }
}
#endif

void LP::set_flags_dirty()
{
  ordering_valid_ = false;
  omega_anisotropic_valid_ = false;
  omega_anisotropic2_valid_ = false;
  omega_isotropic_valid_ = false;
  omega_isotropic_damped_valid_ = false;
  omega_mixed_valid_ = false;
#ifdef LP_MP_PARALLEL
  synchronization_valid_ = false;
#endif
}

std::vector<bool> LP::get_inconsistent_mask(const std::size_t no_fatten_rounds)
{
  std::vector<bool> inconsistent_mask(f_.size(),false);
  
  // check for locally non-optimal factors
  for(std::size_t i=0; i<f_.size(); ++i) {
    auto* f = f_[i];
    assert(f->EvaluatePrimal() < std::numeric_limits<REAL>::infinity());
    if(f->LowerBound() < f->EvaluatePrimal() - eps) {
      inconsistent_mask[i] = true;
    }
  }

  // check for violated messages
  for(auto* f : f_) {
      if(!f->check_primal_consistency()) {
          auto f_index = factor_address_to_index_[f];
          inconsistent_mask[f_index] = true;
      }
  }

  // fatten the region
  auto fatten = [&]() {
    for(auto m : m_) {
      auto* l = m.left;
      auto l_index = factor_address_to_index_[l];
      auto* r = m.right;
      auto r_index = factor_address_to_index_[r];

      if(inconsistent_mask[l_index] == true || inconsistent_mask[r_index] == true) {
        inconsistent_mask[l_index] = true;
        inconsistent_mask[r_index] = true;
      }
    }
  };

  for(INDEX iter=0; iter<no_fatten_rounds; ++iter) {
    fatten();
  }

  if(debug()) {
    std::cout << "\% inconsistent factors = " << std::count(inconsistent_mask.begin(), inconsistent_mask.end(), true)/REAL(f_.size()) << "\n";
  }

  return inconsistent_mask;
}

template<typename FACTOR_ITERATOR, typename FACTOR_MASK_ITERATOR>
std::vector<FactorTypeAdapter*> LP::get_masked_factors(
    FACTOR_ITERATOR factor_begin, FACTOR_ITERATOR factor_end,
    FACTOR_MASK_ITERATOR factor_mask_begin, FACTOR_MASK_ITERATOR factor_mask_end)
{
  assert(std::distance(factor_mask_begin, factor_mask_end) == f_.size());
  
  std::vector<FactorTypeAdapter*> factors;
  for(auto f_it=factor_begin; f_it!=factor_end; ++f_it) {
    const auto f_index = factor_address_to_index_[*f_it];
    if(factor_mask_begin[f_index]) {
      factors.push_back(*f_it);
    }
  }

  return factors;
}

inline void LP::reduce_optimization_factors()
{
  factor_mask = get_inconsistent_mask();
  set_flags_dirty();
}


} // end namespace LP_MP

#endif // LP_MP_MAIN


