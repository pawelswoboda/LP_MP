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
   virtual INDEX size() const = 0;
   virtual INDEX PrimalSize() const = 0;
   MessageIterator begin(); 
   MessageIterator end();
   virtual INDEX GetNoMessages() const = 0;
   virtual MessageTypeAdapter* GetMessage(const INDEX n) const = 0;
   virtual FactorTypeAdapter* GetConnectedFactor(const INDEX i) const = 0;
   virtual bool CanSendMessage(const INDEX i) const = 0;
   virtual REAL LowerBound() const = 0;
   virtual std::vector<REAL> GetReparametrizedPotential() const = 0;
   virtual void init_primal() = 0;
   virtual void MaximizePotentialAndComputePrimal() = 0;
   //virtual PrimalSolutionStorageAdapter* AllocatePrimalSolutionStorage() const = 0;
   //virtual bool CanComputePrimalSolution() const = 0;
   // the offset in the primal storage
   virtual void SetPrimalOffset(const INDEX) = 0; // do zrobienia: delete
   virtual INDEX GetPrimalOffset() const = 0; // do zrobienia: delete
   
   virtual void SetAuxOffset(const INDEX n) = 0; // do zrobienia: delete
   virtual INDEX GetAuxOffset() const = 0; // do zrobienia: delete
   
   // do zrobienia: this function is not needed. Evaluation can be performed automatically
   virtual REAL EvaluatePrimal() const = 0;
   // do zrobienia: this is not needed as well and could be automated. Possibly it is good to keep this to enable solution rewriting.
   //virtual void WritePrimal(PrimalSolutionStorage::Element primalSolution, std::ofstream& fs) const = 0;

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
   {
      std::cout << "kwaskwaskwas3\n";
   }
   ~LP() 
   {
      for(INDEX i=0; i<m_.size(); i++) { delete m_[i]; }
      for(INDEX i=0; i<f_.size(); i++) { delete f_[i]; }
   }

   // make a deep copy of factors and messages. Adjust pointers to messages and factors
   LP(const LP& o) 
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
      omega_isotropic_valid_ = o.omega_isotropic_valid_ ;
      omega_isotropic_damped_valid_ = o.omega_isotropic_damped_valid_ ;

      omegaForwardAnisotropic_ = o.omegaForwardAnisotropic_; 
      omegaBackwardAnisotropic_ = o.omegaBackwardAnisotropic_;
      omegaForwardIsotropic_ = o.omegaForwardIsotropic_; 
      omegaBackwardIsotropic_ = o.omegaBackwardIsotropic_;
      omegaForwardIsotropicDamped_ = o.omegaForwardIsotropicDamped_; 
      omegaBackwardIsotropicDamped_ = o.omegaBackwardIsotropicDamped_;

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
      //factorRel_.push_back(std::make_pair(f1,f2));
      ForwardPassFactorRelation(f1,f2);
      BackwardPassFactorRelation(f1,f2);
   }

   void ForwardPassFactorRelation(FactorTypeAdapter* f1, FactorTypeAdapter* f2) { set_flags_dirty(); assert(f1!=f2); forward_pass_factor_rel_.push_back({f1,f2}); }
   void BackwardPassFactorRelation(FactorTypeAdapter* f1, FactorTypeAdapter* f2) { set_flags_dirty(); assert(f1!=f2); backward_pass_factor_rel_.push_back({f1,f2}); }


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
         std::vector<FactorTypeAdapter*>& update_ordering
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

      f_sorted_ = g.topologicalSort();
      //std::vector<INDEX> sortedIndices = g.topologicalSort();
      assert(f_sorted_.size() == f_.size());

      std::vector<FactorTypeAdapter*> fSorted;
      fSorted.reserve(f_.size());
      for(INDEX i=0; i<f_sorted_.size(); i++) {
         fSorted.push_back( f_[ f_sorted_[i] ] );//indexToFactor[sortedIndices[i]] );
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
      // check whether sorting was suffessful
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

      SortFactors(forward_pass_factor_rel_, forwardOrdering_, forwardUpdateOrdering_);
      SortFactors(backward_pass_factor_rel_, backwardOrdering_, backwardUpdateOrdering_);
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
   void ComputeUniformWeights();
   void ComputeDampedUniformWeights();
   template<typename FACTOR_ITERATOR>
   void ComputeUniformWeights(FACTOR_ITERATOR factorIt, FACTOR_ITERATOR factorEndIt, two_dim_variable_array<REAL>& omega, const REAL leave_weight); // do zrobienia: rename to isotropic weights

   REAL LowerBound() const
   {
      REAL lb = 0.0;
      for(auto fIt=f_.begin(); fIt!=f_.end(); fIt++) {
         lb += (*fIt)->LowerBound();
         assert( (*fIt)->LowerBound() > -10000000.0);
         assert(std::isfinite(lb));
      }
      return lb;
   }

   void InitializePrimalVector(PrimalSolutionStorage& p) { InitializePrimalVector(f_.begin(), f_.end(), p); }
   template<typename FACTOR_ITERATOR>
      void InitializePrimalVector(FACTOR_ITERATOR factorIt, const FACTOR_ITERATOR factorEndIt, PrimalSolutionStorage& v);
   bool CheckPrimalConsistency() const
   {
      return CheckPrimalConsistency(f_.begin(), f_.end());
   }
   template<typename FACTOR_ITERATOR>
      bool CheckPrimalConsistency(FACTOR_ITERATOR factorIt, const FACTOR_ITERATOR factorEndIt) const;

   template<typename FACTOR_ITERATOR>
   REAL EvaluatePrimal(cereal::BinaryInputArchive& primal, FACTOR_ITERATOR factor_begin, FACTOR_ITERATOR factor_end) // read in primal solution from primal archive and evaluate cost
   {
      assert(std::distance(factor_begin, factor_end) == f_.size());

   }
   REAL EvaluatePrimal() {
      return EvaluatePrimal(f_.begin(), f_.end());
   }
   template<typename FACTOR_ITERATOR>
   REAL EvaluatePrimal(FACTOR_ITERATOR factorIt, const FACTOR_ITERATOR factorEndIt) const;

   void UpdateFactor(FactorTypeAdapter* f, const weight_vector& omega) // perform one block coordinate step for factor f
   {
      f->UpdateFactor(omega);
   }
   void ComputePass();

   void ComputeForwardPassAndPrimal(const INDEX iteration)
   {
      const auto omega = get_omega();
      ComputePassAndPrimal(forwardUpdateOrdering_.begin(), forwardUpdateOrdering_.end(), omega.forward.begin(), 2*iteration+1); // timestamp must be > 0, otherwise in the first iteration primal does not get initialized
      const REAL forward_cost = EvaluatePrimal();
      std::cout << "forward cost = " << forward_cost << "\n";
   }
   void ComputeBackwardPassAndPrimal(const INDEX iteration)
   {
      const auto omega = get_omega();
      ComputePassAndPrimal(forwardUpdateOrdering_.rbegin(), forwardUpdateOrdering_.rend(), omega.backward.begin(), 2*iteration + 2); 
      const REAL backward_cost = EvaluatePrimal();
      std::cout << "backward cost = " << backward_cost << "\n";
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
      omega_isotropic_valid_ = false;
      omega_isotropic_damped_valid_ = false;
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
      } else if(repamMode_ == LPReparametrizationMode::Uniform) {
         ComputeUniformWeights();
         return omega_storage{omegaForwardIsotropic_, omegaBackwardIsotropic_};
      } else if(repamMode_ == LPReparametrizationMode::DampedUniform) {
         ComputeDampedUniformWeights();
         return omega_storage{omegaForwardIsotropicDamped_, omegaBackwardIsotropicDamped_};
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
   bool omega_isotropic_valid_ = false;
   two_dim_variable_array<REAL> omegaForwardIsotropic_, omegaBackwardIsotropic_;
   bool omega_isotropic_damped_valid_ = false;
   two_dim_variable_array<REAL> omegaForwardIsotropicDamped_, omegaBackwardIsotropicDamped_;

   std::vector<std::pair<FactorTypeAdapter*, FactorTypeAdapter*> > forward_pass_factor_rel_, backward_pass_factor_rel_; // factor ordering relations. First factor must come before second factor. factorRel_ must describe a DAG

   
   std::unordered_map<FactorTypeAdapter*,INDEX> factor_address_to_index_;
   std::vector<INDEX> f_sorted_; // sorted indices in factor vector f_ 

   LPReparametrizationMode repamMode_ = LPReparametrizationMode::Undefined;
};

template<typename BASE_LP_CLASS>
class LP_concurrent : public BASE_LP_CLASS {
public:
  LP_concurrent(TCLAP::CmdLine& cmd) 
    : BASE_LP_CLASS(cmd),
    num_lp_threads_arg_("","numLpThreads","number of threads for message passing, default = 1",false,1,&positiveIntegerConstraint,cmd)
  {}
  void Begin()
  {
    threads_ = decltype(threads_)(num_lp_threads_arg_.getValue());
    std::cout << "number of threads = " << threads_.size() << "\n";
    BASE_LP_CLASS::Begin();
  }

  //LP_concurrent(const INDEX no_threads = std::thread::hardware_concurrency())
  //  : threads_(no_threads)
  //{
  //  std::cout << "number of threads = " << threads_.size() << "\n";
  //  assert(no_threads >= 1);
  //}

  template<typename FACTOR_ITERATOR, typename OMEGA_ITERATOR>
   void ComputePass(FACTOR_ITERATOR factorIt, const FACTOR_ITERATOR factorItEnd, OMEGA_ITERATOR omegaIt)
   {

     // idea for load-balancing: measure time a thread needs to finish. Adjust numbers of factors given to each thread based on this.

     auto worker = [this] (auto factor_begin, auto factor_end, auto omega_it, const INDEX allocator_index ) {
       stack_allocator_index = allocator_index;
       for(; factor_begin!=factor_end; ++factor_begin, ++omega_it) {
         this->UpdateFactor(*factor_begin, *omega_it);
       }
     };

     const int grainsize = std::distance(factorIt, factorItEnd) / threads_.size();

     INDEX c=0;
     for(auto it = std::begin(threads_); it != std::end(threads_) - 1; ++it) {
       *it = std::thread(worker, factorIt, factorIt + grainsize, omegaIt, c%global_real_block_allocator_array.size());
       factorIt += grainsize;
       omegaIt += grainsize;
       ++c;
     }
     threads_.back() = std::thread(worker, factorIt, factorItEnd, omegaIt, c%global_real_block_allocator_array.size());

     for(auto&& i : threads_) {
       i.join();
     } 
   }

  void ComputePass()
  {
     const auto omega = this->get_omega();
     this->ComputePass(this->forwardUpdateOrdering_.begin(), this->forwardUpdateOrdering_.end(), omega.forward.begin());
     this->ComputePass(this->forwardUpdateOrdering_.rbegin(), this->forwardUpdateOrdering_.rend(), omega.backward.begin());
  }

private:

  template<typename LAMBDA, typename FACTOR_ITERATOR>
  void iterate_over_factors(LAMBDA& f, FACTOR_ITERATOR factor_begin, FACTOR_ITERATOR factor_end)
  {
     const int grainsize = std::distance(factor_begin, factor_end) / threads_.size();

     for(auto it = std::begin(threads_); it != std::end(threads_) - 1; ++it) {
       *it = std::thread(f, factor_begin, factor_begin + grainsize);
       factor_begin += grainsize;
     }
     threads_.back() = std::thread(f, factor_begin, factor_end);

     for(auto&& i : threads_) {
       i.join();
     }
  }

  mutable std::vector<std::thread> threads_; // thread pool to be reused

  TCLAP::ValueArg<INDEX> num_lp_threads_arg_;
};

template<typename BASE_LP_CLASS>
class LP_sat : public BASE_LP_CLASS
{
private:
  struct sat_th {
    sat_th() 
      : th(0.01),
      feasible_bound(0.1),
      infeasible_bound(0.0)
    {}

    void adjust_th(const bool feasible)
    {
      assert(infeasible_bound <= th && th <= feasible_bound);
      if(feasible) {
        if(th >= feasible_bound - 0.01) {
          infeasible_bound *= 0.8;
        }
        feasible_bound = std::min(feasible_bound, th);
        th += (infeasible_bound - th)/2.0; 
      } else {
        if(th <= infeasible_bound + 0.001) {
          feasible_bound *= 2.0;
        }
        infeasible_bound = std::max(infeasible_bound, th);
        th += (feasible_bound - th)/2.0;
      } 
    }

    REAL th, feasible_bound, infeasible_bound;
  };


public:
  using LP_type = LP_sat<BASE_LP_CLASS>;

  LP_sat(TCLAP::CmdLine& cmd) : BASE_LP_CLASS(cmd)
  {
     std::cout << "kwaskwaskwas1\n";
     sat_ = lglinit();
     assert(sat_ != nullptr);
     //sat_.set_no_simplify(); // seems to make solver much faster
  }

  ~LP_sat()
  {
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


  virtual INDEX AddFactor(FactorTypeAdapter* f) 
  {
     sat_var_.push_back( lglmaxvar(sat_) );
     //sat_var_.push_back(sat_.nVars());
     //std::cout << "number of variables in sat_ = " << sat_var_[sat_var_.size()-1] << "\n";
     f->construct_sat_clauses(sat_);
     INDEX n = BASE_LP_CLASS::AddFactor(f);
     return n;
  }

   virtual INDEX AddMessage(MessageTypeAdapter* m)
   {
      assert(this->factor_address_to_index_.find(m->GetLeftFactor()) != this->factor_address_to_index_.end());
      assert(this->factor_address_to_index_.find(m->GetRightFactor()) != this->factor_address_to_index_.end());
      const INDEX left_factor_number = this->factor_address_to_index_[m->GetLeftFactor()];
      const INDEX right_factor_number = this->factor_address_to_index_[m->GetRightFactor()];
      m->construct_sat_clauses(sat_, sat_var_[left_factor_number], sat_var_[right_factor_number]);

      return BASE_LP_CLASS::AddMessage(m);
   }

   static bool solve_sat_problem(LP_type* c, sat_vec<sat_literal> assumptions, sat_th* th)
   {
     for(INDEX i=0; i<assumptions.size(); ++i) {
       lglassume( c->sat_, assumptions[i] );
     }
     for(INDEX i=0; i<lglmaxvar(c->sat_); ++i) {
       lglfreeze(c->sat_, to_literal(i));
     }
     const int sat_ret = lglsat(c->sat_);
     std::cout << "solved sat " << sat_ret << "\n";

     const bool feasible = sat_ret == LGL_SATISFIABLE;
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
         compute_pass_reduce_sat(this->forwardUpdateOrdering_.rbegin(), this->forwardUpdateOrdering_.rend(), omega.backward.begin(), backward_sat_th_);
         cur_sat_reduction_direction_ = Direction::forward; 
      } else {
         this->ComputePass(this->backwardUpdateOrdering_.begin(), this->backwardUpdateOrdering_.end(), omega.backward.begin()); 
      }
   }
   void ComputePassAndPrimal(const INDEX iteration)
   {
      assert(false);
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
      sat_vec<sat_literal> assumptions;
      for(auto it=factor_begin; it!=factor_end; ++it, ++omega_begin) {
         const INDEX factor_number = this->factor_address_to_index_[*it];
         (*it)->UpdateFactorSAT(*omega_begin, th.th, sat_var_[factor_number], assumptions);
      }

      // run sat solver on reduced problem asynchronuously
      if(!sat_handle_.valid()) { 
         std::cout << "start sat calculation\n";
         sat_handle_ = std::async(std::launch::async, solve_sat_problem, this, assumptions, &th);
      }

      std::cout << "collect sat result\n";
      const bool feasible = sat_handle_.get();
      if(feasible) {
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
         REAL primal_cost = this->EvaluatePrimal(this->f_.begin(), this->f_.end());
         std::cout << "sat solution cost = " << primal_cost << ", sat threshold = " << th.th << "\n"; 
      } else {
         std::cout << "sat not feasible with current threshold = " << th.th << "\n";
      }

      std::cout << "restart sat calculation\n";
      sat_handle_ = std::async(std::launch::async, solve_sat_problem, this, assumptions, &th);
   }
private:

   std::vector<sat_var> sat_var_;
   //Glucose::SimpSolver sat_;
   //CMSat::SATSolver sat_;
   LGL* sat_;

   sat_th forward_sat_th_, backward_sat_th_;

   decltype(std::async(std::launch::async, solve_sat_problem, nullptr, sat_vec<sat_literal>{}, nullptr)) sat_handle_;
   Direction cur_sat_reduction_direction_ = Direction::forward;
};



// factors are arranged in trees.
class LP_tree
{
public:
   void AddMessage(MessageTypeAdapter* m, Chirality c) // chirality denotes which factor is upper
   {
      tree_messages_.push_back(std::make_tuple(m, c));
   }
   void compute_subgradient()
   {
      // send messages up the tree
      for(auto it = tree_messages_.begin(); it!= tree_messages_.end(); ++it) {
         auto* m = std::get<0>(*it);
         Chirality c = std::get<1>(*it);
         m->send_message_up(c);
      }
      // compute primal for topmost factor
      // also init primal for top factor, all other primals were initialized already by send_message_up
      if(std::get<1>(tree_messages_.back()) == Chirality::right) {
         // init primal for right factor!
         std::get<0>(tree_messages_.back())->GetRightFactorTypeAdapter()->init_primal();
         std::get<0>(tree_messages_.back())->GetRightFactorTypeAdapter()->MaximizePotentialAndComputePrimal();
      } else {
         std::get<0>(tree_messages_.back())->GetLeftFactorTypeAdapter()->init_primal(); 
         std::get<0>(tree_messages_.back())->GetLeftFactorTypeAdapter()->MaximizePotentialAndComputePrimal(); 
      }
      // track down optimal primal solution
      for(auto it = tree_messages_.rbegin(); it!= tree_messages_.rend(); ++it) {
         auto* m = std::get<0>(*it);
         Chirality c = std::get<1>(*it);
         m->track_solution_down(c);
      } 
   }

   REAL lower_bound() const 
   {
      REAL lb = 0.0;
      for(auto it = tree_messages_.begin(); it!= tree_messages_.end(); ++it) {
         auto* m = std::get<0>(*it);
         Chirality c = std::get<1>(*it);
         assert(false); // why is m and c not used?
         if(std::get<1>(tree_messages_.back()) == Chirality::right) {
            lb += std::get<0>(tree_messages_.back())->GetLeftFactorTypeAdapter()->LowerBound(); 
         } else {
            lb += std::get<0>(tree_messages_.back())->GetRightFactorTypeAdapter()->LowerBound();
         }
      }
      if(std::get<1>(tree_messages_.back()) == Chirality::right) {
         lb += std::get<0>(tree_messages_.back())->GetRightFactorTypeAdapter()->LowerBound();
      } else {
         lb += std::get<0>(tree_messages_.back())->GetLeftFactorTypeAdapter()->LowerBound(); 
      }
      return lb;
      
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
protected:
   std::vector< std::tuple<MessageTypeAdapter*, Chirality>> tree_messages_; // messages forming a tree. Chirality says which side comprises the lower  factor
   // subgradient information = primal solution to tree
   // for sending messages down we need to know to how many lower factors an upper factor is connected and then we need to average messages sent down appropriately.
};


class LP_with_trees : public LP
{

protected:
   std::vector<LP_tree> trees_;
};

inline void LP::Begin()
{
   CalculatePrimalOffsets();

   repamMode_ = LPReparametrizationMode::Undefined;
   assert(f_.size() > 1); // otherwise we need not perform optimization: Just MaximizePotential f_[0]
}

inline void LP::ComputePass()
{
   const auto omega = get_omega();
   assert(forwardUpdateOrdering_.size() == omega.forward.size());
   assert(forwardUpdateOrdering_.size() == omega.backward.size());
   ComputePass(forwardUpdateOrdering_.begin(), forwardUpdateOrdering_.end(), omega.forward.begin());
   ComputePass(forwardUpdateOrdering_.rbegin(), forwardUpdateOrdering_.rend(), omega.backward.begin());
}

template<typename FACTOR_ITERATOR, typename OMEGA_ITERATOR>
void LP::ComputePass(FACTOR_ITERATOR factorIt, const FACTOR_ITERATOR factorItEnd, OMEGA_ITERATOR omegaIt)
{
   for(; factorIt!=factorItEnd; ++factorIt, ++omegaIt) {
      UpdateFactor(*factorIt, *omegaIt);
   }
}

inline void LP::ComputeAnisotropicWeights()
{
   if(!omega_anisotropic_valid_) {
      omega_anisotropic_valid_ = true;
      auto forward = std::async(std::launch::async, [&](){  
            ComputeAnisotropicWeights(forwardOrdering_.begin(), forwardOrdering_.end(), f_sorted_.begin(), f_sorted_.end(), omegaForwardAnisotropic_);
      });
      ComputeAnisotropicWeights(backwardOrdering_.rbegin(), backwardOrdering_.rend(), f_sorted_.rbegin(), f_sorted_.rend(), omegaBackwardAnisotropic_);
      forward.wait();
   }
}

inline void LP::ComputeUniformWeights()
{
   if(!omega_isotropic_valid_) {
      omega_isotropic_valid_ = true;
      auto forward = std::async(std::launch::async, [&](){  
            ComputeUniformWeights(forwardOrdering_.begin(), forwardOrdering_.end(), omegaForwardIsotropic_, 0);
      });
      ComputeUniformWeights(backwardOrdering_.rbegin(), backwardOrdering_.rend(), omegaBackwardIsotropic_, 0);
      forward.wait();
   };
}

inline void LP::ComputeDampedUniformWeights()
{
   if(!omega_isotropic_damped_valid_) {
      omega_isotropic_damped_valid_ = true;
      auto forward = std::async(std::launch::async, [&](){  
            ComputeUniformWeights(forwardOrdering_.begin(), forwardOrdering_.end(), omegaForwardIsotropicDamped_, 1.0);
      });
      ComputeUniformWeights(backwardOrdering_.rbegin(), backwardOrdering_.rend(), omegaBackwardIsotropicDamped_, 1.0);
      forward.wait();
   };
}

// Here we check whether messages constraints are satisfied
template<typename FACTOR_ITERATOR>
bool LP::CheckPrimalConsistency(FACTOR_ITERATOR factorIt, const FACTOR_ITERATOR factorEndIt) const
{
   for(auto msgIt = m_.begin(); msgIt!=m_.end(); ++msgIt) {
      if(!(*msgIt)->CheckPrimalConsistency()) {
         std::cout << "message constraints are not fulfilled by primal solution\n";
         return false;
      }
   }
   return true;
}

template<typename FACTOR_ITERATOR>
REAL LP::EvaluatePrimal(FACTOR_ITERATOR factorIt, const FACTOR_ITERATOR factorEndIt) const
{
   const bool consistent = CheckPrimalConsistency(factorIt, factorEndIt);
   if(consistent == false) return std::numeric_limits<REAL>::infinity();

   REAL cost = 0.0;
   for(; factorIt!=factorEndIt; ++factorIt) {
      cost += (*factorIt)->EvaluatePrimal();
   }
   std::cout << "primal cost = " << cost << "\n";
   return cost;
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

// do zrobienia: possibly templatize this for use with iterators
// note: this function is not working properly. We should only compute factors for messages which can actually send
template<typename FACTOR_ITERATOR, typename FACTOR_SORT_ITERATOR>
void LP::ComputeAnisotropicWeights(
      FACTOR_ITERATOR factorIt, FACTOR_ITERATOR factorEndIt, // sorted pointers to factors
      FACTOR_SORT_ITERATOR factor_sort_begin, FACTOR_SORT_ITERATOR factor_sort_end, // sorted factor indices in f_
      two_dim_variable_array<REAL>& omega)
{

   std::vector<INDEX> f_sorted_inverse(std::distance(factor_sort_begin, factor_sort_end));
   for(INDEX i=0; i<std::distance(factor_sort_begin, factor_sort_end); ++i) {
      f_sorted_inverse[ factor_sort_begin[i] ] = i;
   }

   assert(std::distance(factorIt,factorEndIt) == f_.size());
   assert(std::distance(factor_sort_begin, factor_sort_end) == f_.size());

   
   //std::map<FactorTypeAdapter*, INDEX> factorToIndex;
   //std::map<INDEX, FactorTypeAdapter*> indexToFactor;
   //BuildIndexMaps(factorIt, factorEndIt, factorToIndex, indexToFactor);
   

   // compute the following numbers: 
   // 1) #{factors after current one, to which messages are sent from current factor}
   // 2) #{factors after current one, which receive messages from current one}
//#ifdef LP_MP_PARALLEL
//   std::unique_ptr<std::atomic<INDEX>[]> no_send_factors(new std::atomic<INDEX>[f_.size()]);
//   std::fill(no_send_factors.get(), no_send_factors.get() + f_.size(), 0);
//   std::unique_ptr<std::atomic<INDEX>[]> no_send_factors_later(new std::atomic<INDEX>[f_.size()]);
//   std::fill(no_send_factors_later.get(), no_send_factors_later.get() + f_.size(), 0);
//   std::unique_ptr<std::atomic<INDEX>[]> no_receiving_factors_later(new std::atomic<INDEX>[f_.size()]);
//   std::fill(no_receiving_factors_later.get(), no_receiving_factors_later.get() + f_.size(), 0);
//   std::unique_ptr<std::atomic<INDEX>[]> last_receiving_factor(new std::atomic<INDEX>[f_.size()]); // what is the last (in the order given by factor iterator) factor that receives a message?
//   std::fill(last_receiving_factor.get(), last_receiving_factor.get() + f_.size(), 0);
//#else
   std::vector<INDEX> no_send_factors(f_.size(),0);
   std::vector<INDEX> no_send_factors_later(f_.size(),0);
   std::vector<INDEX> no_receiving_factors_later(f_.size(),0);
   std::vector<INDEX> last_receiving_factor(f_.size(), 0);
//#endif

   // do zrobienia: if factor is not visited at all, then omega is not needed for that entry. We must filter out such entries still
//#pragma omp parallel for schedule(guided)
   for(INDEX i=0; i<m_.size(); ++i) {
      auto* f_left = m_[i]->GetLeftFactor();
      const INDEX f_index_left = factor_address_to_index_[f_left];
      const INDEX index_left = f_sorted_inverse[f_index_left];
      //assert(index_left == factorToIndex[f_left]);
      //const INDEX index_left = factorToIndex[f_left];
      auto* f_right = m_[i]->GetRightFactor();
      const INDEX f_index_right = factor_address_to_index_[f_right];
      const INDEX index_right = f_sorted_inverse[f_index_right];
      //assert(index_right == factorToIndex[f_right]);
      
      if(m_[i]->ReceivesMessageFromLeft()) {
         if(index_left < index_right) {
            no_receiving_factors_later[index_left]++;
         }
//#ifdef LP_MP_PARALLEL
//         INDEX old_val = last_receiving_factor[index_left];
//         const INDEX new_val = std::max(old_val, index_left);
//         while(old_val < new_val && !last_receiving_factor[index_left].compare_exchange_weak(old_val, new_val)) ;
//#else
         last_receiving_factor[index_left] = std::max(last_receiving_factor[index_left], index_right);
//#endif
      }

      if(m_[i]->ReceivesMessageFromRight()) {
         if(index_left > index_right) {
            no_receiving_factors_later[index_right]++;
         }
//#ifdef LP_MP_PARALLEL
//         INDEX old_val = last_receiving_factor[index_right];
//         const INDEX new_val = std::max(old_val, index_right);
//         while(old_val < new_val && !last_receiving_factor[index_right].compare_exchange_weak(old_val, new_val)) ;
//#else
         last_receiving_factor[index_right] = std::max(last_receiving_factor[index_right], index_left);
//#endif
      }
   }

//#pragma omp parallel for schedule(guided)
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
     std::vector<std::vector<FactorTypeAdapter*> > fc = ComputeSendFactorConnection(factorIt,factorEndIt); // possibly not the most efficient way!
     auto fcIt = fc.begin();
     for(auto it=factorIt; it!=factorEndIt; ++it, ++fcIt) {
       if((*it)->FactorUpdated()) {
         omega_size[c] = fcIt->size();
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
       //assert(i == factorToIndex[*it]);
       assert(i == f_sorted_inverse[ factor_address_to_index_[*it] ]);
       if((*it)->FactorUpdated()) {
         INDEX k=0;
         for(auto mIt=(*it)->begin(); mIt!=(*it)->end(); ++mIt) {
           if(mIt.CanSendMessage()) {
             auto* f_connected = mIt.GetConnectedFactor();
             //const INDEX j = factorToIndex[ f_connected ];
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

   // rewrite omega to be a fixed two dimensional array

   //auto last_valid_omega = std::remove_if(omega.begin(), omega.end(), [&](auto& weights) { 
   //      const INDEX i = &weights - &*omega.begin();
   //      return !(*(factorIt+i))->FactorUpdated();
   //      });
   //omega.resize(std::distance(omega.begin(), last_valid_omega)); 
   //omega.shrink_to_fit();

   /*

   assert(std::distance(factorIt,factorEndIt) == f_.size());
   std::vector<std::vector<FactorTypeAdapter*> > fc = ComputeSendFactorConnection(factorIt, factorEndIt);
   std::vector<std::vector<bool> > fcAccessedLater(f_.size());

   auto fcIt = fc.begin();
   auto factorItTmp = factorIt;
   auto fcAccessedLaterIt = fcAccessedLater.begin();
   for(; factorItTmp != factorEndIt; ++factorItTmp, ++fcAccessedLaterIt, ++fcIt) {

      INDEX no_outgong_factors_later = 0;
      INDEX no_outgong_factors_earlier = 0;
      INDEX no_incoming_factors_later = 0;

      assert(factorToIndex.find(*factorItTmp) != factorToIndex.end());
      const INDEX index1 = factorToIndex[*factorItTmp];
      assert(index1 == factorItTmp - factorIt);
      //(*fcAccessedLaterIt).resize(fcIt->size(),false);
      for(INDEX j=0; j<(*fcIt).size(); j++) {
         const INDEX index2 = factorToIndex[ (*fcIt)[j] ];
         if() {
            if(index1 < index2) {
               ++no_outgoing_factors_later;
            } else {
               ++no_outgoing_factors_earlier;
            }
         }
         if( && index1 < index2) {
            ++no_incoming_factors_later;
         }
         assert(factorToIndex.find((*fcIt)[j]) != factorToIndex.end());
         //std:: cout << index2 << ", ";
         bool intermedFactor = false;
         // do zrobienia: this will take extremely long for factors connected to very many other factors.
         
         for(auto fIt2=fc[index2].begin(); fIt2!=fc[index2].end(); ++fIt2) { // note that this is extremely slow for min cost flow factor which is connected to all unary factors!
            assert(factorToIndex.find(*fIt2) != factorToIndex.end());
            const INDEX index3 = factorToIndex[*fIt2];
            if(index1 < index3 && index1 != index3)
               intermedFactor = true;
         }
         
         if(index1 < index2 || intermedFactor == true) {
         //if(index1 < index2 ) {
            (*fcAccessedLaterIt)[j] = true;
         }
      }
   }
   
   omega.clear();
   omega.reserve(factorEndIt - factorIt);
   if(m_.size() == 0) { 
      std::cout << "no messages in problem\n"; 
      //assert(false);
      omega.resize(factorEndIt - factorIt);
      return;
   }

   fcIt = fc.begin();
   fcAccessedLaterIt = fcAccessedLater.begin();
   for(;factorIt != factorEndIt; ++factorIt, ++fcIt, ++fcAccessedLaterIt) {
      if((*factorIt)->FactorUpdated()) {
         assert(fcIt->size() == fcAccessedLaterIt->size());
         omega.push_back(std::vector<REAL>(fcIt->size(), 0.0));
         INDEX noFactorsAccessedLater = 0;
         for(INDEX j=0;j<fcIt->size(); j++) {
            if((*fcAccessedLaterIt)[j]) {
               noFactorsAccessedLater++;
            }
         }
         //const INDEX numberActiveMessages = std::count(fcAccessedLaterIt->begin(), fcAccessedLaterIt->end(), true);

         //const REAL weight = 1.0 / REAL(noFactorsAccessedLater);
         //const REAL weight = 0.8*1.0 / REAL(noFactorsAccessedLater);
         // do zrobienia: not the traditional way
         const REAL weight = 1.0 / (std::max(REAL(noFactorsAccessedLater), REAL(fcAccessedLaterIt->size() - noFactorsAccessedLater)) );
         //const REAL weight = 0.1; // 0.5 works well for pure assignment with equality messages
         for(INDEX j=0; j<fcIt->size(); j++) {
            if((*fcAccessedLaterIt)[j]) {
               omega.back()[j] = weight; 
            } else {
               omega.back()[j] = 0.0;
            }
         }
         assert( std::accumulate(omega.back().begin(), omega.back().end(),0.0) <= 1.0 + eps);
      }
   }
  */
   /*
   for(INDEX i=0; i<omega.size(); ++i) {
      for(INDEX j=0; j<omega[i].size(); ++j) {
         std::cout << omega[i][j] << ", ";
      }
      if(omega[i].size() > 0) {
         std::cout << "\n";
      }
   }
   std::cout << "\n";
   omega.shrink_to_fit();
   */
}

// compute uniform weights so as to help decoding for obtaining primal solutions
// do zrobienia: make omega a TwoDimVariableArray, same in ComputeAnisotropicWeights
// leave_weight signals how much weight to leave in sending factor. Important for rounding and tightening
template<typename FACTOR_ITERATOR>
void LP::ComputeUniformWeights(FACTOR_ITERATOR factorIt, FACTOR_ITERATOR factorEndIt, two_dim_variable_array<REAL>& omega, const REAL leave_weight)
{
   assert(leave_weight >= 0.0 && leave_weight <= 1.0);
   assert(factorEndIt - factorIt == f_.size());
   std::vector<std::vector<FactorTypeAdapter*> > fc = ComputeSendFactorConnection(factorIt,factorEndIt);
   assert(f_.size() == fc.size());

   std::vector<INDEX> omega_size(std::distance(factorIt, factorEndIt),0);
   auto fcIt = fc.begin();
   INDEX c=0;
   for(; factorIt != factorEndIt; ++factorIt, ++fcIt) {
      if((*factorIt)->FactorUpdated()) {
        omega_size[c] = fcIt->size();
        ++c;
      }
   }
   omega_size.resize(c);
   omega = two_dim_variable_array<REAL>(omega_size);
   for(INDEX i=0; i<omega.size(); ++i) {
     for(INDEX j=0; j<omega_size[i]; ++j) {
       omega[i][j] = 1.0/REAL( omega_size[i] + leave_weight );
     }
   }
}

template<typename FACTOR_ITERATOR, typename OMEGA_ITERATOR>
void LP::ComputePassAndPrimal(FACTOR_ITERATOR factorIt, const FACTOR_ITERATOR factorEndIt, OMEGA_ITERATOR omegaIt, INDEX iteration)
{
   //assert(false); // initialize primal before going over it
   for(auto factorItTmp = factorIt; factorItTmp!=factorEndIt; ++factorItTmp, ++omegaIt) {
      UpdateFactorPrimal(*factorItTmp, *omegaIt, iteration);
   }
}

template<typename FACTOR_ITERATOR>
void LP::InitializePrimalVector(FACTOR_ITERATOR factorIt, const FACTOR_ITERATOR factorEndIt, PrimalSolutionStorage& v)
{
   //v.Init(factorIt, factorEndIt);
   v = PrimalSolutionStorage(factorIt, factorEndIt);
}

//inline const PrimalSolutionStorage& LP::GetBestPrimal() const
//{
//   if(bestForwardPrimalCost_ < bestBackwardPrimalCost_) {
//      return bestForwardPrimal_;
//   } else {
//      return bestBackwardPrimal_;
//   }
//}




} // end namespace LP_MP

#endif // LP_MP_MAIN

