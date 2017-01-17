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

namespace LP_MP {

// forward declaration
class MessageTypeAdapter;
class MessageIterator;

// pure virtual base class for factor container used by LP class
class FactorTypeAdapter
{
public:
   virtual ~FactorTypeAdapter() {}
   virtual void UpdateFactor(const std::vector<REAL>& omega) = 0;
   virtual void UpdateFactorPrimal(const std::vector<REAL>& omega, const INDEX iteration) = 0;
   virtual bool FactorUpdated() const = 0; // does calling UpdateFactor do anything? If no, it need not be called while in ComputePass, saving time.
   virtual INDEX size() const = 0;
   virtual INDEX PrimalSize() const = 0;
   MessageIterator begin(); 
   MessageIterator end();
   virtual const INDEX GetNoMessages() const = 0;
   virtual MessageTypeAdapter* GetMessage(const INDEX n) const = 0;
   virtual FactorTypeAdapter* GetConnectedFactor(const INDEX i) const = 0;
   virtual bool CanSendMessage(const INDEX i) const = 0;
   virtual REAL LowerBound() const = 0;
   virtual std::vector<REAL> GetReparametrizedPotential() const = 0;
   //virtual PrimalSolutionStorageAdapter* AllocatePrimalSolutionStorage() const = 0;
   //virtual bool CanComputePrimalSolution() const = 0;
   // the offset in the primal storage
   virtual void SetPrimalOffset(const INDEX) = 0;
   virtual INDEX GetPrimalOffset() const = 0;
   
   virtual void SetAuxOffset(const INDEX n) = 0;
   virtual INDEX GetAuxOffset() const = 0;
   
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
   virtual FactorTypeAdapter* GetLeftFactor() const = 0;
   virtual FactorTypeAdapter* GetRightFactor() const = 0;
   //virtual bool CheckPrimalConsistency(typename PrimalSolutionStorage::Element left, typename PrimalSolutionStorage::Element right) const = 0;
   virtual bool SendsMessageToLeft() const = 0;
   virtual bool SendsMessageToRight() const = 0;
   virtual bool ReceivesMessageFromLeft() const = 0;
   virtual bool ReceivesMessageFromRight() const = 0;
   virtual bool CheckPrimalConsistency() const = 0;

   //virtual void send_message_up(FactorTypeAdapter* lower, FactorTypeAdapter* upper) = 0;
   
   // Also true, if SendMessagesTo{Left|Right} is active. Used for weight computation. Disregard message in weight computation if it does not send messages at all
   // do zrobienia: throw them out again
   //virtual bool CanSendMessageToLeft() const = 0;
   //virtual bool CanSendMessageToRight() const = 0;

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
inline MessageIterator FactorTypeAdapter::end()  { return MessageIterator(this, GetNoMessages()); };


class LP {
public:
   LP() {}
   ~LP() 
   {
      for(INDEX i=0; i<m_.size(); i++) { delete m_[i]; }
      for(INDEX i=0; i<f_.size(); i++) { delete f_[i]; }
   }
   // make a deep copy of factors and messages. Adjust pointers of messages
   LP(const LP& o) {
      f_.reserve(o.f_.size());
      assert(false);
      for(auto& f : o.f_) {

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
   }

   INDEX AddFactor(FactorTypeAdapter* f)
   {
      //INDEX primalOffset;
      //if(f_.size() ==0) {
      //   primalOffset = 0;
      //} else {
      //   primalOffset = f_.back()->GetPrimalOffset() + f_.back()->PrimalSize();
      //}
      //std::cout << "primal offset = " << primalOffset << "\n";
      //f->SetPrimalOffset(primalOffset);
      f_.push_back(f);
      return f_.size() - 1;

   }
   INDEX GetNumberOfFactors() const { return f_.size(); }
   FactorTypeAdapter* GetFactor(const INDEX i) const { return f_[i]; }
   INDEX size() const { INDEX size=0; for(auto* f : f_) { size += f->size(); } return size; }
   INDEX AddMessage(MessageTypeAdapter* m)
   {
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

   void ForwardPassFactorRelation(FactorTypeAdapter* f1, FactorTypeAdapter* f2) { assert(f1!=f2); forward_pass_factor_rel_.push_back({f1,f2}); }
   void BackwardPassFactorRelation(FactorTypeAdapter* f1, FactorTypeAdapter* f2) { assert(f1!=f2); backward_pass_factor_rel_.push_back({f1,f2}); }


   void Init(); // must be called after all messages and factors have been added
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

      std::map<FactorTypeAdapter*,INDEX> factorToIndex; // possibly do it with a hash_map for speed
      std::map<INDEX,FactorTypeAdapter*> indexToFactor; // do zrobienia: need not be map, oculd be vector!
      BuildIndexMaps(f_.begin(), f_.end(),factorToIndex,indexToFactor);

      for(auto fRelIt=factor_rel.begin(); fRelIt!=factor_rel.end(); fRelIt++) {
         // do zrobienia: why do these asserts fail?
         assert(factorToIndex.find(fRelIt->first ) != factorToIndex.end());
         assert(factorToIndex.find(fRelIt->second) != factorToIndex.end());
         INDEX f1 = factorToIndex[fRelIt->first];
         INDEX f2 = factorToIndex[fRelIt->second];
         g.addEdge(f1,f2);
      }

      std::vector<INDEX> sortedIndices = g.topologicalSort();
      assert(sortedIndices.size() == f_.size());

      std::vector<FactorTypeAdapter*> fSorted;
      fSorted.reserve(f_.size());
      for(INDEX i=0; i<sortedIndices.size(); i++) {
         fSorted.push_back( indexToFactor[sortedIndices[i]] );
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
      std::map<FactorTypeAdapter*, INDEX> factorToIndexSorted;
      std::map<INDEX, FactorTypeAdapter*> indexToFactorSorted;
      BuildIndexMaps(ordering.begin(), ordering.end(), factorToIndexSorted, indexToFactorSorted);
      for(auto rel : factor_rel) {
         const INDEX index_left = factorToIndexSorted[ std::get<0>(rel) ];
         const INDEX index_right = factorToIndexSorted[ std::get<1>(rel) ];
         assert(index_left < index_right);
      }
   }
   void SortFactors()
   {
      SortFactors(forward_pass_factor_rel_, forwardOrdering_, forwardUpdateOrdering_);
      SortFactors(backward_pass_factor_rel_, backwardOrdering_, backwardUpdateOrdering_);
   }

   template<typename FACTOR_ITERATOR>
      std::vector<std::vector<FactorTypeAdapter*> > ComputeFactorConnection(FACTOR_ITERATOR factorIt, FACTOR_ITERATOR factorEndIt);
   template<typename FACTOR_ITERATOR>
      std::vector<std::vector<FactorTypeAdapter*> > ComputeSendFactorConnection(FACTOR_ITERATOR factorIt, FACTOR_ITERATOR factorEndIt);

   void ComputeWeights(const LPReparametrizationMode m);
   void ComputeAnisotropicWeights();
   template<typename FACTOR_ITERATOR>
   void ComputeAnisotropicWeights(FACTOR_ITERATOR factorIt, FACTOR_ITERATOR factorItEnd, std::vector<std::vector<REAL> >& omega); 
   void ComputeUniformWeights();
   template<typename FACTOR_ITERATOR>
   void ComputeUniformWeights(FACTOR_ITERATOR factorIt, FACTOR_ITERATOR factorEndIt, std::vector<std::vector<REAL> >& omega, const REAL leave_weight = 0.01); // do zrobienia: rename to isotropic weights

   REAL LowerBound() const
   {
      REAL lb = 0.0;
      for(auto fIt=f_.begin(); fIt!=f_.end(); fIt++) {
         lb += (*fIt)->LowerBound();
         assert( (*fIt)->LowerBound() > -10000000.0);
      }
      return lb;
   }

   void InitializePrimalVector(PrimalSolutionStorage& p) { InitializePrimalVector(f_.begin(), f_.end(), p); }
   template<typename FACTOR_ITERATOR>
      void InitializePrimalVector(FACTOR_ITERATOR factorIt, const FACTOR_ITERATOR factorEndIt, PrimalSolutionStorage& v);
   template<typename FACTOR_ITERATOR>
      bool CheckPrimalConsistency(FACTOR_ITERATOR factorIt, const FACTOR_ITERATOR factorEndIt) const;

   REAL EvaluatePrimal() {
      return EvaluatePrimal(f_.begin(), f_.end());
   }
   template<typename FACTOR_ITERATOR>
   REAL EvaluatePrimal(FACTOR_ITERATOR factorIt, const FACTOR_ITERATOR factorEndIt) const;

   void UpdateFactor(FactorTypeAdapter* f, const std::vector<REAL>& omega) // perform one block coordinate step for factor f
   {
      f->UpdateFactor(omega);
   }
   void ComputePass();

   void ComputePassAndPrimal(const INDEX iteration);

   template<typename FACTOR_ITERATOR, typename OMEGA_ITERATOR>
   void ComputePassAndPrimal(FACTOR_ITERATOR factorIt, const FACTOR_ITERATOR factorEndIt, OMEGA_ITERATOR omegaIt, const INDEX iteration);

   template<typename FACTOR_ITERATOR, typename OMEGA_ITERATOR>
   void ComputePass(FACTOR_ITERATOR factorIt, const FACTOR_ITERATOR factorItEnd, OMEGA_ITERATOR omegaIt);
   void UpdateFactorPrimal(FactorTypeAdapter* f, const std::vector<REAL>& omega, const INDEX iteration)
   {
      f->UpdateFactorPrimal(omega, iteration);
   }

   const PrimalSolutionStorage& GetBestPrimal() const;
   //REAL BestPrimalBound() const { return std::min(bestForwardPrimalCost_, bestBackwardPrimalCost_); }
   //REAL BestLowerBound() const { return bestLowerBound_; }
   //REAL CurrentLowerBound() const { return currentLowerBound_; }
   template<typename FACTOR_ITERATOR>
      void WritePrimal(FACTOR_ITERATOR factorIt, const FACTOR_ITERATOR factorEndIt, std::ofstream& fs) const;
      void WritePrimal(const INDEX factorIndexBegin, const INDEX factorIndexEnd, std::ofstream& fs) const;

   LPReparametrizationMode GetRepamMode() const { return repamMode_; }
private:
   // do zrobienia: possibly hold factors and messages in shared_ptr?
   std::vector<FactorTypeAdapter*> f_; // note that here the factors are stored in the original order they were given. They will be output in this order as well, e.g. by problemDecomposition
   std::vector<MessageTypeAdapter*> m_;
   std::vector<FactorTypeAdapter*> forwardOrdering_, backwardOrdering_; // separate forward and backward ordering are not needed: Just store factorOrdering_ and generate forward order by begin() and backward order by rbegin().
   std::vector<FactorTypeAdapter*> forwardUpdateOrdering_, backwardUpdateOrdering_; // like forwardOrdering_, but includes only those factors where UpdateFactor actually does something
   std::vector<std::vector<REAL> > omegaForward_, omegaBackward_;
   std::vector<std::pair<FactorTypeAdapter*, FactorTypeAdapter*> > forward_pass_factor_rel_, backward_pass_factor_rel_; // factor ordering relations. First factor must come before second factor. factorRel_ must describe a DAG

   
   //REAL bestLowerBound_ = -std::numeric_limits<REAL>::infinity();
   //REAL currentLowerBound_ = -std::numeric_limits<REAL>::infinity();

   LPReparametrizationMode repamMode_ = LPReparametrizationMode::Undefined;
};

// factors are arranged in trees.
/*
class LP_tree
{
public:
   void AddMessage(MessageTypeAdapter* m, FactorTypeAdapter* lower_factor, FactorTypeAdapter* upper_factor) 
   {
      if(lower_factor == m->GetLeftFactor()) {
         tree_messages_.push_back(std::make_tuple(m, Chirality::left));
      } else {
         tree_messages_.push_back(std::make_tuple(m, Chirality::right));
      }
   }
   void compute_tree_subgradient()
   {
      // reparametrize everything into upmost factor
      for(auto it = tree_messages_.begin(); it!= tree_messages_.end(); ++it) {
         auto* m = std::get<0>(*it);
         Chirality c = std::get<1>(*it);
         m->send_message_up(c);
      }
      // compute primal for topmost factor
      if(std::get<1>(tree_messages_.back()) == Chirality::left) {
         tree_messages_.back()->RightFactor()->MaximizePotentialAndComputePrimal();
      } else {
         tree_messages_.back()->LeftFactor()->MaximizePotentialAndComputePrimal(); 
      }
      // track down optimal primal solution
      for(auto it = tree_messages_.rbegin(); it!= tree_messages_.rend(); ++it) {
         auto* m = std::get<0>(*it);
         Chirality c = std::get<1>(*it);
         m->trace_solution_down_tree(c);
      } 
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
*/

inline void LP::Init()
{
   //std::cout << "Determining factor ordering." << std::endl;
   SortFactors();
   CalculatePrimalOffsets();
   // recalculation of primal offsets is needed whenever factors have been changed (factors can change size!) or have been added.

   // initialize three arrays of primal solutions corresponding to factors, one computed in the forward pass, one computed in the backward pass, and one with the best solution obtained so far
   //InitializePrimalVector(forwardOrdering_.begin(), forwardOrdering_.end(), forwardPrimal_);
   //InitializePrimalVector(forwardOrdering_.begin(), forwardOrdering_.end(), bestForwardPrimal_);
   //InitializePrimalVector(forwardOrdering_.begin(), forwardOrdering_.end(), backwardPrimal_);
   //InitializePrimalVector(forwardOrdering_.begin(), forwardOrdering_.end(), bestBackwardPrimal_);
   //InitializePrimalVector(backwardOrdering_.begin(), backwardOrdering_.end(), backwardPrimal_);
   //InitializePrimalVector(backwardOrdering_.begin(), backwardOrdering_.end(), bestBackwardPrimal_);

   repamMode_ = LPReparametrizationMode::Undefined;
   assert(f_.size() > 1); // otherwise we need not perform optimization: Just MaximizePotential f_[0]
}

inline void LP::ComputePass()
{
   assert(forwardUpdateOrdering_.size() == omegaForward_.size());
   assert(forwardUpdateOrdering_.size() == omegaBackward_.size());
   ComputePass(forwardUpdateOrdering_.begin(), forwardUpdateOrdering_.end(), omegaForward_.begin());
   ComputePass(forwardUpdateOrdering_.rbegin(), forwardUpdateOrdering_.rend(), omegaBackward_.begin());
}

template<typename FACTOR_ITERATOR, typename OMEGA_ITERATOR>
void LP::ComputePass(FACTOR_ITERATOR factorIt, const FACTOR_ITERATOR factorItEnd, OMEGA_ITERATOR omegaIt)
{
   for(; factorIt!=factorItEnd; ++factorIt, ++omegaIt) {
      UpdateFactor(*factorIt, *omegaIt);
   }
}

//inline void LP::ComputeLowerBound()
//{
//   currentLowerBound_ = LowerBound();
//   bestLowerBound_ = std::max(currentLowerBound_,bestLowerBound_);
//}

inline void LP::ComputeWeights(const LPReparametrizationMode m)
{
   assert(m != LPReparametrizationMode::Undefined);
   if(repamMode_ != m) {
      if(m == LPReparametrizationMode::Anisotropic) {
         ComputeAnisotropicWeights();
      } else if(m == LPReparametrizationMode::Uniform) {
         ComputeUniformWeights();
      } else {
         throw std::runtime_error("unknown repam mode");
      }
   }
}
inline void LP::ComputeAnisotropicWeights()
{
   if(repamMode_ != LPReparametrizationMode::Anisotropic) {
      ComputeAnisotropicWeights(forwardOrdering_.begin(), forwardOrdering_.end(), omegaForward_);
      ComputeAnisotropicWeights(backwardOrdering_.rbegin(), backwardOrdering_.rend(), omegaBackward_);
      repamMode_ = LPReparametrizationMode::Anisotropic;
   }
}

inline void LP::ComputeUniformWeights()
{
   if(repamMode_ != LPReparametrizationMode::Uniform) {
      ComputeUniformWeights(forwardOrdering_.begin(), forwardOrdering_.end(), omegaForward_);
      ComputeUniformWeights(backwardOrdering_.rbegin(), backwardOrdering_.rend(), omegaBackward_);
      //assert(omegaForward_.size() == forwardOrdering_.size()); // need not be true for factors that are not sending/receiving
      //assert(omegaBackward_.size() == forwardOrdering_.size());
      repamMode_ = LPReparametrizationMode::Uniform;
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
   return cost;
}


// write primal solutions in bounds [factorIndexBegin,factorIndexEnd) to filestream
template<typename FACTOR_ITERATOR>
void LP::WritePrimal(FACTOR_ITERATOR factorIt, const FACTOR_ITERATOR factorEndIt, std::ofstream& fs) const
{
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
template<typename FACTOR_ITERATOR>
void LP::ComputeAnisotropicWeights(FACTOR_ITERATOR factorIt, FACTOR_ITERATOR factorEndIt, std::vector<std::vector<REAL> >& omega)
{
   assert(std::distance(factorIt,factorEndIt) == f_.size());
   std::map<FactorTypeAdapter*, INDEX> factorToIndex;
   std::map<INDEX, FactorTypeAdapter*> indexToFactor;
   BuildIndexMaps(factorIt, factorEndIt, factorToIndex, indexToFactor);

   // compute the following numbers: 
   // 1) #{factors after current one, to which messages are sent from current factor}
   // 2) #{factors after current one, which receive messages from current one}
   std::vector<INDEX> no_send_factors(f_.size());
   std::vector<INDEX> no_send_factors_later(f_.size(),0);
   std::vector<INDEX> no_receiving_factors_later(f_.size(),0);

   // what is the last (in the order given by factor iterator) factor that receives a message?
   std::vector<INDEX> last_receiving_factor(f_.size(), 0);

   // do zrobienia: if factor is not visited at all, then omega is not needed for that entry. We must filter out such entries still
   for(auto* m : m_) {
      auto* f_left = m->GetLeftFactor();
      const INDEX index_left = factorToIndex[f_left];
      auto* f_right = m->GetRightFactor();
      const INDEX index_right = factorToIndex[f_right];
      
      if(m->ReceivesMessageFromLeft()) {
         if(index_left < index_right) {
            no_receiving_factors_later[index_left]++;
         }
         last_receiving_factor[index_left] = std::max(last_receiving_factor[index_left], index_right);
      }

      if(m->ReceivesMessageFromRight()) {
         if(index_left > index_right) {
            no_receiving_factors_later[index_right]++;
         }
         last_receiving_factor[index_right] = std::max(last_receiving_factor[index_right], index_left);
      }
   }

   for(auto* m : m_) {
      auto* f_left = m->GetLeftFactor();
      const INDEX index_left = factorToIndex[f_left];
      auto* f_right = m->GetRightFactor();
      const INDEX index_right = factorToIndex[f_right];

      if(m->SendsMessageToRight()) {
         no_send_factors[index_left]++;
         if(index_left < index_right || last_receiving_factor[index_right] > index_left) {
            no_send_factors_later[index_left]++;
         }
      }
      if(m->SendsMessageToLeft()) {
         no_send_factors[index_right]++;
         if(index_right < index_left || last_receiving_factor[index_left] > index_right) {
            no_send_factors_later[index_right]++;
         }
      }
   }


   omega.clear();
   omega.resize(std::distance(factorIt, factorEndIt),std::vector<REAL>(0));
   for(auto it=factorIt; it!=factorEndIt; ++it) {
      const INDEX i = factorToIndex[*it];
      for(auto mIt=(*it)->begin(); mIt!=(*it)->end(); ++mIt) {
         if(mIt.CanSendMessage()) {
            auto* f_connected = mIt.GetConnectedFactor();
            const INDEX j = factorToIndex[ f_connected ];
            if(i<j || last_receiving_factor[j] > i) {
               omega[i].push_back(1.0/REAL(no_receiving_factors_later[i] + std::max(no_send_factors_later[i], no_send_factors[i] - no_send_factors_later[i])));
            } else {
               omega[i].push_back(0.0);
            } 
         }
      }
   }


   /*
   for(auto* m : m_) {
      auto* f_left = m->GetLeftFactor();
      const INDEX index_left = factorToIndex[f_left];
      auto* f_right = m->GetRightFactor();
      const INDEX index_right = factorToIndex[f_right];
      assert(index_left != index_right);

      if(m->SendsMessageToRight()) {
         if(index_left < index_right) {
            omega[index_left].push_back(1.0/REAL(no_receiving_factors_later[index_left] + std::max(no_send_factors_later[index_left], no_send_factors[index_left] - no_send_factors_later[index_left])));
         } else {
            omega[index_left].push_back(0.0);
         }
      }
      if(m->SendsMessageToLeft()) {
         if(index_right < index_left) {
            omega[index_right].push_back(1.0/REAL(no_receiving_factors_later[index_right] + std::max(no_send_factors_later[index_right], no_send_factors[index_right] - no_send_factors_later[index_right])));
         } else {
            omega[index_right].push_back(0.0);
         }
      }
   }
   */

   // check whether all messages were added to m_. Possibly, this can be automated: Traverse all factors, get all messages, add them to m_ and avoid duplicates along the way.
   assert(2*m_.size() == std::accumulate(f_.begin(), f_.end(), 0, [](INDEX sum, auto* f){ return sum + f->GetNoMessages(); }));
   assert(HasUniqueValues(m_));
   for(INDEX i=0; i<f_.size(); ++i) {
      assert(omega[i].size() <= (*(factorIt+i))->GetNoMessages());
      assert(std::accumulate(omega[i].begin(), omega[i].end(), 0.0) <= 1.0 + eps);
   }

   std::vector<std::vector<REAL>> omega_filtered;
   for(INDEX i=0; i<omega.size(); ++i) {
      if( (*(factorIt+i))->FactorUpdated() ) {
         omega_filtered.push_back(omega[i]);
      }
   }
   std::swap(omega,omega_filtered);
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
void LP::ComputeUniformWeights(FACTOR_ITERATOR factorIt, FACTOR_ITERATOR factorEndIt, std::vector<std::vector<REAL> >& omega, const REAL leave_weight)
{
   assert(leave_weight >= 0.0 && leave_weight <= 1.0);
   assert(factorEndIt - factorIt == f_.size());
   std::vector<std::vector<FactorTypeAdapter*> > fc = ComputeSendFactorConnection(factorIt,factorEndIt);
   assert(f_.size() == fc.size());

   omega.clear();
   omega.reserve(factorEndIt-factorIt);
   auto omegaIt = omega.begin();
   auto fcIt = fc.begin();
   for(; factorIt != factorEndIt; ++factorIt, ++fcIt) {
      // do zrobienia: let the factor below be variable and specified on the command line
      // better for rounding
      if((*factorIt)->FactorUpdated()) {
         omega.push_back(std::vector<REAL>(fcIt->size(), 1.0/REAL(2.0*fcIt->size() + leave_weight) ));
      }
      //(*omegaIt) = std::vector<REAL>(fcIt->size(), 1.0/REAL(fcIt->size() + 1.0) );
      // better for dual convergence
      //(*omegaIt) = std::vector<REAL>(fcIt->size(), 1.0/REAL(fcIt->size()) );
   }
   omega.shrink_to_fit();
}

inline void LP::ComputePassAndPrimal(const INDEX iteration)
{
   ComputePassAndPrimal(forwardUpdateOrdering_.begin(), forwardUpdateOrdering_.end(), omegaForward_.begin(), 2*iteration+1); // timestamp must be > 0, otherwise in the first iteration primal does not get initialized
   const REAL forward_cost = EvaluatePrimal();
   std::cout << "forward cost = " << forward_cost << "\n";
   
   ComputePassAndPrimal(forwardUpdateOrdering_.rbegin(), forwardUpdateOrdering_.rend(), omegaBackward_.begin(), 2*iteration + 2); 
   const REAL backward_cost = EvaluatePrimal();
   std::cout << "backward cost = " << backward_cost << "\n";
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

