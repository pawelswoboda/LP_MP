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
#include "lp_interface/lp_interface.h"
#include "primal_solution_storage.hxx"

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
   virtual void UpdateFactor(const std::vector<REAL>& omega, typename PrimalSolutionStorage::Element primalIt) = 0;
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
   // do zrobienia: this function is not needed. Evaluation can be performed automatically
   virtual REAL EvaluatePrimal(typename PrimalSolutionStorage::Element primalSolution) const = 0;
   // do zrobienia: this is not needed as well and could be automated. Possibly it is good to keep this to enable solution rewriting.
   //virtual void WritePrimal(PrimalSolutionStorage::Element primalSolution, std::ofstream& fs) const = 0;

   // for the LP interface
   virtual void CreateConstraints(LpInterfaceAdapter* lpInterface) const = 0;
};

class MessageTypeAdapter
{
public:
   virtual ~MessageTypeAdapter() {}
   virtual FactorTypeAdapter* GetLeftFactor() const = 0;
   virtual FactorTypeAdapter* GetRightFactor() const = 0;
   //virtual bool CheckPrimalConsistency(typename PrimalSolutionStorage::Element left, typename PrimalSolutionStorage::Element right) const = 0;
   virtual bool CheckPrimalConsistency(typename PrimalSolutionStorage::Element primal) const = 0;
   virtual void SetMessage(const std::valarray<REAL>& m) = 0; // do zrobienia: function is not used anywhere currently // do zrobienia: change to vector
   virtual const std::valarray<REAL> GetMessage() const = 0; // do zrobienia: function is not used anywhere currently // do zrobienia: change to vector
   
   // Also true, if SendMessagesTo{Left|Right} is active. Used for weight computation. Disregard message in weight computation if it does not send messages at all
   // do zrobienia: throw them out again
   //virtual bool CanSendMessageToLeft() const = 0;
   //virtual bool CanSendMessageToRight() const = 0;

   // possibly remove these functions again. They are not used anymore
   virtual INDEX GetMessageNumber() const = 0; // give message number as specified in MessageList meta::list
   //virtual REAL GetMessageWeightToRight() const = 0;
   //virtual REAL GetMessageWeightToLeft() const = 0;

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


class LP
{
public:
   LP() {}
   ~LP() 
   {
      for(INDEX i=0; i<m_.size(); i++) { delete m_[i]; }
      for(INDEX i=0; i<f_.size(); i++) { delete f_[i]; }
   }
   INDEX AddFactor(FactorTypeAdapter* f)
   {
      INDEX primalOffset;
      if(f_.size() ==0) {
         primalOffset = 0;
      } else {
         primalOffset = f_.back()->GetPrimalOffset() + f_.back()->PrimalSize();
      }
      //std::cout << "primal offset = " << primalOffset << "\n";
      f_.push_back(f);
      f->SetPrimalOffset(primalOffset);
      return f_.size() - 1;

   }
   INDEX GetNumberOfFactors() const { return f_.size(); }
   FactorTypeAdapter* GetFactor(const INDEX i) const { return f_[i]; }
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
      factorRel_.push_back(std::make_pair(f1,f2));
   }


   void Init(); // must be called after all messages and factors have been added
   void SortFactors()
   {
      // assume that factorRel_ describe a DAG. Compute topological sorting
      Topological_Sort::Graph g(f_.size());

      std::map<FactorTypeAdapter*,INDEX> factorToIndex; // possibly do it with a hash_map for speed
      std::map<INDEX,FactorTypeAdapter*> indexToFactor; // do zrobienia: need not be map, oculd be vector!
      BuildIndexMaps(f_.begin(), f_.end(),factorToIndex,indexToFactor);

      for(auto fRelIt=factorRel_.begin(); fRelIt!=factorRel_.end(); fRelIt++) {
         // do zrobienia: why do these asserts fail?
         //assert(factorToIndex.find(fRelIt->first ) != factorToIndex.end());
         //assert(factorToIndex.find(fRelIt->second) != factorToIndex.end());
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
      forwardOrdering_  = fSorted;
      forwardUpdateOrdering_.clear();
      for(auto f : forwardOrdering_) {
         if(f->FactorUpdated()) {
            forwardUpdateOrdering_.push_back(f);
         }
      }
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
   void ComputeUniformWeights(FACTOR_ITERATOR factorIt, FACTOR_ITERATOR factorItEnd, std::vector<std::vector<REAL> >& omega); // do zrobienia: rename to isotropic weights

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
   template<typename FACTOR_ITERATOR, typename PRIMAL_ITERATOR>
      bool CheckPrimalConsistency(FACTOR_ITERATOR factorIt, const FACTOR_ITERATOR factorEndIt, PRIMAL_ITERATOR primalIt) const;

   REAL EvaluatePrimal(PrimalSolutionStorage::Element primal) {
      return EvaluatePrimal(f_.begin(), f_.end(), primal);
   }
   /*
   REAL EvaluatePrimal(PrimalSolutionStorage::Element primal, const REAL primalBound  = std::numeric_limits<REAL>::infinity())
   {
      return EvaluatePrimal(primal, [](PrimalSolutionStorage::Element primal) { return true; }, primalBound);
   }
   */
   template<typename FACTOR_ITERATOR, typename PRIMAL_STORAGE_ITERATOR>
      REAL EvaluatePrimal(FACTOR_ITERATOR factorIt, const FACTOR_ITERATOR factorEndIt, PRIMAL_STORAGE_ITERATOR primalIt) const;
   void UpdateFactor(FactorTypeAdapter* f, const std::vector<REAL>& omega) // perform one block coordinate step for factor f
   {
      f->UpdateFactor(omega);
   }
   void ComputePass();

   void ComputePassAndPrimal(LPReparametrizationMode repam, PrimalSolutionStorage& forwardPrimal, PrimalSolutionStorage& backwardPrimal);

   template<typename FACTOR_ITERATOR, typename OMEGA_ITERATOR, typename PRIMAL_SOLUTION_STORAGE_ITERATOR>
   void ComputePassAndPrimal(FACTOR_ITERATOR factorIt, const FACTOR_ITERATOR factorEndIt, OMEGA_ITERATOR omegaIt, PRIMAL_SOLUTION_STORAGE_ITERATOR primalIt);

   template<typename FACTOR_ITERATOR, typename OMEGA_ITERATOR>
   void ComputePass(FACTOR_ITERATOR factorIt, const FACTOR_ITERATOR factorItEnd, OMEGA_ITERATOR omegaIt);
   void UpdateFactor(FactorTypeAdapter* f, const std::vector<REAL>& omega, typename PrimalSolutionStorage::Element primal)
   {
      f->UpdateFactor(omega, primal);
   }

   const PrimalSolutionStorage& GetBestPrimal() const;
   //REAL BestPrimalBound() const { return std::min(bestForwardPrimalCost_, bestBackwardPrimalCost_); }
   //REAL BestLowerBound() const { return bestLowerBound_; }
   //REAL CurrentLowerBound() const { return currentLowerBound_; }
   template<typename FACTOR_ITERATOR, typename PRIMAL_ITERATOR>
      void WritePrimal(FACTOR_ITERATOR factorIt, const FACTOR_ITERATOR factorEndIt, PRIMAL_ITERATOR primalIt, std::ofstream& fs) const;
   template<typename PRIMAL_ITERATOR>
      void WritePrimal(const INDEX factorIndexBegin, const INDEX factorIndexEnd, PRIMAL_ITERATOR primalIt, std::ofstream& fs) const;

   LPReparametrizationMode GetRepamMode() const { return repamMode_; }
private:
   // do zrobienia: possibly hold factors and messages in shared_ptr?
   std::vector<FactorTypeAdapter*> f_; // note that here the factors are stored in the original order they were given. They will be output in this order as well, e.g. by problemDecomposition
   std::vector<MessageTypeAdapter*> m_;
   std::vector<FactorTypeAdapter*> forwardOrdering_;//, backwardOrdering_; // separate forward and backward ordering are not needed: Just store factorOrdering_ and generate forward order by begin() and backward order by rbegin().
   std::vector<FactorTypeAdapter*> forwardUpdateOrdering_; // like forwardOrdering_, but includes only those factors where UpdateFactor actually does something
   std::vector<std::vector<REAL> > omegaForward_, omegaBackward_;
   std::vector<std::pair<FactorTypeAdapter*, FactorTypeAdapter*> > factorRel_; // factor ordering relations. First factor must come before second factor. factorRel_ must describe a DAG

   
   //REAL bestLowerBound_ = -std::numeric_limits<REAL>::infinity();
   //REAL currentLowerBound_ = -std::numeric_limits<REAL>::infinity();

   LPReparametrizationMode repamMode_ = LPReparametrizationMode::Undefined;
};

inline void LP::Init()
{
   //std::cout << "Determining factor ordering." << std::endl;
   SortFactors();

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
      ComputeAnisotropicWeights(forwardOrdering_.rbegin(), forwardOrdering_.rend(), omegaBackward_);
      repamMode_ = LPReparametrizationMode::Anisotropic;
   }
}

inline void LP::ComputeUniformWeights()
{
   if(repamMode_ != LPReparametrizationMode::Uniform) {
      ComputeUniformWeights(forwardOrdering_.begin(), forwardOrdering_.end(), omegaForward_);
      ComputeUniformWeights(forwardOrdering_.rbegin(), forwardOrdering_.rend(), omegaBackward_);
      repamMode_ = LPReparametrizationMode::Uniform;
   };
}

/*
template<typename VISITOR, typename PRIMAL_CONSISTENCY_CHECK_FCT>
SIGNED_INDEX LP::Solve(VISITOR& v, PRIMAL_CONSISTENCY_CHECK_FCT primalCheck)
{
   std::cout << "Current SIMD implementation: ";
   switch( Vc::CurrentImplementation::current()) {
      case Vc::Implementation::ScalarImpl:
         std::cout << "scalar"; break;
      case Vc::Implementation::SSE2Impl:
            std::cout << "SSE2"; break;
      case Vc::Implementation::SSSE3Impl:
            std::cout << "SSSE3"; break;
      case Vc::Implementation::SSE41Impl:
            std::cout << "SSE41"; break;
      case Vc::Implementation::SSE42Impl:
            std::cout << "SSE42"; break;
      case Vc::Implementation::AVXImpl:
            std::cout << "AVX"; break;
      case Vc::Implementation::AVX2Impl:
            std::cout << "AVX2"; break;
      default:
            std::cout << "unknown";
   }
   std::cout << "\n";
   std::cout << "SIMD vector size = " << REAL_SIMD::Size << "\n";
   // note: currently init is called before solve. If repeated solve is done, then we do not to call init, unless factors or messages have changed
   Init();

   LPVisitorReturnType s = v.begin(this);

   while(true) {
      //std::cout << "repam mode = ";
      //if(repamMode_ == LPReparametrizationMode::Anisotropic) {
      //   std::cout << " anisotropic";
      //} else if(repamMode_ == LPReparametrizationMode::Uniform) {
      //   std::cout << " uniform";
      //} else {
      //   std::cout << "undefined";
      //}
      //std::cout << "\n";
      
      switch(s) {
         case LPVisitorReturnType::ReparametrizeUniform:
            ComputeUniformWeights();
            ComputePass();
            s = v.template visit<LPVisitorReturnType::ReparametrizeUniform>(this);
            break;
         case LPVisitorReturnType::ReparametrizeLowerBoundUniform:
            ComputeUniformWeights();
            ComputePass();
            ComputeLowerBound();
            s = v.template visit<LPVisitorReturnType::ReparametrizeLowerBoundUniform>(this);
            break;
         case LPVisitorReturnType::ReparametrizeLowerBoundPrimalUniform:
            ComputeUniformWeights();
            ComputePassAndPrimal(primalCheck);
            ComputeLowerBound();
            s = v.template visit<LPVisitorReturnType::ReparametrizeLowerBoundPrimalUniform>(this);
            break;
         case LPVisitorReturnType::ReparametrizePrimalUniform:
            ComputeUniformWeights();
            ComputePassAndPrimal(primalCheck);
            s = v.template visit<LPVisitorReturnType::ReparametrizePrimalUniform>(this);
            break;

         case LPVisitorReturnType::ReparametrizeAnisotropic:
            ComputeAnisotropicWeights();
            ComputePass();
            s = v.template visit<LPVisitorReturnType::ReparametrizeAnisotropic>(this);
            break;
         case LPVisitorReturnType::ReparametrizeLowerBoundAnisotropic:
            ComputeAnisotropicWeights();
            ComputePass();
            ComputeLowerBound();
            s = v.template visit<LPVisitorReturnType::ReparametrizeLowerBoundAnisotropic>(this);
            break;
         case LPVisitorReturnType::ReparametrizeLowerBoundPrimalAnisotropic:
            ComputeAnisotropicWeights();
            ComputePassAndPrimal(primalCheck);
            ComputeLowerBound();
            s = v.template visit<LPVisitorReturnType::ReparametrizeLowerBoundPrimalAnisotropic>(this);
            break;
         case LPVisitorReturnType::ReparametrizePrimalAnisotropic:
            ComputeAnisotropicWeights();
            ComputePassAndPrimal(primalCheck);
            s = v.template visit<LPVisitorReturnType::ReparametrizePrimalAnisotropic>(this);
            break;

         case LPVisitorReturnType::Break:
            s = v.template visit<LPVisitorReturnType::Break>(this);
            return 0;
            break;
         case LPVisitorReturnType::Error:
            s = v.template visit<LPVisitorReturnType::Error>(this);
            return -1;
            break;
      }
   }
}
*/

// Here we check whether messages constraints are satisfied
template<typename FACTOR_ITERATOR, typename PRIMAL_ITERATOR>
bool LP::CheckPrimalConsistency(FACTOR_ITERATOR factorIt, const FACTOR_ITERATOR factorEndIt, PRIMAL_ITERATOR primalIt) const
{
   for(auto msgIt = m_.begin(); msgIt!=m_.end(); ++msgIt) {
      if(!(*msgIt)->CheckPrimalConsistency(primalIt)) {
         std::cout << "message constraints are not fulfilled by primal solution\n";
         return false;
      }
   }
   return true;
}

template<typename FACTOR_ITERATOR, typename PRIMAL_STORAGE_ITERATOR>
REAL LP::EvaluatePrimal(FACTOR_ITERATOR factorIt, const FACTOR_ITERATOR factorEndIt, PRIMAL_STORAGE_ITERATOR primalIt) const
{
   const bool consistent = CheckPrimalConsistency(factorIt, factorEndIt, primalIt);
   if(consistent == false) return std::numeric_limits<REAL>::infinity();

   REAL cost = 0.0;
   for(; factorIt!=factorEndIt; ++factorIt) {
      cost += (*factorIt)->EvaluatePrimal(primalIt);
   }
   return cost;
}


// write primal solutions in bounds [factorIndexBegin,factorIndexEnd) to filestream
template<typename FACTOR_ITERATOR, typename PRIMAL_ITERATOR>
void LP::WritePrimal(FACTOR_ITERATOR factorIt, const FACTOR_ITERATOR factorEndIt, PRIMAL_ITERATOR primalIt, std::ofstream& fs) const
{
   for(; factorIt!=factorEndIt; ++factorIt, ++primalIt) {
      (*factorIt)->WritePrimal(*primalIt,fs);
   }
}
// write primal solutions in bounds [factorIndexBegin,factorIndexEnd) to filestream, factor order is given by f_
template<typename PRIMAL_ITERATOR>
void LP::WritePrimal(const INDEX factorIndexBegin, const INDEX factorIndexEnd, PRIMAL_ITERATOR primalIt, std::ofstream& fs) const
{
   WritePrimal(f_.begin() + factorIndexBegin, f_.begin() + factorIndexEnd, primalIt, fs);
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
   assert(factorEndIt - factorIt == f_.size());
   std::vector<std::vector<FactorTypeAdapter*> > fc = ComputeSendFactorConnection(factorIt, factorEndIt);
   std::vector<std::vector<bool> > fcAccessedLater(f_.size());

   std::map<FactorTypeAdapter*, INDEX> factorToIndex;
   std::map<INDEX, FactorTypeAdapter*> indexToFactor;
   BuildIndexMaps(factorIt, factorEndIt, factorToIndex, indexToFactor);

   //for(INDEX i=0; i<f.size(); i++) {
   auto fcIt = fc.begin();
   auto factorItTmp = factorIt;
   auto fcAccessedLaterIt = fcAccessedLater.begin();
   for(; factorItTmp != factorEndIt; ++factorItTmp, ++fcAccessedLaterIt, ++fcIt) {
      assert(factorToIndex.find(*factorItTmp) != factorToIndex.end());
      const INDEX index1 = factorToIndex[*factorItTmp];
      assert(index1 == factorItTmp - factorIt);
      (*fcAccessedLaterIt).resize(fcIt->size(),false);
      for(INDEX j=0; j<(*fcIt).size(); j++) {
         assert(factorToIndex.find((*fcIt)[j]) != factorToIndex.end());
         const INDEX index2 = factorToIndex[ (*fcIt)[j] ];
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
   for(INDEX i=0; i<omega.size(); ++i) {
      for(INDEX j=0; j<omega[i].size(); ++j) {
         //std::cout << omega[i][j] << ", ";
      }
      if(omega[i].size() > 0) {
         //std::cout << "\n";
      }
   }
   //std::cout << "\n";
   omega.shrink_to_fit();
}

// compute uniform weights so as to help decoding for obtaining primal solutions
// do zrobienia: make omega a TwoDimVariableArray, same in ComputeAnisotropicWeights
template<typename FACTOR_ITERATOR>
void LP::ComputeUniformWeights(FACTOR_ITERATOR factorIt, FACTOR_ITERATOR factorEndIt, std::vector<std::vector<REAL> >& omega)
{
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
         omega.push_back(std::vector<REAL>(fcIt->size(), 1.0/REAL(fcIt->size() + 0.001) ));
      }
      //(*omegaIt) = std::vector<REAL>(fcIt->size(), 1.0/REAL(fcIt->size() + 1.0) );
      // better for dual convergence
      //(*omegaIt) = std::vector<REAL>(fcIt->size(), 1.0/REAL(fcIt->size()) );
   }
   omega.shrink_to_fit();
}

inline void LP::ComputePassAndPrimal(LPReparametrizationMode repam, PrimalSolutionStorage& forwardPrimal, PrimalSolutionStorage& backwardPrimal)
{
   if(repam == LPReparametrizationMode::Anisotropic) {
      ComputeAnisotropicWeights();
   } else if(repam == LPReparametrizationMode::Uniform) {
      ComputeUniformWeights();
   } else {
      throw std::runtime_error("repam mode not recognized");
   }
   forwardPrimal.Initialize();
   ComputePassAndPrimal(forwardUpdateOrdering_.begin(), forwardUpdateOrdering_.end(), omegaForward_.begin(), forwardPrimal.begin()); 
   
   backwardPrimal.Initialize();
   ComputePassAndPrimal(forwardUpdateOrdering_.rbegin(), forwardUpdateOrdering_.rend(), omegaBackward_.begin(), backwardPrimal.begin()); 
}

template<typename FACTOR_ITERATOR, typename OMEGA_ITERATOR, typename PRIMAL_SOLUTION_STORAGE_ITERATOR>
void LP::ComputePassAndPrimal(FACTOR_ITERATOR factorIt, const FACTOR_ITERATOR factorEndIt, OMEGA_ITERATOR omegaIt, PRIMAL_SOLUTION_STORAGE_ITERATOR primalIt)
{
   for(auto factorItTmp = factorIt; factorItTmp!=factorEndIt; ++factorItTmp, ++omegaIt) {
      UpdateFactor(*factorItTmp, *omegaIt, primalIt);
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

