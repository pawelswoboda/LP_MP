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

namespace LP_MP {

// forward declaration
class MessageTypeAdapter;
class MessageIterator;

// here we store the primal solution as a std::vector<bool>
// Each factor is allocated factor.size() entries. Offsets are stored in an auxiliary std::vector<INDEX>.
// Do zrobienia: offsets could be stored once for multiple primal solutions. Then iterator based on offsets would be more complicated, though.
class PrimalSolutionStorage {
public:
   // revert to bool
   using Element = std::vector<bool>::iterator;
   using ConstElement = std::vector<bool>::const_iterator;

   PrimalSolutionStorage() {}

   void Initialize()
   {
      std::fill(primal_.begin(), primal_.end(), true);
      /*
      for(INDEX i=0; i<primal_.size(); ++i) {
         primal_[i] = true;
      }
      */
   }
   // we must make a deep copy for the offset variables, as they point directly to entried of primal_
   void operator=(const PrimalSolutionStorage& rhs)
   {
      *this = PrimalSolutionStorage(rhs);
   }
   PrimalSolutionStorage(const PrimalSolutionStorage& rhs)
   {
      primal_ = rhs.primal_;
      offset_.resize(rhs.offset_.size());
      for(INDEX i=0; i<offset_.size(); ++i) {
         offset_[i] = primal_.begin() + (offset_[i] - offset_[0]);
      }
   }
   template<typename FACTOR_ITERATOR>
   PrimalSolutionStorage(FACTOR_ITERATOR factorItBegin, FACTOR_ITERATOR factorItEnd)
   {
      Init(factorItBegin, factorItEnd);
   }
   template<typename FACTOR_ITERATOR>
   void Init(FACTOR_ITERATOR factorItBegin, FACTOR_ITERATOR factorItEnd)  
   {
      primal_.clear();
      offset_.clear();

      auto factorItTmp = factorItBegin;
      INDEX size = 0;
      while(factorItTmp != factorItEnd) {
         size += (*factorItTmp)->size();
         ++factorItTmp;
      }
      primal_.resize(size,true);
      offset_.reserve(factorItEnd - factorItBegin);
      INDEX offset = 0;
      while(factorItBegin != factorItEnd) {
         offset_.push_back(primal_.begin() + offset);
         offset += (*factorItBegin)->size();
         ++factorItBegin;
      }
   }
   INDEX size() const { return offset_.size(); }
   Element operator[](const INDEX i) { assert(i < offset_.size()); return offset_[i]; }
   ConstElement operator[](const INDEX i) const { assert(i < offset_.size()); return offset_[i]; }
   std::vector<Element>::iterator begin() { return offset_.begin(); }
   std::vector<Element>::iterator end() { return offset_.end(); }
   std::vector<Element>::reverse_iterator rbegin() { return offset_.rbegin(); }
   std::vector<Element>::reverse_iterator rend() { return offset_.rend(); }

private:
   std::vector<bool> primal_;
   std::vector<Element> offset_; // order of factors
   std::vector<std::vector<Element>> adjacentOffset_; // adjacent factors for every factor for which ReceiveRestrictedMessages can be called.
};


// pure virtual base class for factors used by LP class
class FactorTypeAdapter
{
public:
   virtual ~FactorTypeAdapter() {};
   virtual void UpdateFactor(const std::vector<REAL>& omega) = 0;
   virtual void UpdateFactor(const std::vector<REAL>& omega, typename PrimalSolutionStorage::Element primalIt) = 0;
   virtual INDEX size() const = 0;
   MessageIterator begin(); 
   MessageIterator end();
   virtual const INDEX GetNoMessages() const = 0;
   virtual MessageTypeAdapter* GetMessage(const INDEX n) const = 0;
   virtual FactorTypeAdapter* GetConnectedFactor(const INDEX i) const = 0;
   virtual bool CanSendMessage(const INDEX i) const = 0;
   virtual REAL LowerBound() const = 0;
   virtual std::vector<REAL> GetReparametrizedPotential() const = 0;
   //virtual PrimalSolutionStorageAdapter* AllocatePrimalSolutionStorage() const = 0;
   virtual bool CanComputePrimalSolution() const = 0;
   virtual void ComputePrimalThroughMessages(typename PrimalSolutionStorage::Element& primalSolution, std::vector<typename PrimalSolutionStorage::Element>& connectedPrimalSolution) const = 0;
   // do zrobienia: this function is not needed. Evaluation can be performed automatically
   virtual REAL EvaluatePrimal(typename PrimalSolutionStorage::Element primalSolution) const = 0;
   // do zrobienia: this is not needed as well and could be automated. Possibly it is good to keep this to enable solution rewriting.
   virtual void WritePrimal(PrimalSolutionStorage::Element primalSolution, std::ofstream& fs) const = 0;
};

class MessageTypeAdapter
{
public:
   virtual ~MessageTypeAdapter() {};
   virtual FactorTypeAdapter* GetLeftFactor() const = 0;
   virtual FactorTypeAdapter* GetRightFactor() const = 0;
   virtual bool CheckPrimalConsistency(typename PrimalSolutionStorage::Element left, typename PrimalSolutionStorage::Element right) const = 0;
   virtual void SetMessage(const std::valarray<REAL>& m) = 0; // do zrobienia: change to vector
   virtual const std::valarray<REAL> GetMessage() const = 0; // do zrobienia: change to vector
   
   // Also true, if SendMessagesTo{Left|Right} is active. Used for weight computation. Disregard message in weight computation if it does not send messages at all
   // do zrobienia: throw them out again
   //virtual bool CanSendMessageToLeft() const = 0;
   //virtual bool CanSendMessageToRight() const = 0;

   // possibly remove these functions again. They are not used anymore
   virtual INDEX GetMessageNumber() const = 0; // give message number as specified in MessageList meta::list
   virtual REAL GetMessageWeightToRight() const = 0;
   virtual REAL GetMessageWeightToLeft() const = 0;

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


// steers optimization of LP solver.
//enum class LPVisitorReturnType {
//   SetRoundingReparametrization, SetAnisotropicReparametrization, Reparametrize, ReparametrizeAndComputePrimal, Break, Error
//};
enum class LPVisitorReturnType {
   ReparametrizeUniform,ReparametrizeLowerBoundUniform,ReparametrizeLowerBoundPrimalUniform,ReparametrizePrimalUniform,
   ReparametrizeAnisotropic,ReparametrizeLowerBoundAnisotropic,ReparametrizeLowerBoundPrimalAnisotropic,ReparametrizePrimalAnisotropic,
   Break,Error
};
enum class LPReparametrizationMode {Anisotropic, Uniform, Undefined};

class LP
{
public:
   LP();
   ~LP();
   INDEX AddFactor(FactorTypeAdapter* f);
   INDEX GetNumberOfFactors() const;
   FactorTypeAdapter* GetFactor(const INDEX i) const;
   INDEX AddMessage(MessageTypeAdapter* m);
   MessageTypeAdapter* GetMessage(const INDEX i) const;
   INDEX GetNumberOfMessages() const;
   void AddFactorRelation(FactorTypeAdapter* f1, FactorTypeAdapter* f2); // indicate that factor f1 comes before factor f2

   void Init(); // must be called after all messages and factors have been added
   void SortFactors();
   template<typename FACTOR_ITERATOR>
      std::vector<std::vector<FactorTypeAdapter*> > ComputeFactorConnection(FACTOR_ITERATOR factorIt, FACTOR_ITERATOR factorEndIt);
   template<typename FACTOR_ITERATOR>
      std::vector<std::vector<FactorTypeAdapter*> > ComputeSendFactorConnection(FACTOR_ITERATOR factorIt, FACTOR_ITERATOR factorEndIt);

   void ComputeAnisotropicWeights();
   template<typename FACTOR_ITERATOR>
   void ComputeAnisotropicWeights(FACTOR_ITERATOR factorIt, FACTOR_ITERATOR factorItEnd, std::vector<std::vector<REAL> >& omega); 
   void ComputeUniformWeights();
   template<typename FACTOR_ITERATOR>
   void ComputeUniformWeights(FACTOR_ITERATOR factorIt, FACTOR_ITERATOR factorItEnd, std::vector<std::vector<REAL> >& omega); // do zrobienia: rename to isotropic weights

   REAL LowerBound() const;
   template<typename FACTOR_ITERATOR, typename PRIMAL_STORAGE_ITERATOR>
      REAL EvaluatePrimal(FACTOR_ITERATOR factorIt, const FACTOR_ITERATOR factorEndIt, PRIMAL_STORAGE_ITERATOR primalIt) const;
   void UpdateFactor(FactorTypeAdapter* f, const std::vector<REAL>& omega); // perform one block coordinate step for factor f
   void ComputePass();
   void ComputePassAndPrimal();
   void ComputeLowerBound();
   template<typename FACTOR_ITERATOR, typename OMEGA_ITERATOR>
      void ComputePass(FACTOR_ITERATOR factorIt, const FACTOR_ITERATOR factorItEnd, OMEGA_ITERATOR omegaIt);
   void UpdateFactor(FactorTypeAdapter* f, const std::vector<REAL>& omega, typename PrimalSolutionStorage::Element primal);
   template<typename FACTOR_ITERATOR, typename OMEGA_ITERATOR, typename PRIMAL_SOLUTION_STORAGE_ITERATOR>
   void ComputePassAndPrimal(FACTOR_ITERATOR factorIt, const FACTOR_ITERATOR factorItEnd, OMEGA_ITERATOR omegaIt, 
         PRIMAL_SOLUTION_STORAGE_ITERATOR primalIt, PRIMAL_SOLUTION_STORAGE_ITERATOR bestPrimalIt, REAL& bestPrimalcost);
   template<typename VISITOR>
   SIGNED_INDEX Solve(VISITOR& v);

   std::vector<std::vector<REAL> > GetReparametrizedModel() const;
   std::vector<std::vector<REAL> > GetReparametrizedModel(const INDEX begin, const INDEX end) const;

   //PrimalSolutionStorage GetBestPrimal() const;
   REAL BestPrimalBound() const { return std::min(bestForwardPrimalCost_, bestBackwardPrimalCost_); }
   REAL BestLowerBound() const { return bestLowerBound_; }
   REAL CurrentLowerBound() const { return currentLowerBound_; }
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
   std::vector<std::vector<REAL> > omegaForward_, omegaBackward_;
   std::vector<std::pair<FactorTypeAdapter*, FactorTypeAdapter*> > factorRel_; // factor ordering relations. First factor must come before second factor. factorRel_ must describe a DAG

   template<typename FACTOR_ITERATOR>
      void InitializePrimalVector(FACTOR_ITERATOR factorIt, const FACTOR_ITERATOR factorEndIt, PrimalSolutionStorage& v);
   template<typename FACTOR_ITERATOR, typename PRIMAL_STORAGE_ITERATOR>
      void ComputePrimalThroughMessages(FACTOR_ITERATOR factorIt, const FACTOR_ITERATOR factorEndIt, PRIMAL_STORAGE_ITERATOR primalIt);
   template<typename FACTOR_ITERATOR, typename PRIMAL_ITERATOR>
      bool CheckPrimalConsistency(FACTOR_ITERATOR factorIt, const FACTOR_ITERATOR factorEndIt, PRIMAL_ITERATOR primalIt) const;

   REAL bestForwardPrimalCost_ = std::numeric_limits<REAL>::max(), 
        bestBackwardPrimalCost_ = std::numeric_limits<REAL>::max();

   REAL bestLowerBound_ = -std::numeric_limits<REAL>::max();
   REAL currentLowerBound_ = -std::numeric_limits<REAL>::max();

   PrimalSolutionStorage forwardPrimal_, backwardPrimal_, bestForwardPrimal_, bestBackwardPrimal_; // note: these vectors are stored in the order of forwardOrdering_

   LPReparametrizationMode repamMode_ = LPReparametrizationMode::Undefined;
};

inline void LP::Init()
{
   std::cout << "Determining factor ordering." << std::endl;
   SortFactors();

   // initialize three arrays of primal solutions corresponding to factors, one computed in the forward pass, one computed in the backward pass, and one with the best solution obtained so far
   InitializePrimalVector(forwardOrdering_.begin(), forwardOrdering_.end(), forwardPrimal_);
   InitializePrimalVector(forwardOrdering_.begin(), forwardOrdering_.end(), bestForwardPrimal_);
   InitializePrimalVector(forwardOrdering_.begin(), forwardOrdering_.end(), backwardPrimal_);
   InitializePrimalVector(forwardOrdering_.begin(), forwardOrdering_.end(), bestBackwardPrimal_);
   //InitializePrimalVector(backwardOrdering_.begin(), backwardOrdering_.end(), backwardPrimal_);
   //InitializePrimalVector(backwardOrdering_.begin(), backwardOrdering_.end(), bestBackwardPrimal_);

   repamMode_ = LPReparametrizationMode::Undefined;
   /*
   if(repamMode_ == LPReparametrizationMode::Anisotropic) {
      ComputeAnisotropicWeights(forwardOrdering_, omegaForward_);
      ComputeAnisotropicWeights(backwardOrdering_, omegaBackward_);
   } else if(repamMode_ == LPReparametrizationMode::Rounding) {
      ComputeUniformWeights(forwardOrdering_, omegaForward_);
      ComputeUniformWeights(backwardOrdering_, omegaBackward_);
   }
   */
   assert(f_.size() > 1); // otherwise we need not perform optimization: Just MaximizePotential f_[0]
}

inline void LP::ComputePass()
{
   ComputePass(forwardOrdering_.begin(), forwardOrdering_.end(), omegaForward_.begin());
   ComputePass(forwardOrdering_.rbegin(), forwardOrdering_.rend(), omegaBackward_.begin());
}

template<typename FACTOR_ITERATOR, typename OMEGA_ITERATOR>
void LP::ComputePass(FACTOR_ITERATOR factorIt, const FACTOR_ITERATOR factorItEnd, OMEGA_ITERATOR omegaIt)
{
   for(; factorIt!=factorItEnd; ++factorIt, ++omegaIt) {
      UpdateFactor(*factorIt, *omegaIt);
   }
}

inline void LP::ComputePassAndPrimal()
{
   forwardPrimal_.Initialize();
   assert(forwardPrimal_.size() == forwardOrdering_.size());
   ComputePassAndPrimal(forwardOrdering_.begin(), forwardOrdering_.end(), omegaForward_.begin(), forwardPrimal_.begin(), bestForwardPrimal_.begin(), bestForwardPrimalCost_); 
   backwardPrimal_.Initialize();
   assert(backwardPrimal_.size() == forwardOrdering_.size());
   ComputePassAndPrimal(forwardOrdering_.rbegin(), forwardOrdering_.rend(), omegaBackward_.begin(), backwardPrimal_.rbegin(), bestBackwardPrimal_.rbegin(), bestBackwardPrimalCost_); 
}

inline void LP::ComputeLowerBound()
{
   currentLowerBound_ = LowerBound();
   bestLowerBound_ = std::max(currentLowerBound_,bestLowerBound_);
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
   }
}

template<typename FACTOR_ITERATOR, typename OMEGA_ITERATOR, typename PRIMAL_SOLUTION_STORAGE_ITERATOR>
void LP::ComputePassAndPrimal(FACTOR_ITERATOR factorIt, const FACTOR_ITERATOR factorEndIt, OMEGA_ITERATOR omegaIt, 
      PRIMAL_SOLUTION_STORAGE_ITERATOR primalIt, PRIMAL_SOLUTION_STORAGE_ITERATOR bestPrimalIt, REAL& bestPrimalCost)
{
   // first we need to set primal storage iterator to zero
   //ResetPrimalStorage(primalIt, primalIt + INDEX(factorEndIt-factorIt));
   PRIMAL_SOLUTION_STORAGE_ITERATOR primalItTmp = primalIt;
   FACTOR_ITERATOR factorItTmp = factorIt;
   for(; factorItTmp!=factorEndIt; ++factorItTmp, ++omegaIt, ++primalItTmp) {
      UpdateFactor(*factorItTmp, *omegaIt, *primalItTmp);
   }
   ComputePrimalThroughMessages(factorIt, factorEndIt, primalIt);
   const REAL currentPrimalCost = EvaluatePrimal(factorIt,factorEndIt,primalIt);
   if(currentPrimalCost < bestPrimalCost || bestPrimalCost >= std::numeric_limits<REAL>::max()) { // the second case occurs whenever we have some infeasible primal solution. Possibly, the newer solution, although infeasible, are not as infeasible as the old one.
      bestPrimalCost = currentPrimalCost;
      std::swap(*(primalIt),*(bestPrimalIt));
   }
}


template<typename VISITOR>
SIGNED_INDEX LP::Solve(VISITOR& v)
{
/*
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
*/
   // note: currently init is called before solve. If repeated solve is done, then we do not to call init, unless factors or messages have changed
   Init();

   LPVisitorReturnType s = v.begin(this);

   while(true) {
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
            ComputePassAndPrimal();
            ComputeLowerBound();
            s = v.template visit<LPVisitorReturnType::ReparametrizeLowerBoundPrimalUniform>(this);
            break;
         case LPVisitorReturnType::ReparametrizePrimalUniform:
            ComputeUniformWeights();
            ComputePassAndPrimal();
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
            ComputePassAndPrimal();
            ComputeLowerBound();
            s = v.template visit<LPVisitorReturnType::ReparametrizeLowerBoundPrimalAnisotropic>(this);
            break;
         case LPVisitorReturnType::ReparametrizePrimalAnisotropic:
            ComputeAnisotropicWeights();
            ComputePassAndPrimal();
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




            /*
         case LPVisitorReturnType::SetAnisotropicReparametrization:
            if(repamMode_ != LPReparametrizationMode::Anisotropic) {
               ComputeAnisotropicWeights(forwardOrdering_.begin(), forwardOrdering_.end(), omegaForward_);
               ComputeAnisotropicWeights(forwardOrdering_.rbegin(), forwardOrdering_.rend(), omegaBackward_);
               repamMode_ = LPReparametrizationMode::Anisotropic;
            }
            s = v.template visit<LPVisitorReturnType::SetAnisotropicReparametrization>(this);
            break;
         case LPVisitorReturnType::SetRoundingReparametrization:
            if(repamMode_ != LPReparametrizationMode::Rounding) {
               ComputeUniformWeights(forwardOrdering_.begin(), forwardOrdering_.end(), omegaForward_);
               ComputeUniformWeights(forwardOrdering_.rbegin(), forwardOrdering_.rend(), omegaBackward_);
               repamMode_ = LPReparametrizationMode::Rounding;
            }
            s = v.template visit<LPVisitorReturnType::SetRoundingReparametrization>(this);
            break;
         case LPVisitorReturnType::Reparametrize:
            //std::cout << "reparametrize in " << (repamMode_ == LPReparametrizationMode::Rounding ? " rounding " : " anisotropic ") << "mode\n";
            ComputePass(forwardOrdering_.begin(), forwardOrdering_.end(), omegaForward_.begin());
            ComputePass(forwardOrdering_.rbegin(), forwardOrdering_.rend(), omegaBackward_.begin());
            s = v.template visit<LPVisitorReturnType::Reparametrize>(this);
            break;
         case LPVisitorReturnType::ReparametrizeAndComputePrimal:
            //std::cout << "reparametrize in " << (repamMode_ == LPReparametrizationMode::Rounding ? " rounding " : " anisotropic ") << "mode\n";
            forwardPrimal_.Initialize();
            assert(forwardPrimal_.size() == forwardOrdering_.size());
            ComputePassAndPrimal(forwardOrdering_.begin(), forwardOrdering_.end(), omegaForward_.begin(), forwardPrimal_.begin(), bestForwardPrimal_.begin(), bestForwardPrimalCost_); 
            backwardPrimal_.Initialize();
            assert(backwardPrimal_.size() == forwardOrdering_.size());
            ComputePassAndPrimal(forwardOrdering_.rbegin(), forwardOrdering_.rend(), omegaBackward_.begin(), backwardPrimal_.rbegin(), bestBackwardPrimal_.rbegin(), bestBackwardPrimalCost_); 
            s = v.template visit<LPVisitorReturnType::ReparametrizeAndComputePrimal>(this);
            break;
         case LPVisitorReturnType::Break:
            s = v.template visit<LPVisitorReturnType::Break>(this);
            return 0;
            break;
         case LPVisitorReturnType::Error:
            return -1;
            break;
            */
      }
   }
}

// do zrobienia: maybe const possible for argument v
// note: this function propagates primals in bfs-manner beginning from factors, for which primal was computed in the update pass.
// Especially, the propagation in the factor message network must enable this and ensure, that no primal conflict can occur
template<typename FACTOR_ITERATOR, typename PRIMAL_STORAGE_ITERATOR>
void LP::ComputePrimalThroughMessages(FACTOR_ITERATOR factorIt, const FACTOR_ITERATOR factorEndIt, PRIMAL_STORAGE_ITERATOR primalIt)
{
   std::map<FactorTypeAdapter*, INDEX> factorToIndex;
   std::map<INDEX, FactorTypeAdapter*> indexToFactor; // do zrobienia: this should equal forwardFactors_. Modify BuildIndexMaps
   BuildIndexMaps(factorIt,factorEndIt,factorToIndex,indexToFactor);

   std::vector<bool> visited(factorEndIt - factorIt, false);

   std::queue<INDEX> q;
   for(INDEX i=0; factorIt+i!=factorEndIt; ++i) {
      const FactorTypeAdapter* f = *(factorIt+i);
      if(f->CanComputePrimalSolution()) {
         q.push(i); // we assume that the corresponding primal solutions in v have been computed
      }
   }
   assert(q.size() > 0); // otherwise primal solution cannot be propagated

   while(!q.empty()) {
      const INDEX i = q.front();
      q.pop();
      if(visited[i] == true) { continue; }
      else { visited[i] = true; }

      FactorTypeAdapter* const f = *(factorIt+i);
      // build vector of PrimalSolutionStorageAdapter* such that the factors referenced by the current factors messages are in there
      std::vector<typename PrimalSolutionStorage::Element> cur_primal_adjacent;
      cur_primal_adjacent.reserve(f->GetNoMessages());
      for(MessageIterator it=f->begin(); it!=f->end(); ++it) {
         const INDEX j = factorToIndex[it.GetConnectedFactor()]; // do zrobienia: is this the right dereferencing operator definition? Possibly change in LP_MP.h
         cur_primal_adjacent.push_back(*(primalIt+j));
         q.push(j);
      }
      f->ComputePrimalThroughMessages(*(primalIt+i), cur_primal_adjacent);
   }

   assert(std::count(visited.begin(), visited.end(),false) == 0); // primal solutions have been computed for every factor
}

// Here we check whether messages do not
template<typename FACTOR_ITERATOR, typename PRIMAL_ITERATOR>
bool LP::CheckPrimalConsistency(FACTOR_ITERATOR factorIt, const FACTOR_ITERATOR factorEndIt, PRIMAL_ITERATOR primalIt) const
{
   std::map<FactorTypeAdapter*, INDEX> factorToIndex;
   std::map<INDEX, FactorTypeAdapter*> indexToFactor; // do zrobienia: this should equal f_. Modify BuildIndexMaps
   BuildIndexMaps(factorIt,factorEndIt,factorToIndex,indexToFactor); // do zrobienia: precompute BuildIndexMaps by hashing it with pair(factorIt,factorEndIt). Also this can be done faster by noting that factorIt is either forward or backward order of primalIt, hence it is just simple index calculations everywhere

   for(auto msgIt = m_.begin(); msgIt!=m_.end(); ++msgIt) {
      FactorTypeAdapter* const leftFactor = (*msgIt)->GetLeftFactor();
      FactorTypeAdapter* const rightFactor = (*msgIt)->GetRightFactor();
      const INDEX leftFactorIndex = factorToIndex[leftFactor];
      const INDEX rightFactorIndex = factorToIndex[rightFactor];
      assert(leftFactorIndex < f_.size() && rightFactorIndex < f_.size());
      const bool consistent = (*msgIt)->CheckPrimalConsistency(*(primalIt+leftFactorIndex), *(primalIt+rightFactorIndex));
      if(consistent == false) return false;
   }

   return true;
}

template<typename FACTOR_ITERATOR, typename PRIMAL_STORAGE_ITERATOR>
REAL LP::EvaluatePrimal(FACTOR_ITERATOR factorIt, const FACTOR_ITERATOR factorEndIt, PRIMAL_STORAGE_ITERATOR primalIt) const
{
   const bool consistent = CheckPrimalConsistency(factorIt, factorEndIt, primalIt);
   if(consistent == false) return std::numeric_limits<REAL>::max();

   REAL cost = 0.0;
   for(; factorIt!=factorEndIt; ++factorIt, ++primalIt) {
      cost += (*factorIt)->EvaluatePrimal(*primalIt);
   }
   return cost;
}

template<typename FACTOR_ITERATOR>
void LP::InitializePrimalVector(FACTOR_ITERATOR factorIt, const FACTOR_ITERATOR factorEndIt, PrimalSolutionStorage& v)
{
   v.Init(factorIt, factorEndIt);
}

// return best primal solution in the order of f_. Note: this is not a deep copy, but just a vector of pointers to solutions stored in the LP class
/*
inline PrimalSolutionStorage LP::GetBestPrimal() const
{
   const bool forward = bestForwardPrimalCost_ < bestBackwardPrimalCost_;
   //const std::vector<PrimalSolutionStorageAdapter*>& bestPrimal = forward ? bestForwardPrimal_ : bestBackwardPrimal_;
   //const std::vector<FactorTypeAdapter*>& f = forward ? forwardOrdering_ : backwardOrdering_;

   assert(false);
   PrimalSolutionStorage primalSorted(f_.begin(),f_.end());

   //auto perm = ComputePermutation(f_,f);
   //for(INDEX i=0; i<perm.size(); ++i) {
   //   primalSorted[i] = bestPrimal[perm[i]];
   //}

   return primalSorted;
}
*/

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
   omega.resize(factorEndIt - factorIt);
   if(m_.size() == 0) { 
      std::cout << "no messages in problem\n"; 
      assert(false);
      return;
   }

   fcIt = fc.begin();
   fcAccessedLaterIt = fcAccessedLater.begin();
   for(auto omegaIt = omega.begin(); factorIt != factorEndIt; ++factorIt, ++omegaIt, ++fcIt, ++fcAccessedLaterIt) {
      assert(fcIt->size() == fcAccessedLaterIt->size());
      (*omegaIt).resize(fcIt->size(), 0.0);
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
            (*omegaIt)[j] = weight; 
         } else {
            (*omegaIt)[j] = 0.0;
         }
      }
      assert( std::accumulate(omegaIt->begin(), omegaIt->end(),0.0) <= 1.0 + eps);
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
   omega.resize(factorEndIt-factorIt);
   auto omegaIt = omega.begin();
   auto fcIt = fc.begin();
   for(; factorIt != factorEndIt; ++factorIt, ++omegaIt, ++fcIt) {
      // better for rounding
      (*omegaIt) = std::vector<REAL>(fcIt->size(), 1.0/REAL(fcIt->size() + 0.1) );
      // better for dual convergence
      //(*omegaIt) = std::vector<REAL>(fcIt->size(), 1.0/REAL(fcIt->size()) );
   }
}


} // end namespace LP_MP

#endif // LP_MP_MAIN

