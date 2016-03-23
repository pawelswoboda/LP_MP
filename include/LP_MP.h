#ifndef LP_MP_MAIN
#define LP_MP_MAIN

#include "config.hxx"
//#include "instances.inc"
#include <vector>
#include <valarray>
#include <map>
#include <list>
#include <iostream>
#include <numeric>
#include <algorithm>
#include <functional>
#include <utility>
#include <limits>
#include <exception>
#include "template_utilities.hxx"
#include <assert.h>
//#include "help_functions.hxx"
#include "topological_sort.hxx"
#include <memory>

// do zrobienia: make those functions const which are logically so

namespace LP_MP {

// forward declaration
class MessageTypeAdapter;
class MessageIterator;
class PrimalSolutionStorageAdapter;

// virtual base class for factors used by LP class
class FactorTypeAdapter
{
public:
   virtual ~FactorTypeAdapter() {};
   virtual void UpdateFactor(const std::vector<REAL>& omega) = 0;
   virtual void UpdateFactor(const std::vector<REAL>& omega, PrimalSolutionStorageAdapter* const primal) = 0;
   MessageIterator begin(); 
   MessageIterator end();
   virtual const INDEX GetNoMessages() const = 0;
   virtual MessageTypeAdapter* GetMessage(const INDEX n) const = 0;
   virtual FactorTypeAdapter* GetConnectedFactor(const INDEX i) const = 0;
   virtual bool CanSendMessage(const INDEX i) const = 0;
   virtual REAL LowerBound() const = 0;
   virtual std::vector<REAL> GetReparametrizedPotential() const = 0;
   virtual PrimalSolutionStorageAdapter* AllocatePrimalSolutionStorage() const = 0;
   virtual bool CanComputePrimalSolution() const = 0;
   virtual void ComputePrimalThroughMessages(PrimalSolutionStorageAdapter* primalSolution, std::vector<PrimalSolutionStorageAdapter*>& connectedPrimalSolution) const = 0;
   virtual REAL EvaluatePrimal(PrimalSolutionStorageAdapter* primalSolution) const = 0;
   virtual void WritePrimal(PrimalSolutionStorageAdapter* primalSolution, std::ofstream& fs) const = 0;
};

class MessageTypeAdapter
{
public:
   virtual ~MessageTypeAdapter() {};
   virtual FactorTypeAdapter* GetLeftFactor() const = 0;
   virtual FactorTypeAdapter* GetRightFactor() const = 0;
   virtual bool CheckPrimalConsistency(PrimalSolutionStorageAdapter* left, PrimalSolutionStorageAdapter* right) const = 0;
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

// abstract base class for PrimalSolutionStorage classes
class PrimalSolutionStorageAdapter
{
public:
   virtual void reset() = 0;
   virtual ~PrimalSolutionStorageAdapter() {}
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
enum class LPVisitorReturnType {
   SetRoundingReparametrization, SetAnisotropicReparametrization, Reparametrize, ReparametrizeAndComputePrimal, Break, Error
};
enum class LPReparametrizationMode {Anisotropic, Rounding, Undefined};

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
   std::vector<std::vector<FactorTypeAdapter*> > ComputeFactorConnection(const std::vector<FactorTypeAdapter* >& f);
   std::vector<std::vector<FactorTypeAdapter*> > ComputeSendFactorConnection(const std::vector<FactorTypeAdapter* >& f);

   void ComputeAnisotropicWeights(const std::vector<FactorTypeAdapter*>& f, std::vector<std::vector<REAL> >& omega); 
   void ComputeUniformWeights(const std::vector<FactorTypeAdapter*>& f, std::vector<std::vector<REAL> >& omega); // do zrobienia: rename to isotropic weights

   REAL LowerBound() const;
   template<typename FACTOR_ITERATOR, typename PRIMAL_STORAGE_ITERATOR>
      REAL EvaluatePrimal(FACTOR_ITERATOR factorIt, const FACTOR_ITERATOR factorEndIt, PRIMAL_STORAGE_ITERATOR primalIt) const;
   REAL UpdateFactor(FactorTypeAdapter* f, const std::vector<REAL>& omega); // perform one block coordinate step for factor f
   template<typename FACTOR_ITERATOR, typename OMEGA_ITERATOR>
      REAL ComputePass(FACTOR_ITERATOR factorIt, const FACTOR_ITERATOR factorItEnd, OMEGA_ITERATOR omegaIt);
   REAL UpdateFactor(FactorTypeAdapter* f, const std::vector<REAL>& omega, PrimalSolutionStorageAdapter* primal);
   template<typename FACTOR_ITERATOR, typename OMEGA_ITERATOR, typename PRIMAL_SOLUTION_STORAGE_ITERATOR>
   std::pair<REAL,REAL> ComputePassAndPrimal(FACTOR_ITERATOR factorIt, const FACTOR_ITERATOR factorItEnd, OMEGA_ITERATOR omegaIt, 
         PRIMAL_SOLUTION_STORAGE_ITERATOR primalIt, PRIMAL_SOLUTION_STORAGE_ITERATOR bestPrimalIt, REAL& bestPrimalcost);
   template<typename VISITOR>
   SIGNED_INDEX Solve(VISITOR& v);

   std::vector<std::vector<REAL> > GetReparametrizedModel() const;
   std::vector<std::vector<REAL> > GetReparametrizedModel(const INDEX begin, const INDEX end) const;

   std::vector<PrimalSolutionStorageAdapter*> GetBestPrimal() const;
   REAL BestPrimalBound() const { return std::min(bestForwardPrimalCost_, bestBackwardPrimalCost_); }
   template<typename FACTOR_ITERATOR, typename PRIMAL_ITERATOR>
      void WritePrimal(FACTOR_ITERATOR factorIt, const FACTOR_ITERATOR factorEndIt, PRIMAL_ITERATOR primalIt, std::ofstream& fs) const;
   template<typename PRIMAL_ITERATOR>
      void WritePrimal(const INDEX factorIndexBegin, const INDEX factorIndexEnd, PRIMAL_ITERATOR primalIt, std::ofstream& fs) const;

   LPReparametrizationMode GetRepamMode() const { return repamMode_; }
private:
   // do zrobienia: possibly hold factors and messages in shared_ptr?
   std::vector<FactorTypeAdapter*> f_; // note that here the factors are stored in the original order they were given. They will be output in this order as well, e.g. by problemDecomposition
   std::vector<MessageTypeAdapter*> m_;
   std::vector<FactorTypeAdapter*> forwardOrdering_, backwardOrdering_; // separate forward and backward ordering are not needed: Just store factorOrdering_ and generate forward order by begin() and backward order by rbegin().
   std::vector<std::vector<REAL> > omegaForward_, omegaBackward_;
   std::vector<std::pair<FactorTypeAdapter*, FactorTypeAdapter*> > factorRel_; // factor ordering relations. First factor must come before second factor. factorRel_ must describe a DAG

   template<typename FACTOR_ITERATOR>
      void InitializePrimalVector(FACTOR_ITERATOR factorIt, const FACTOR_ITERATOR factorEndIt, std::vector<PrimalSolutionStorageAdapter*>& v);
   template<typename FACTOR_ITERATOR, typename PRIMAL_STORAGE_ITERATOR>
      void ComputePrimalThroughMessages(FACTOR_ITERATOR factorIt, const FACTOR_ITERATOR factorEndIt, PRIMAL_STORAGE_ITERATOR primalIt);
   template<typename FACTOR_ITERATOR, typename PRIMAL_ITERATOR>
      bool CheckPrimalConsistency(FACTOR_ITERATOR factorIt, const FACTOR_ITERATOR factorEndIt, PRIMAL_ITERATOR primalIt) const;

   std::vector<PrimalSolutionStorageAdapter*> forwardPrimal_, backwardPrimal_, bestForwardPrimal_, bestBackwardPrimal_; // note: these vectors are stored in the order of forwardOrdering_ and backwardOrdering_ respectively
   REAL bestForwardPrimalCost_ = std::numeric_limits<REAL>::max(), 
        bestBackwardPrimalCost_ = std::numeric_limits<REAL>::max();

   LPReparametrizationMode repamMode_ = LPReparametrizationMode::Undefined;
};

inline void LP::Init()
{
   std::cout << "Determining factor ordering." << std::endl;
   SortFactors();

   // initialize three arrays of primal solutions corresponding to factors, one computed in the forward pass, one computed in the backward pass, and one with the best solution obtained so far
   InitializePrimalVector(forwardOrdering_.begin(), forwardOrdering_.end(), forwardPrimal_);
   InitializePrimalVector(forwardOrdering_.begin(), forwardOrdering_.end(), bestForwardPrimal_);
   InitializePrimalVector(backwardOrdering_.begin(), backwardOrdering_.end(), backwardPrimal_);
   InitializePrimalVector(backwardOrdering_.begin(), backwardOrdering_.end(), bestBackwardPrimal_);

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


template<typename FACTOR_ITERATOR, typename OMEGA_ITERATOR>
REAL LP::ComputePass(FACTOR_ITERATOR factorIt, const FACTOR_ITERATOR factorItEnd, OMEGA_ITERATOR omegaIt)
{
   REAL lb = 0.0;
   for(; factorIt!=factorItEnd; ++factorIt, ++omegaIt) {
      lb += UpdateFactor(*factorIt, *omegaIt);
   }
   return lb;
}

template<typename PRIMAL_SOLUTION_STORAGE_ITERATOR>
void ResetPrimalStorage(PRIMAL_SOLUTION_STORAGE_ITERATOR primalIt, const PRIMAL_SOLUTION_STORAGE_ITERATOR primalEndIt)
{
   for(; primalIt!=primalEndIt; ++primalIt) {
      (*primalIt)->reset();
   }
}

template<typename FACTOR_ITERATOR, typename OMEGA_ITERATOR, typename PRIMAL_SOLUTION_STORAGE_ITERATOR>
std::pair<REAL,REAL> LP::ComputePassAndPrimal(FACTOR_ITERATOR factorIt, const FACTOR_ITERATOR factorEndIt, OMEGA_ITERATOR omegaIt, 
      PRIMAL_SOLUTION_STORAGE_ITERATOR primalIt, PRIMAL_SOLUTION_STORAGE_ITERATOR bestPrimalIt, REAL& bestPrimalCost)
{
   // first we need to set primal storage iterator to zero
   ResetPrimalStorage(primalIt, primalIt + INDEX(factorEndIt-factorIt));
   REAL lowerBound = 0.0;
   PRIMAL_SOLUTION_STORAGE_ITERATOR primalItTmp = primalIt;
   FACTOR_ITERATOR factorItTmp = factorIt;
   for(; factorItTmp!=factorEndIt; ++factorItTmp, ++omegaIt, ++primalItTmp) {
      lowerBound += UpdateFactor(*factorItTmp, *omegaIt, *primalItTmp);
   }
   ComputePrimalThroughMessages(factorIt, factorEndIt, primalIt);
   const REAL currentPrimalCost = EvaluatePrimal(factorIt,factorEndIt,primalIt);
   if(currentPrimalCost < bestPrimalCost || bestPrimalCost >= std::numeric_limits<REAL>::max()) { // the second case occurs whenever we have some infeasible primal solution. Possibly, the newer solution, although infeasible, are not as infeasible as the old one.
      bestPrimalCost = currentPrimalCost;
      for(INDEX i=0; i<factorEndIt-factorIt; ++i, ++primalIt, ++bestPrimalIt) {
         std::swap(*(primalIt),*(bestPrimalIt));
      }
   }

   return std::make_pair(lowerBound,bestPrimalCost);
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
         case LPVisitorReturnType::SetAnisotropicReparametrization:
            if(repamMode_ != LPReparametrizationMode::Anisotropic) {
               ComputeAnisotropicWeights(forwardOrdering_, omegaForward_);
               ComputeAnisotropicWeights(backwardOrdering_, omegaBackward_);
               repamMode_ = LPReparametrizationMode::Anisotropic;
            }
            s = v.template visit<LPVisitorReturnType::SetAnisotropicReparametrization>(this);
            break;
         case LPVisitorReturnType::SetRoundingReparametrization:
            if(repamMode_ != LPReparametrizationMode::Rounding) {
               ComputeUniformWeights(forwardOrdering_, omegaForward_);
               ComputeUniformWeights(backwardOrdering_, omegaBackward_);
               repamMode_ = LPReparametrizationMode::Rounding;
            }
            s = v.template visit<LPVisitorReturnType::SetRoundingReparametrization>(this);
            break;
         case LPVisitorReturnType::Reparametrize:
            //std::cout << "reparametrize in " << (repamMode_ == LPReparametrizationMode::Rounding ? " rounding " : " anisotropic ") << "mode\n";
            ComputePass(forwardOrdering_.begin(), forwardOrdering_.end(), omegaForward_.begin());
            ComputePass(backwardOrdering_.begin(), backwardOrdering_.end(), omegaBackward_.begin());
            s = v.template visit<LPVisitorReturnType::Reparametrize>(this);
            break;
         case LPVisitorReturnType::ReparametrizeAndComputePrimal:
            //std::cout << "reparametrize in " << (repamMode_ == LPReparametrizationMode::Rounding ? " rounding " : " anisotropic ") << "mode\n";
            ComputePassAndPrimal(forwardOrdering_.begin(), forwardOrdering_.end(), omegaForward_.begin(), forwardPrimal_.begin(), bestForwardPrimal_.begin(), bestForwardPrimalCost_); 
            ComputePassAndPrimal(backwardOrdering_.begin(), backwardOrdering_.end(), omegaBackward_.begin(), backwardPrimal_.begin(), bestBackwardPrimal_.begin(), bestBackwardPrimalCost_); 
            s = v.template visit<LPVisitorReturnType::ReparametrizeAndComputePrimal>(this);
            break;
         case LPVisitorReturnType::Break:
            s = v.template visit<LPVisitorReturnType::Break>(this);
            return 0;
            break;
         case LPVisitorReturnType::Error:
            return -1;
            break;
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
   std::map<INDEX, FactorTypeAdapter*> indexToFactor; // do zrobienia: this should equal f_. Modify BuildIndexMaps
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
      std::vector<PrimalSolutionStorageAdapter*> cur_primal_adjacent;
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
   BuildIndexMaps(factorIt,factorEndIt,factorToIndex,indexToFactor);

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
void LP::InitializePrimalVector(FACTOR_ITERATOR factorIt, const FACTOR_ITERATOR factorEndIt, std::vector<PrimalSolutionStorageAdapter*>& v)
{
   v.clear(); 
   v.reserve(factorEndIt - factorIt);
   for(; factorIt!=factorEndIt; factorIt++) {
      v.push_back((*factorIt)->AllocatePrimalSolutionStorage());
   }
}

// return best primal solution in the order of f_. Note: this is not a deep copy, but just a vector of pointers to solutions stored in the LP class
inline std::vector<PrimalSolutionStorageAdapter*> LP::GetBestPrimal() const
{
   const bool forward = bestForwardPrimalCost_ < bestBackwardPrimalCost_;
   const std::vector<PrimalSolutionStorageAdapter*>& bestPrimal = forward ? bestForwardPrimal_ : bestBackwardPrimal_;
   const std::vector<FactorTypeAdapter*>& f = forward ? forwardOrdering_ : backwardOrdering_;

   std::vector<PrimalSolutionStorageAdapter*> primalSorted(bestPrimal.size());

   auto perm = ComputePermutation(f_,f);
   for(INDEX i=0; i<perm.size(); ++i) {
      primalSorted[i] = bestPrimal[perm[i]];
   }

   return primalSorted;
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


} // end namespace LP_MP

#endif // LP_MP_MAIN

