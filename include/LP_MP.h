#ifndef LP_MP_MAIN
#define LP_MP_MAIN

#include "instances.inc"
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
#include "help_functions.hxx"
#include "topological_sort.hxx"

// do zrobienia: make those functions const which are logically so
// do zrobienia: introduce MessageConstraINDEX and FactorConstraINDEX
// do zrobienia: introduce LP_MP with MessageTypeCRTP

namespace LP_MP {

// forward declaration
class MessageIterator;

// virtual base class for factors used by LP class
class FactorTypeAdapter
{
public:
   virtual ~FactorTypeAdapter() {};
   virtual void UpdateFactor(const std::vector<REAL>& omega) = 0; // do zrobienia: make this a std::vector<REAL>& omega or some other container
   MessageIterator begin(); 
   MessageIterator end();
   virtual const INDEX GetNoMessages() const = 0;
   virtual FactorTypeAdapter* GetConnectedFactor(const INDEX i) const = 0;
   virtual REAL LowerBound() const = 0;
   virtual std::vector<REAL> GetReparametrizedPotential() const = 0;
};

class MessageTypeAdapter
{
public:
   virtual ~MessageTypeAdapter() {};
   virtual FactorTypeAdapter* GetLeftFactor() const = 0;
   virtual FactorTypeAdapter* GetRightFactor() const = 0;
   virtual void SetMessage(const std::valarray<REAL>& m) = 0; // do zrobienia: change to vector
   virtual const std::valarray<REAL> GetMessage() const = 0; // do zrobienia: change to vector
};

// very primitive iterator class. Access may be slow. A more direct implementation woulb be more complicated, though
// used in main LP class
class MessageIterator //: public std::iterator<random_access_iterator>
{
public:
   MessageIterator(FactorTypeAdapter* const factor, const INDEX msg_idx) : factor_(factor), msg_idx_(msg_idx) {}
   MessageIterator* operator++() { ++msg_idx_; return this; }
   bool operator==(const MessageIterator& rhs) { return (factor_ == rhs.factor_ && msg_idx_ == rhs.msg_idx_); }
   bool operator!=(const MessageIterator& rhs) { return !operator==(rhs); }
   FactorTypeAdapter& operator*() { return *(factor_->GetConnectedFactor(msg_idx_)); }
   FactorTypeAdapter* operator->() { return factor_->GetConnectedFactor(msg_idx_); }
private:
   FactorTypeAdapter* const factor_;
   INDEX msg_idx_;
};

inline MessageIterator FactorTypeAdapter::begin() { return MessageIterator(this,0); }
inline MessageIterator FactorTypeAdapter::end()  { return MessageIterator(this, GetNoMessages()); };


// do zrobienia: return this for AddFactor, AddMessage and where appropriate
class LP
{
public:
   LP();
   LP(const INDEX noFactors, const INDEX noMessages);
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

   void ComputeWeights(const std::vector<FactorTypeAdapter*>& f, std::vector<std::vector<REAL> >& omega);
   void ComputeUniformWeights(const std::vector<FactorTypeAdapter*>& f, std::vector<std::vector<REAL> >& omega);

   REAL LowerBound();
   REAL UpdateFactor(FactorTypeAdapter* f, const std::vector<REAL>& omega); // perform one block coordinate step for factor f
   REAL ComputePass(const std::vector<FactorTypeAdapter*>& f, const std::vector<std::vector<REAL> >& omega);
   void Solve(const INDEX noIter);

   std::vector<std::vector<REAL> > GetReparametrizedModel() const;
private:
   // do zrobienia: possibly hold them in shared_ptr?
   std::vector<FactorTypeAdapter*> f_;
   std::vector<MessageTypeAdapter*> m_;
   std::vector<FactorTypeAdapter*> forwardOrdering_, backwardOrdering_;
   std::vector<std::vector<REAL> > omegaForward_, omegaBackward_;
   std::vector<std::pair<FactorTypeAdapter*, FactorTypeAdapter*> > factorRel_; // factor ordering relations
};

} // end namespace LP_MP


#endif // LP_MP_MAIN

