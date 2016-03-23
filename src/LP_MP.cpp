#include "LP_MP.h"
#include "tolerance.hxx"

namespace LP_MP {

LP::LP() 
{}

LP::~LP()
{
   for(INDEX i=0; i<m_.size(); i++) { delete m_[i]; }
   for(INDEX i=0; i<f_.size(); i++) { delete f_[i]; }
   for(INDEX i=0; i<forwardPrimal_.size(); i++) { delete forwardPrimal_[i]; }
   for(INDEX i=0; i<backwardPrimal_.size(); i++) { delete backwardPrimal_[i]; }
   for(INDEX i=0; i<bestForwardPrimal_.size(); i++) { delete bestForwardPrimal_[i]; }
   for(INDEX i=0; i<bestBackwardPrimal_.size(); i++) { delete bestBackwardPrimal_[i]; }
}

INDEX LP::AddFactor(FactorTypeAdapter* f) 
{ 
   f_.push_back(f);
   return f_.size() - 1;
}

INDEX LP::GetNumberOfFactors() const { return f_.size(); }
FactorTypeAdapter* LP::GetFactor(const INDEX i) const { return f_[i]; }

INDEX LP::AddMessage(MessageTypeAdapter* m) 
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

MessageTypeAdapter* LP::GetMessage(const INDEX i) const { return m_[i]; }
INDEX LP::GetNumberOfMessages() const { return m_.size(); }

void LP::AddFactorRelation(FactorTypeAdapter* f1, FactorTypeAdapter* f2)
{
   factorRel_.push_back(std::make_pair(f1,f2));
}

void LP::SortFactors()
{
   // assume that factorRel_ describe a DAG. Compute topological sorting
   Topological_Sort::Graph g(f_.size());

   std::map<FactorTypeAdapter*,INDEX> factorToIndex; // possibly do it with a hash_map for speed
   std::map<INDEX,FactorTypeAdapter*> indexToFactor;
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
   assert(HasUniqueValues(fSorted));
   forwardOrdering_  = fSorted;
   backwardOrdering_ = fSorted;
   std::reverse(backwardOrdering_.begin(), backwardOrdering_.end());
}

// only compute factors adjacent to which also messages can be send
std::vector<std::vector<FactorTypeAdapter*> > LP::ComputeSendFactorConnection(const std::vector<FactorTypeAdapter* >& f)
{
   std::vector<std::vector<FactorTypeAdapter*> > fc(f.size());
   for(INDEX i=0; i<f.size(); ++i) {
      for(auto mIt=f[i]->begin(); mIt!=f[i]->end(); ++mIt) {
         if(mIt.CanSendMessage()) {
            fc[i].push_back( mIt.GetConnectedFactor() );
         }
      }
   }

   return fc;
}

std::vector<std::vector<FactorTypeAdapter*> > LP::ComputeFactorConnection(const std::vector<FactorTypeAdapter* >& f)
{
   std::vector<std::vector<FactorTypeAdapter*> > fc(f.size());
   for(INDEX i=0; i<f.size(); ++i) {
      for(auto mIt=f[i]->begin(); mIt!=f[i]->end(); ++mIt) {
         fc[i].push_back( mIt.GetConnectedFactor() );
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
void LP::ComputeAnisotropicWeights(const std::vector<FactorTypeAdapter*>& f, std::vector<std::vector<REAL> >& omega)
{
   assert(f.size() == f_.size());
   std::vector<std::vector<FactorTypeAdapter*> > fc = ComputeSendFactorConnection(f);
   std::vector<std::vector<bool> > fcAccessedLater(f.size());

   std::map<FactorTypeAdapter*, INDEX> factorToIndex;
   std::map<INDEX, FactorTypeAdapter*> indexToFactor;
   BuildIndexMaps(f.begin(), f.end(), factorToIndex,indexToFactor);

   for(INDEX i=0; i<f.size(); i++) {
      assert(factorToIndex.find(f[i]) != factorToIndex.end());
      const INDEX index1 = factorToIndex[f[i]];
      assert(INDEX(index1) == i);
      fcAccessedLater[i].resize(fc[i].size(),false);
      for(INDEX j=0; j<fc[i].size(); j++) {
         assert(factorToIndex.find(fc[i][j]) != factorToIndex.end());
         const INDEX index2 = factorToIndex[ fc[i][j] ];
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
            fcAccessedLater[i][j] = true;
         }
      }
   }
   
   omega.clear();
   omega.resize(f.size());
   if(m_.size() == 0) { 
      std::cout << "no messages in problem\n"; 
      return;
   }
   for(INDEX i=0; i<f.size(); i++) {
      omega[i].resize(fc[i].size(), 0.0);
      INDEX noFactorsAccessedLater = 0;
      for(INDEX j=0;j<fc[i].size(); j++) {
         if(fcAccessedLater[i][j]) {
            noFactorsAccessedLater++;
         }
      }
      // this was for custom defined message weights. Currently not supported.
      //std::vector<INDEX> messageClassNumber(omega[i].size());
      //for(INDEX j=0; j<messageClassNumber.size(); ++j) {
      //   messageClassNumber[j] = f[i]->GetMessage(j)->GetMessageNumber();
      //}
      //std::vector<INDEX> numberMessagesOfClass(*std::max_element(messageClassNumber.begin(), messageClassNumber.end())+1, 0);
      //for(INDEX j=0; j<messageClassNumber.size(); ++j) {
      //   if(fcAccessedLater[i][j]) {
      //      numberMessagesOfClass[messageClassNumber[j]] += 1;
      //   }
      //}
      //std::vector<REAL> minimumWeight(*std::max_element(messageClassNumber.begin(), messageClassNumber.end())+1, 0.0);
      //for(INDEX j=0; j<messageClassNumber.size(); ++j) {
      //   if(fcAccessedLater[i][j]) {
      //      const INDEX c = messageClassNumber[j];
      //      auto m = f[i]->GetMessage(j);
      //      if(f[i] == m->GetLeftFactor()) {
      //         minimumWeight[c] = m->GetMessageWeightToRight();
      //      } else {
      //         assert(f[i] == m->GetRightFactor());
      //         minimumWeight[c] = m->GetMessageWeightToLeft();
      //      }
      //   }
      //}
      //assert(std::accumulate(minimumWeight.begin(), minimumWeight.end(), 0.0) <= 1.0 + eps);
      // now distribute slack weights evenly across messages
      //const REAL slackWeight = 1.0 - std::accumulate(minimumWeight.begin(), minimumWeight.end(), 0.0); // do zrobienia: take into account weight construction below of traditional weight as in SRMP
      const INDEX numberActiveMessages = std::count(fcAccessedLater[i].begin(), fcAccessedLater[i].end(), true);

      std::vector<REAL> weights(omega[i].size(),0.0);
      //const REAL weight = 1.0 / REAL(noFactorsAccessedLater);
      //const REAL weight = 0.8*1.0 / REAL(fcAccessedLater[i].size());
      const REAL weight = 1.0 / (std::max(REAL(noFactorsAccessedLater), REAL(fcAccessedLater[i].size() - noFactorsAccessedLater)) );
      //const REAL weight = 0.1; // 0.5 works well for pure assignment with equality messages
      for(INDEX j=0; j<fc[i].size(); j++) {
         if(fcAccessedLater[i][j]) {
            //const INDEX c = messageClassNumber[j];
            // do zrobienia: denominator should count how many messages of this class are accessed later
            //omega[i][j] = minimumWeight[c]/REAL(numberMessagesOfClass[c]) + slackWeight/REAL(numberActiveMessages);
            omega[i][j] = weight; 
         } else {
            omega[i][j] = 0.0;
         }
      }
      assert( std::accumulate(omega[i].begin(), omega[i].end(),0.0) <= 1.0 + eps);
   }
}

// compute uniform weights so as to help decoding for obtaining primal solutions
void LP::ComputeUniformWeights(const std::vector<FactorTypeAdapter*>& f, std::vector<std::vector<REAL> >& omega)
{
   assert(f.size() == f_.size());
   std::vector<std::vector<FactorTypeAdapter*> > fc = ComputeSendFactorConnection(f);
   assert(f.size() == fc.size());

   omega.clear();
   omega.resize(f.size());
   for(INDEX i=0; i<f.size(); i++) {
      omega[i] = std::vector<REAL>(fc[i].size(), 1.0/REAL(fc[i].size() + 1.0) );
   }
}

REAL LP::LowerBound() const
{
   REAL lb = 0.0;
   for(auto fIt=f_.begin(); fIt!=f_.end(); fIt++) {
      lb += (*fIt)->LowerBound();
      assert( (*fIt)->LowerBound() > -10000000.0);
   }
   return lb;
}

REAL LP::UpdateFactor(FactorTypeAdapter* f, const std::vector<REAL>& omega)
{
   f->UpdateFactor(omega);
   return 0.0;
}

REAL LP::UpdateFactor(FactorTypeAdapter* f, const std::vector<REAL>& omega, PrimalSolutionStorageAdapter* primal)
{
   f->UpdateFactor(omega, primal);
   return 0.0;
}

std::vector<std::vector<REAL> > LP::GetReparametrizedModel(const INDEX begin, const INDEX end) const
{
   assert(end >= begin);
   assert(end <= f_.size());
   std::vector<std::vector<REAL> > repam(end-begin);
   for(INDEX i=begin; i<end; i++) {
      repam[i-begin] = f_[i]->GetReparametrizedPotential();
   }
   return repam;
}
std::vector<std::vector<REAL> > LP::GetReparametrizedModel() const
{
   return GetReparametrizedModel(0,f_.size());
}

} // end namespace LP_MP

