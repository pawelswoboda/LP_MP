#include "LP_MP.h"

namespace LP_MP {

LP::LP() 
{}

LP::~LP()
{
   for(INDEX i=0; i<m_.size(); i++) { delete m_[i]; }
   for(INDEX i=0; i<f_.size(); i++) { delete f_[i]; }
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
   //backwardOrdering_ = fSorted;
   //std::reverse(backwardOrdering_.begin(), backwardOrdering_.end());
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

void LP::UpdateFactor(FactorTypeAdapter* f, const std::vector<REAL>& omega)
{
   f->UpdateFactor(omega);
}

void LP::UpdateFactor(FactorTypeAdapter* f, const std::vector<REAL>& omega, PrimalSolutionStorage::Element primalIt)
{
   f->UpdateFactor(omega, primalIt);
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

