#include "LP_MP.h"
#include <chrono>

namespace LP_MP {

LP::LP() 
   : f_(0), m_(0), forwardOrdering_(0), backwardOrdering_(0), omegaForward_(0), omegaBackward_(0), factorRel_(0)
{}

LP::LP(const INDEX noFactors, const INDEX noMessages) 
   : f_(0), m_(0), forwardOrdering_(0), backwardOrdering_(0), omegaForward_(0), omegaBackward_(0), factorRel_(0)
{
   f_.reserve(noFactors);
   m_.reserve(noMessages);
   forwardOrdering_.reserve(noFactors);
   backwardOrdering_.reserve(noFactors);
   omegaForward_.reserve(noFactors);
   omegaBackward_.reserve(noFactors);
}

LP::~LP()
{
   for(INDEX i=0; i<m_.size(); i++) {
      delete m_[i];
   }
   for(INDEX i=0; i<f_.size(); i++) {
      delete f_[i];
   }
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
   BuildIndexMaps(f_,factorToIndex,indexToFactor);

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

std::vector<std::vector<FactorTypeAdapter*> > LP::ComputeFactorConnection(const std::vector<FactorTypeAdapter* >& f)
{
   std::vector<std::vector<FactorTypeAdapter*> > fc(f.size());
   for(INDEX i=0; i<f.size(); ++i) {
      for(auto mIt=f[i]->begin(); mIt!=f[i]->end(); ++mIt) {
         fc[i].push_back( &*mIt ); // dereference iterator to factor and then take address of factor
      }
   }

   for(INDEX i=0; i<fc.size(); ++i) {
      assert(HasUniqueValues(fc[i]));
      for(INDEX j=0; j<fc[i].size(); ++j) {
         assert(fc[i][j] != f[i]);
      }
   }

   return fc;
}

void LP::ComputeWeights(const std::vector<FactorTypeAdapter*>& f, std::vector<std::vector<REAL> >& omega)
{
   assert(f.size() == f_.size());
   std::vector<std::vector<FactorTypeAdapter*> > fc = ComputeFactorConnection(f);
   std::vector<std::vector<bool> > fcAccessedLater(f.size());

   std::map<FactorTypeAdapter*, INDEX> factorToIndex;
   std::map<INDEX, FactorTypeAdapter*> indexToFactor;
   BuildIndexMaps(f,factorToIndex,indexToFactor);

   for(INDEX i=0; i<f.size(); i++) {
      assert(factorToIndex.find(f[i]) != factorToIndex.end());
      const INDEX index1 = factorToIndex[f[i]];
      assert(INDEX(index1) == i);
      for(INDEX j=0; j<fc[i].size(); j++) {
         fcAccessedLater[i].resize(fc[i].size(),false);
         assert(factorToIndex.find(fc[i][j]) != factorToIndex.end());
         const INDEX index2 = factorToIndex[ fc[i][j] ];
         bool INDEXermedFactor = false;
         for(auto fIt2=fc[index2].begin(); fIt2!=fc[index2].end(); ++fIt2) {
            assert(factorToIndex.find(*fIt2) != factorToIndex.end());
            const INDEX index3 = factorToIndex[*fIt2];
            if(index1 < index3 && index1 != index3)
               INDEXermedFactor = true;
         }
         if(index1 < index2 || INDEXermedFactor == true) {
         //if(index1 < index2 ) {
            fcAccessedLater[i][j] = true;
         }
      }
   }

   omega.clear();
   omega.resize(f.size(),std::vector<REAL>(0));
   for(INDEX i=0; i<f.size(); i++) {
      omega[i].resize(fc[i].size(), 0.0);
      INDEX noFactorsAccessedLater = 0;
      for(INDEX j=0;j<fc[i].size(); j++) {
         if(fcAccessedLater[i][j]) {
            noFactorsAccessedLater++;
         }
      }
      //const REAL weight = 1.0 / REAL(noFactorsAccessedLater);
      //const REAL weight = 0.8*1.0 / REAL(fcAccessedLater[i].size());
      const REAL weight = 1.0 / (std::max(REAL(noFactorsAccessedLater), REAL(fcAccessedLater[i].size() - noFactorsAccessedLater)) );
      //const REAL weight = 0.1; // 0.5 works well for pure assignment with equality messages
      for(INDEX j=0; j<fc[i].size(); j++) {
         if(fcAccessedLater[i][j]) {
            omega[i][j] = weight;
         } else {
            omega[i][j] = 0.0;
         }
      }
      assert( std::accumulate(omega[i].begin(), omega[i].end(),0.0) < 1.0 + 1e-7);
   }
}

// compute uniform weights so as to help decoding for obtaining primal solutions
void LP::ComputeUniformWeights(const std::vector<FactorTypeAdapter*>& f, std::vector<std::vector<REAL> >& omega)
{
   assert(f.size() == f_.size());
   std::vector<std::vector<FactorTypeAdapter*> > fc = ComputeFactorConnection(f);

   omega.clear();
   omega.resize(f.size());
   for(INDEX i=0; i<f.size(); i++) {
      omega[i] = std::vector<REAL>(fc[i].size(), 1.0/REAL(2*fc[i].size() + 1.0) );
   }
}

void LP::Init()
{
   std::cout << "Determining factor ordering." << std::endl;
   SortFactors();
   std::cout << "Computing weights." << std::endl;
   ComputeWeights(forwardOrdering_, omegaForward_);
   ComputeWeights(backwardOrdering_, omegaBackward_);

   std::cout << "Initial lower bound before optimizing = " << LowerBound() << std::endl;
}

REAL LP::LowerBound()
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

REAL LP::ComputePass(const std::vector<FactorTypeAdapter*>& f, const std::vector<std::vector<REAL> >& omega)
{
   REAL lb = 0.0;
   for(INDEX i=0; i<f.size(); i++) {
      lb += UpdateFactor(f[i], omega[i]);
   }
   return lb;
}

void LP::Solve(const INDEX noIter)
{
   std::cout << "Current vectorization implementation ";
   switch(VC_IMPL) {
      case Vc::Implementation::ScalarImpl:
         std::cout << "scalar"; break;
      case Vc::Implementation::SSE3Impl:
            std::cout << "SSE3"; break;
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
            std::cout << "vectorization support unknown";
   }

   std::cout << "\n";
   std::cout << "vector size = " << REAL_SIMD::Size << "\n";
   Init();

   auto t_begin = std::chrono::steady_clock::now();

   REAL prevLowerBound = -std::numeric_limits<REAL>::max();
   INDEX iter;
   for(iter=0; iter<noIter; iter++) {
      REAL lbForward = ComputePass(forwardOrdering_, omegaForward_);
      REAL lbBackward = ComputePass(backwardOrdering_, omegaBackward_);
      if(iter%10 == 0) {
         const REAL lowerBound = LowerBound();
         std::cout << "Iteration = " << iter << "\n"; // ", lower bound after forward move = " << lbForward << ", lower bound after backward move = " << lbBackward << std::endl;
         std::cout << "true lower bound = " << lowerBound << std::endl;
         // early stopping because of lacking progress
         if( std::abs(lowerBound - prevLowerBound)/(std::abs(0.001*lowerBound) + 1.0) < 0.0001 ) {
            std::cout << "No progress, early stopping" << std::endl; 
            break;
         }
         prevLowerBound = lowerBound;
      }
   }

   auto t_end = std::chrono::steady_clock::now();
   std::cout << "Optimization took " << std::chrono::duration_cast<std::chrono::milliseconds>(t_end - t_begin).count() << " milliseconds and " << iter << " iterations\n";

   ComputeUniformWeights(forwardOrdering_, omegaForward_);
   ComputeUniformWeights(backwardOrdering_, omegaBackward_);
   std::cout << "Compute rounding reparametrization" << std::endl;
   for(INDEX iter=0; iter<1; iter++) {
      ComputePass(forwardOrdering_, omegaForward_);
      ComputePass(backwardOrdering_, omegaBackward_);
   }
   std::cout << "true lower bound = " << LowerBound() << std::endl;


}

std::vector<std::vector<REAL> > LP::GetReparametrizedModel() const
{
   std::vector<std::vector<REAL> > repam(f_.size());
   for(INDEX i=0; i<f_.size(); i++) {
      repam[i] = f_[i]->GetReparametrizedPotential();
   }
   return repam;
}

} // end namespace LP_MP

