#ifndef LP_MP_MULTICUT_GLOBAL_FACTOR_HXX
#define LP_MP_MULTICUT_GLOBAL_FACTOR_HXX

#include "LP_MP.h"
#include "union_find.hxx"

namespace LP_MP {

// this factor is a dummy global factor used for checking whether all constraints in multicut hold. Used for primal solution computation.
// Also better rounding schemes can be put here. Look into UltrametricRounding
// dummy reparametrization. Size is used to adjust the primal solution length. Not the nicest way to approach the problem
template<typename FACTOR_CONTAINER> // do zrobienia: assert that factor is multicut
class MulticutGlobalRepamStorage {
public:
   template<typename FACTOR_TYPE>
   MulticutGlobalRepamStorage(const FACTOR_TYPE&, const INDEX n) : N(n) {}
   const REAL operator[](const INDEX i) const { assert(false); return 0.0; }
   /*
   const REAL operator[](const INDEX i) const {
      assert(false);
   }
   REAL& operator[](const INDEX i) {
      assert(false);
   }
   */
   INDEX size() const { return N; }
   void ResizeRepam(const INDEX n) { N = n; }
private:
   INDEX N;
};

class MulticutGlobalFactor 
{
public:
   MulticutGlobalFactor() {};

   INDEX AddEdge(const INDEX i1, const INDEX i2) 
   {
      edges_.push_back(std::make_tuple(i1,i2));
      noNodes_ = std::max(noNodes_, std::max(i1,i2)+1);
      return edges_.size()-1;
   }

   template<typename REPAM_ARRAY>
   void MaximizePotential(const REPAM_ARRAY& repam) {};
   template<typename REPAM_ARRAY>
   REAL LowerBound(const REPAM_ARRAY& repamPot) const {
      assert(repamPot.size() == edges_.size());
      return 0;
   }

   INDEX size() const { return 0; }

   template<typename REPAM_ARRAY>
   REAL EvaluatePrimal(const REPAM_ARRAY& repam, const PrimalSolutionStorage::Element primal) const
   {
      assert(repam.size() == edges_.size());
      // search for violated cycles. If there is one, return infty, otherwise return 0;
      // have k sets initialy, where k is number of nodes. Join two sets containing these nodes whenever there is a connecting edge with primal value false.
      // go through all edges with primal true and see whether nodes are in same component. If yes, return infinity.
      UnionFind uf(noNodes_);
      for(INDEX e=0; e<edges_.size(); ++e) {
         if(primal[e] == false) {
            const INDEX i1 = std::get<0>(edges_[e]);
            const INDEX i2 = std::get<1>(edges_[e]);
            // connect components with node i1 and with node i2
            uf.merge(i1,i2);
         }
      }
      for(INDEX e=0; e<edges_.size(); ++e) {
         if(primal[e] == true) {
            const INDEX i1 = std::get<0>(edges_[e]);
            const INDEX i2 = std::get<1>(edges_[e]);
            // there may not be a path from i1 to i2 consisting of edges with primal value false only
            if(uf.connected(i1,i2)) {
               return std::numeric_limits<REAL>::max();
            }
         }
      }
      return 0.0;
   }
   void WritePrimal(const PrimalSolutionStorage::Element primal, std::ofstream& fs) const
   {
      //fs << primal;
   }

private:
   std::vector<std::tuple<INDEX,INDEX>> edges_;
   INDEX noNodes_ = 0;
};

class MulticutUnaryGlobalMessage
{
public:
   MulticutUnaryGlobalMessage(const INDEX i) : i_(i) {}
   void ComputeRightFromLeftPrimal(const PrimalSolutionStorage::Element leftPrimal, PrimalSolutionStorage::Element rightPrimal)
   {
      // not nice: rightPrimal may not be of the correct size. This should be better treated in primal initialization, i.e. reset function
      // one possibility is to have dummy reparametrization of correct size which is however not writable and then we can collect the primal potentials as done. Other approach might be multicut constructor handling the constraints
      //assert(false);
      /*
      if(i_ >= rightPrimal.size()) {
         rightPrimal.resize(i_+1);
      }
      */
      rightPrimal[i_] = leftPrimal[0];
   }
   // dummy functions: they do nothing. Allow for no Repam{LEFT|RIGHT} in factor_messages to indicate not reparametrized factor
   template<typename G>
   void RepamLeft(G& repamPot, const REAL msg, const INDEX msg_dim)
   {}
   template<typename G>
   void RepamRight(G& repamPot, const REAL msg, const INDEX msg_dim)
   {}
   

private:
   const INDEX i_; // index of the affected variable in the cycle factor.
};

/*
class LiftedMulticutGlobalFactor : MulticutGlobalFactor{
public:
private:
   enum class EdgeType {Base,Lifted,Tightening};
   struct LiftedMulticutEdge {INDEX i; INDEX j; EdgeType t;}
   std::vector
    std::vector<Edge> liftedEdges_;
};
*/

} // end namespace LP_MP

#endif // LP_MP_MULTICUT_GLOBAL_FACTOR




