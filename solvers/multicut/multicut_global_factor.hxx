#ifndef LP_MP_MULTICUT_GLOBAL_FACTOR_HXX
#define LP_MP_MULTICUT_GLOBAL_FACTOR_HXX

#include "LP_MP.h"
#include "union_find.hxx"

namespace LP_MP {

// this factor is a dummy global factor used for checking whether all constraints in multicut hold. Used for primal solution computation.
// Also better rounding schemes can be put here. Look into UltrametricRounding
class MulticutGlobalFactor 
{
public:
   using PrimalType = std::vector<bool>;
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
      assert(repamPot.size() == 0);
      return 0;
   }

   INDEX size() const { return 0; }

   template<typename REPAM_ARRAY>
   REAL EvaluatePrimal(const REPAM_ARRAY& repam, const std::vector<bool>& primal) const
   {
      assert(primal.size() == edges_.size());
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
   void WritePrimal(const PrimalType& primal, std::ofstream& fs) const
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
   void ComputeRightFromLeftPrimal(const bool leftPrimal, MulticutGlobalFactor::PrimalType& rightPrimal)
   {
      // not nice: rightPrimal may not be of the correct size. This should be better treated in primal initialization, i.e. reset function
      if(i_ >= rightPrimal.size()) {
         rightPrimal.resize(i_+1);
      }
      assert(i_<rightPrimal.size());
      rightPrimal[i_] = leftPrimal;
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

} // end namespace LP_MP

#endif // LP_MP_MULTICUT_GLOBAL_FACTOR




