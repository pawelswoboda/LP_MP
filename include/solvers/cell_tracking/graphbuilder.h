#ifndef GRAPH_BUILDER
#define GRAPH_BUILDER

#include <vector>
#include <map>

namespace dpct
{
  // ----------------------------------------------------------------------------------------
  /**
   * @brief Base class for graph building factories used by the JsonGraphReader
   */
  class GraphBuilder {
    public:
      typedef double ValueType;
      typedef std::vector<ValueType> CostDeltaVector;
      typedef std::map<size_t, size_t> NodeValueMap;
      typedef std::map<size_t, bool> DivisionValueMap;
      typedef std::map<size_t, size_t> AppearanceValueMap;
      typedef std::map<size_t, size_t> DisappearanceValueMap;
      typedef std::map<std::pair<size_t, size_t>, size_t> ArcValueMap;


      /**
       * @brief add a node which can be indexed by its id. Costs 
       */
      virtual void addNode(
          size_t id,
          const CostDeltaVector& detectionCosts,
          const CostDeltaVector& detectionCostDeltas, 
          const CostDeltaVector& appearanceCostDeltas, 
          const CostDeltaVector& disappearanceCostDeltas,
          size_t targetIdx) = 0;

      /**
       * @brief Specify in which timestep a node lies. Must be called before adding the respective node
       */
      void setNodeTimesteps(size_t id, std::pair<size_t, size_t> timesteps)
      {
        idToTimestepsMap_[id] = timesteps;
      }

      /**
       * @brief add a move arc between the given nodes with the specified cost deltas
       */
      virtual void addArc(size_t srcId, size_t destId, const CostDeltaVector& costDeltas) = 0;

      /**
       * @brief allow a node to divide with the specified cost change
       */
      virtual void allowMitosis(size_t id, ValueType divisionCostDelta) = 0;

      /**
       * @brief return a mapping from node id to value in solution
       * @details doesn't necessarily contain all nodes, implicitly assumes that all missing ones have value zero
       * @return copy/move constructed node value map
       */
      virtual NodeValueMap getNodeValues() = 0;

      /**
       * @brief return a mapping from (src node id, target node id) to value in solution
       * @details doesn't necessarily contain all links, implicitly assumes that all missing ones have value zero
       * @return copy/move constructed link value map
       */
      virtual ArcValueMap getArcValues() = 0;

      /**
       * @brief return a mapping from node id to division value in solution
       * @details doesn't necessarily contain all nodes, implicitly assumes that all missing ones are not dividing
       * @return copy/move constructed division value map
       */
      virtual DivisionValueMap getDivisionValues() = 0;

    protected:
      /// mapping from id to timesteps
      std::map<size_t, std::pair<size_t, size_t> > idToTimestepsMap_;
  };

} // end namespace dpct

#endif
