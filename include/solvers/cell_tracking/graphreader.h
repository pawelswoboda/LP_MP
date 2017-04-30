#ifndef GRAPH_READER
#define GRAPH_READER 

#include <vector>
#include <map>

namespace dpct
{

  class GraphBuilder;

  // ----------------------------------------------------------------------------------------
  /**
   * @brief Base class for all graph readers implementing common constants and functions
   * 
   */
  class GraphReader
  {
    public:
      typedef double ValueType;
      typedef std::vector<ValueType> FeatureVector;
      typedef std::vector<FeatureVector> StateFeatureVector;

    protected:
      /// Enumerate the strings for attributes used in the Json files and python dicts
      enum class JsonTypes {Segmentations, 
        Links, 
        Exclusions, 
        LinkResults, 
        DivisionResults,
        DetectionResults,
        SrcId, 
        DestId, 
        Value, 
        DivisionValue, 
        Id, 
        Timestep,
        Features, 
        DivisionFeatures,
        AppearanceFeatures,
        DisappearanceFeatures,
        DisappearanceTarget,
        Weights,
        // settings-related
        Settings,
        StatesShareWeights,
        OptimizerEpGap,
        OptimizerVerbose,
        OptimizerNumThreads,
        AllowPartialMergerAppearance,
        RequireSeparateChildrenOfDivision
      };

      /// mapping from JsonTypes to strings which are used in the Json files
      static std::map<JsonTypes, std::string> JsonTypeNames;

    public:
      /**
       * @return the energy of the initial state of flow based solving: no objects tracked at all
       */
      const double getInitialStateEnergy() const { return initialStateEnergy_; }

    protected:
      GraphReader(GraphBuilder* graphBuilder):
        graphBuilder_(graphBuilder)
    {}

      FeatureVector weightedSumOfFeatures(const StateFeatureVector& stateFeatures, const FeatureVector& weights, size_t offset, bool statesShareWeights);
      FeatureVector costsToScoreDeltas(const FeatureVector& costs);
      ValueType costsToScoreDelta(const FeatureVector& costs);

    protected:

      /// the weight vector loaded from file
      FeatureVector weights_;

      /// store the model's energy of the initial solution = all zeros
      double initialStateEnergy_;

      /// the graph builder instance which provides the graph specific node/arc 
      /// creation methods as well as result parsing
      GraphBuilder* graphBuilder_;
  };

} // end namespace dpct

#endif // GRAPH_READER
