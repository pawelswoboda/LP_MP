#include "solvers/cell_tracking/graphreader.h"
#include <assert.h>
#include <fstream>
#include <iostream>

namespace dpct
{

  std::map<GraphReader::JsonTypes, std::string> GraphReader::JsonTypeNames = 
  {
    {JsonTypes::Segmentations, "segmentationHypotheses"}, 
    {JsonTypes::Links, "linkingHypotheses"}, 
    {JsonTypes::Exclusions, "exclusions"},
    {JsonTypes::LinkResults, "linkingResults"},
    {JsonTypes::DivisionResults, "divisionResults"},
    {JsonTypes::DetectionResults, "detectionResults"},
    {JsonTypes::SrcId, "src"}, 
    {JsonTypes::DestId, "dest"}, 
    {JsonTypes::Value, "value"},
    {JsonTypes::DivisionValue, "divisionValue"},
    {JsonTypes::Id, "id"}, 
    {JsonTypes::Timestep, "timestep"}, 
    {JsonTypes::Features, "features"},
    {JsonTypes::DivisionFeatures, "divisionFeatures"},
    {JsonTypes::AppearanceFeatures, "appearanceFeatures"},
    {JsonTypes::DisappearanceFeatures, "disappearanceFeatures"},
    {JsonTypes::DisappearanceTarget, "disappTarget"},
    {JsonTypes::Weights, "weights"},
    {JsonTypes::StatesShareWeights, "statesShareWeights"},
    {JsonTypes::Settings, "settings"},
    {JsonTypes::OptimizerEpGap, "optimizerEpGap"},
    {JsonTypes::OptimizerVerbose, "optimizerVerbose"},
    {JsonTypes::OptimizerNumThreads, "optimizerNumThreads"},
    {JsonTypes::AllowPartialMergerAppearance, "allowPartialMergerAppearance"},
    {JsonTypes::RequireSeparateChildrenOfDivision, "requireSeparateChildrenOfDivision"}
  };

  GraphReader::FeatureVector GraphReader::weightedSumOfFeatures(
      const StateFeatureVector& stateFeatures, 
      const FeatureVector& weights,
      size_t offset, 
      bool statesShareWeights)
  {
    FeatureVector costPerState(stateFeatures.size(), 0.0);

    size_t weightIdx = offset;
    for(size_t state = 0; state < stateFeatures.size(); state++)
    {
      for(auto f : stateFeatures[state])
      {
        costPerState[state] = f * weights[(weightIdx++)];
      }

      if(statesShareWeights)
        weightIdx = offset;
    }

    initialStateEnergy_ += costPerState[0];

    return costPerState;
  }

  GraphReader::FeatureVector GraphReader::costsToScoreDeltas(const FeatureVector& costs)
  {
    FeatureVector result;
    for(size_t i = 1; i < costs.size(); i++)
    {
      result.push_back(costs[i] - costs[i-1]);
      if(i > 1 && result[i-1] <= result[i-2])
        std::cout << "Warning: found potentially problematic score setup: " << result[i-1] << " <= " << result[i-2] << std::endl;
    }
    return result;
  }

  GraphReader::ValueType GraphReader::costsToScoreDelta(const FeatureVector& costs)
  {
    assert(costs.size() == 2);
    return costs[1] - costs[0];
  }

} // end namespace dpct
