#ifndef JSON_GRAPH_READER
#define JSON_GRAPH_READER 

#include "json/json.h"

#include "graphreader.h"
//#include "trackingalgorithm.h"
//#include "log.h"
#include "graphbuilder.h"

namespace dpct
{

enum class JsonTypes {Segmentations, 
  OptimizerNumThreads,
  AllowPartialMergerAppearance,
  RequireSeparateChildrenOfDivision
};

// ----------------------------------------------------------------------------------------
/**
 * @brief A json graph reader provides functions to read a flow graph from json, 
 * and stores a mapping from JSON ids to graph nodes
 * 
 */
class JsonGraphReader : public GraphReader
{
public:
	/**
	 * @brief Construct a json graph reader that reads from the specified files
	 * @param modelFilename filename of JSON model description
	 * @param weightsFilename filename where the weights are stored as JSON
	 */
	JsonGraphReader(const std::string& modelFilename, const std::string& weightsFilename, GraphBuilder* graphBuilder);

	/**
	 * @brief Add nodes and arcs to the graph builder according to the model file. 
	 * Costs are computed from features times weights
	 */
	void createGraphFromJson();

	/**
	 * @brief Save a resulting flow map back to JSON, the flow is extracted by the graph builder
	 * 
	 * @param filename filename where the results will be stored as Json
	 */
	void saveResultJson(const std::string& filename);

	/**
	 * @return the energy of the initial state of flow based solving: no objects tracked at all
	 */
	//const double getInitialStateEnergy() const { return initialStateEnergy_; }

private:
	StateFeatureVector extractFeatures(const Json::Value& entry, JsonTypes type);
	FeatureVector readWeightsFromJson(const std::string& filename);
	size_t getNumWeights(const Json::Value& jsonHyp, JsonTypes type, bool statesShareWeights);

private:
	/// filename where the model is stored in json
	std::string modelFilename_;

	/// filename where the wheights are stored in json
	std::string weightsFilename_;
};

} // end namespace dpct

#endif // JSON_GRAPH_READER
