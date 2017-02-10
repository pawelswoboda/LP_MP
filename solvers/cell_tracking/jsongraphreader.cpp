#include "json_graph_reader.h"
#include <assert.h>
#include <iostream>
#include <fstream>

namespace dpct
{

JsonGraphReader::JsonGraphReader(const std::string& modelFilename, const std::string& weightsFilename, GraphBuilder* graphBuilder):
	GraphReader(graphBuilder),
	modelFilename_(modelFilename),
	weightsFilename_(weightsFilename)
{
}

size_t JsonGraphReader::getNumWeights(const Json::Value& jsonHyp, JsonTypes type, bool statesShareWeights)
{
	size_t numWeights;
	StateFeatureVector stateFeatVec = extractFeatures(jsonHyp, type);
	
	if(statesShareWeights)
		numWeights = stateFeatVec[0].size();
	else
	{
		numWeights = 0;
		for(auto stateFeats : stateFeatVec)
			numWeights += stateFeats.size();
	}
	return numWeights;
}

void JsonGraphReader::createGraphFromJson()
{
	std::ifstream input(modelFilename_.c_str());
	if(!input.good())
		throw std::runtime_error("Could not open JSON model file for reading: " + modelFilename_);

	Json::Value root;
	input >> root;

	// get flag whether states should share weights or not
	bool statesShareWeights = false;
	if(root.isMember(JsonTypeNames[JsonTypes::Settings]))
	{
		Json::Value& settings = root[JsonTypeNames[JsonTypes::Settings]];
		if(settings.isMember(JsonTypeNames[JsonTypes::StatesShareWeights]))
			statesShareWeights = settings[JsonTypeNames[JsonTypes::StatesShareWeights]].asBool();
	}

	// ------------------------------------------------------------------------------
	// get weight vector and number of weights needed for each different variable type
	FeatureVector weights = readWeightsFromJson(weightsFilename_);
	size_t numDetWeights = 0;
	size_t numDivWeights = 0;
	size_t numAppWeights = 0;
	size_t numDisWeights = 0;
	size_t numLinkWeights = 0;

	const Json::Value segmentationHypotheses = root[JsonTypeNames[JsonTypes::Segmentations]];
	for(int i = 0; i < (int)segmentationHypotheses.size(); i++)
	{
		const Json::Value jsonHyp = segmentationHypotheses[i];
		numDetWeights = getNumWeights(jsonHyp, JsonTypes::Features, statesShareWeights);

		if(jsonHyp.isMember(JsonTypeNames[JsonTypes::DivisionFeatures]))
			numDivWeights = getNumWeights(jsonHyp, JsonTypes::DivisionFeatures, statesShareWeights);

		if(jsonHyp.isMember(JsonTypeNames[JsonTypes::AppearanceFeatures]))
			numAppWeights = getNumWeights(jsonHyp, JsonTypes::AppearanceFeatures, statesShareWeights);

		if(jsonHyp.isMember(JsonTypeNames[JsonTypes::DisappearanceFeatures]))
			numDisWeights = getNumWeights(jsonHyp, JsonTypes::DisappearanceFeatures, statesShareWeights);
	}

	const Json::Value linkingHypotheses = root[JsonTypeNames[JsonTypes::Links]];
	for(int i = 0; i < (int)linkingHypotheses.size(); i++)
	{
		const Json::Value jsonHyp = linkingHypotheses[i];
		if(jsonHyp.isMember(JsonTypeNames[JsonTypes::Features]))
				numLinkWeights = getNumWeights(jsonHyp, JsonTypes::Features, statesShareWeights);
	}

	if(weights.size() != numDetWeights + numDivWeights + numAppWeights + numDisWeights + numLinkWeights)
	{
		std::stringstream s;
		s << "Loaded weights do not meet model requirements! Got " << weights.size() << ", need " 
			<< numDetWeights + numDivWeights + numAppWeights + numDisWeights + numLinkWeights;
		throw std::runtime_error(s.str());
	}

	size_t linkWeightOffset = 0;
	size_t detWeightOffset = linkWeightOffset + numLinkWeights;
	size_t divWeightOffset = detWeightOffset + numDetWeights;
	size_t appWeightOffset = divWeightOffset + numDivWeights;
	size_t disWeightOffset = appWeightOffset + numAppWeights;

	// ------------------------------------------------------------------------------
	// read segmentation hypotheses and add to flowgraph
	std::cout << "\tcontains " << segmentationHypotheses.size() << " segmentation hypotheses" << std::endl;
	
	for(int i = 0; i < (int)segmentationHypotheses.size(); i++)
	{
		const Json::Value jsonHyp = segmentationHypotheses[i];
		
		if(!jsonHyp.isMember(JsonTypeNames[JsonTypes::Id]))
			throw std::runtime_error("Cannot read detection hypothesis without Id!");
		size_t id = jsonHyp[JsonTypeNames[JsonTypes::Id]].asInt();

		if(!jsonHyp.isMember(JsonTypeNames[JsonTypes::Features]))
			throw std::runtime_error("Cannot read detection hypothesis without features!");

		if(jsonHyp.isMember(JsonTypeNames[JsonTypes::Timestep]))
		{
			const Json::Value& timeJson = jsonHyp[JsonTypeNames[JsonTypes::Timestep]];
			if(!timeJson.isArray() || timeJson.size() != 2)
				throw std::runtime_error("Node's Timestep is supposed to be a 2-element array");
			std::pair<int, int> timeRange = std::make_pair(timeJson[0].asInt(), timeJson[1].asInt());
			graphBuilder_->setNodeTimesteps(id, timeRange);
		}

		FeatureVector detCosts = weightedSumOfFeatures(extractFeatures(jsonHyp, JsonTypes::Features), weights, detWeightOffset, statesShareWeights);
		FeatureVector detCostDeltas = costsToScoreDeltas(detCosts);
		FeatureVector appearanceCostDeltas;
		FeatureVector disappearanceCostDeltas;
		if(jsonHyp.isMember(JsonTypeNames[JsonTypes::AppearanceFeatures]))
			appearanceCostDeltas = costsToScoreDeltas(weightedSumOfFeatures(extractFeatures(jsonHyp, JsonTypes::AppearanceFeatures), weights, appWeightOffset, statesShareWeights));

		if(jsonHyp.isMember(JsonTypeNames[JsonTypes::DisappearanceFeatures]))
			disappearanceCostDeltas = costsToScoreDeltas(weightedSumOfFeatures(extractFeatures(jsonHyp, JsonTypes::DisappearanceFeatures), weights, disWeightOffset, statesShareWeights));

		size_t targetIdx = 0;
		if(jsonHyp.isMember(JsonTypeNames[JsonTypes::DisappearanceTarget]))
		{
			targetIdx = jsonHyp[JsonTypeNames[JsonTypes::DisappearanceTarget]].asUInt();
		}

		graphBuilder_->addNode(id, detCosts, detCostDeltas, appearanceCostDeltas, disappearanceCostDeltas, targetIdx);
	}

	// read linking hypotheses
	std::cout << "\tcontains " << linkingHypotheses.size() << " linking hypotheses" << std::endl;
	for(int i = 0; i < (int)linkingHypotheses.size(); i++)
	{
		const Json::Value jsonHyp = linkingHypotheses[i];

		size_t srcId = jsonHyp[JsonTypeNames[JsonTypes::SrcId]].asInt();
		size_t destId = jsonHyp[JsonTypeNames[JsonTypes::DestId]].asInt();
		graphBuilder_->addArc(srcId, destId, costsToScoreDeltas(weightedSumOfFeatures(extractFeatures(jsonHyp, JsonTypes::Features), weights, linkWeightOffset, statesShareWeights)));
	}

	// read divisions
	for(int i = 0; i < (int)segmentationHypotheses.size(); i++)
	{
		const Json::Value jsonHyp = segmentationHypotheses[i];
		size_t id = jsonHyp[JsonTypeNames[JsonTypes::Id]].asInt();

		if(jsonHyp.isMember(JsonTypeNames[JsonTypes::DivisionFeatures]))
		{
			//graphBuilder_->allowMitosis(id, costsToScoreDelta(weightedSumOfFeatures(extractFeatures(jsonHyp, JsonTypes::DivisionFeatures), weights, divWeightOffset, statesShareWeights)));
		}
	}

	// read exclusion constraints between detections
	// const Json::Value exclusions = root[JsonTypeNames[JsonTypes::Exclusions]];
	// std::cout << "\tcontains " << exclusions.size() << " exclusions" << std::endl;
	// for(int i = 0; i < (int)exclusions.size(); i++)
	// {
	// 	const Json::Value jsonExc = exclusions[i];
	// 	// TODO: implement some API for those?!
	// }
	if(root.isMember(JsonTypeNames[JsonTypes::Exclusions]) && root[JsonTypeNames[JsonTypes::Exclusions]].size() > 0)
		throw std::runtime_error("FlowSolver cannot deal with exclusion constraints yet!");
}

void JsonGraphReader::saveResultJson(const std::string& filename)
{
	std::ofstream output(filename.c_str());
	if(!output.good())
		throw std::runtime_error("Could not open JSON result file for saving: " + filename);

	Json::Value root;

	// save links
	Json::Value& linksJson = root[JsonTypeNames[JsonTypes::LinkResults]];
	//GraphBuilder::ArcValueMap arcValues = graphBuilder_->getArcValues();
  /*
	for(auto iter : arcValues)
	{
		int value = iter.second;

		if(value > 0)
		{
			Json::Value val;
			val[JsonTypeNames[JsonTypes::SrcId]] = Json::Value((unsigned int)iter.first.first);
			val[JsonTypeNames[JsonTypes::DestId]] = Json::Value((unsigned int)iter.first.second);
			val[JsonTypeNames[JsonTypes::Value]] = Json::Value((unsigned int)value);
			linksJson.append(val);
		}
	}
  */

	// save divisions
	Json::Value& divisionsJson = root[JsonTypeNames[JsonTypes::DivisionResults]];
	//GraphBuilder::DivisionValueMap divisionValues = graphBuilder_->getDivisionValues();
  /*
	for(auto iter : divisionValues)
	{
		bool value = iter.second;
		
		if(value > 0)
		{
			Json::Value val;
			val[JsonTypeNames[JsonTypes::Id]] = Json::Value((unsigned int)iter.first);
			val[JsonTypeNames[JsonTypes::Value]] = Json::Value(true);
			divisionsJson.append(val);
		}
	}
  */

	// save detections
	Json::Value& detectionsJson = root[JsonTypeNames[JsonTypes::DetectionResults]];
  /*
	GraphBuilder::NodeValueMap nodeValues = graphBuilder_->getNodeValues();
	for(auto iter : nodeValues)
	{
		size_t value = iter.second;
		
		if(value > 0)
		{
			Json::Value val;
			val[JsonTypeNames[JsonTypes::Id]] = Json::Value((unsigned int)iter.first);
			val[JsonTypeNames[JsonTypes::Value]] = Json::Value((unsigned int)value);
			detectionsJson.append(val);
		}
	}
  */

	output << root << std::endl;
}

JsonGraphReader::StateFeatureVector JsonGraphReader::extractFeatures(const Json::Value& entry, JsonTypes type)
{
	StateFeatureVector stateFeatVec;
	if(!entry.isMember(JsonTypeNames[type]))
		throw std::runtime_error("Could not find Json tags for " + JsonTypeNames[type]);

	const Json::Value featuresPerState = entry[JsonTypeNames[type]];

	if(!featuresPerState.isArray())
		throw std::runtime_error(JsonTypeNames[type] + " must be an array");

	if(!featuresPerState.size() > 0)
		throw std::runtime_error("Features may not be empty for " + JsonTypeNames[type]);

	// std::cout << "\tReading features for: " << JsonTypeNames[type] << std::endl;

	// get the features per state
	for(int i = 0; i < (int)featuresPerState.size(); i++)
	{
		// get features for the specific state
		FeatureVector featVec;
		const Json::Value& featuresForState = featuresPerState[i];

		if(!featuresForState.isArray())
			throw std::runtime_error("Expected to find a list of features for each state");

		if(!featuresForState.size() > 0)
		throw std::runtime_error("Features for state may not be empty for " + JsonTypeNames[type]);

		for(int j = 0; j < (int)featuresForState.size(); j++)
		{
			featVec.push_back(featuresForState[j].asDouble());
		}

		// std::cout << "\t\tfound " << featVec.size() << " features for state " << i << std::endl;

		stateFeatVec.push_back(featVec);
	}

	return stateFeatVec;
}

JsonGraphReader::FeatureVector JsonGraphReader::readWeightsFromJson(const std::string& filename)
{
	std::ifstream input(filename.c_str());
	if(!input.good())
		throw std::runtime_error("Could not open JSON weight file for reading: " + filename);

	Json::Value root;
	input >> root;

	if(!root.isMember(JsonTypeNames[JsonTypes::Weights]))
		throw std::runtime_error("Could not find 'Weights' group in JSON file");
	
	const Json::Value entry = root[JsonTypeNames[JsonTypes::Weights]];
	if(!entry.isArray())
		throw std::runtime_error("Cannot extract Weights from non-array JSON entry");

	FeatureVector weights;
	for(int i = 0; i < (int)entry.size(); i++)
	{
		weights.push_back(entry[i].asDouble());
	}
	return weights;
}
	
} // end namespace dpct
