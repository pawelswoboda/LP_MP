#ifndef LP_MP_CONSERVATION_TRACKING_CONSTRUCTOR_HXX
#define LP_MP_CONSERVATION_TRACKING_CONSTRUCTOR_HXX

#include "json_graph_reader.h" // do zrobienia: put all code from Haubold into own directory

namespace LP_MP {

template<typename DETECTION_FACTOR_CONTAINER, typename TRANSITION_MESSAGE_CONTAINER> 
class conservation_tracking_constructor {
public:
  using CONSTRUCTOR = conservation_tracking_constructor<DETECTION_FACTOR_CONTAINER, TRANSITION_MESSAGE_CONTAINER>;
  template<typename SOLVER>
  conservation_tracking_constructor(SOLVER& solver) : lp_(&solver.GetLP()) {} 

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


  /*

  virtual void addNode(
      size_t id,
      const CostDeltaVector& detectionCosts,
      const CostDeltaVector& detectionCostDeltas, 
      const CostDeltaVector& appearanceCostDeltas, 
      const CostDeltaVector& disappearanceCostDeltas,
      size_t targetIdx)
  {
    assert(targetIdx == 0);
    assert(detectionCosts.size() == detectionCostDeltas.size()+1);
    assert(appearanceCostDeltas.size() == detectionCostDeltas.size());
    assert(disappearanceCostDeltas.size() == detectionCostDeltas.size());

    std::cout << "add node with id = " << id << " and target idx = " << targetIdx << "\n";
    assert(id >= detection_factor_stat_.size());
    detection_factor_stat_.resize(id+1,{0,0});

    detection_costs_.resize(id+1);
    detection_costs_[id] = detectionCosts;

    appearance_costs_.resize(id+1);
    appearance_costs_[id] = appearanceCostDeltas;

    disappearance_costs_.resize(id+1);
    disappearance_costs_[id] = disappearanceCostDeltas;
  }

  virtual void addArc(size_t srcId, size_t destId, const CostDeltaVector& costDeltas)
  {
    std::cout << "add arc (" << srcId << "," << destId << ")\n";
    detection_factor_stat_[srcId][1] += costDeltas.size();
    detection_factor_stat_[destId][0] += costDeltas.size();

    //INDEX outgoing_edge_index_start = 0;
    //for(auto& c : detection_factor_stat_[srcId][1]) { outgoing_edge_index_start += c.size(); }
    //detection_factor_stat_[srcId][1].push_back(costDeltas);

    //INDEX incoming_edge_index_start = 0;
    //for(auto& c : detection_factor_stat_[destId][0]) { incoming_edge_index_start += c.size(); }
    //detection_factor_stat_[destId][0].push_back(std::vector<REAL>(costDeltas.size(),0.0));

    for(INDEX i=0; i<costDeltas.size(); ++i) {
      arcs_.push_back({INDEX(srcId), INDEX(destId), costDeltas[i]});
    }
  }

  virtual void allowMitosis(size_t id, ValueType divisionCostDelta) 
  {
    assert(false);
    // do zrobienia: for each edge pair, introduce additional two edges with the division cost
  }

  virtual NodeValueMap getNodeValues() 
  {}

  virtual ArcValueMap getArcValues()
  {}

  virtual DivisionValueMap getDivisionValues() 
  {}
  */

  size_t getNumWeights(const Json::Value& jsonHyp, JsonTypes type, bool statesShareWeights)
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

  dpct::JsonGraphReader::FeatureVector readWeightsFromJson(const std::string& filename)
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

  dpct::JsonGraphReader::StateFeatureVector extractFeatures(const Json::Value& entry, JsonTypes type)
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

  void construct(const std::string filename)
  {
    detection_factor_stat_.clear(); // no incoming edges, no outgoing edges

     std::ifstream input(filename.c_str());
     if(!input.good())
       throw std::runtime_error("Could not open JSON model file for reading: " + filename);

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
     auto weights = readWeightsFromJson(filename);
     size_t numDetWeights = 0;
     size_t numDivWeights = 0;
     size_t numAppWeights = 0;
     size_t numDisWeights = 0;
     size_t numLinkWeights = 0;

     const Json::Value segmentationHypotheses = root[JsonTypeNames[JsonTypes::Segmentations]];
     detection_factor_stat_.resize( segmentationHypotheses.size() );
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
     // read segmentation hypotheses 
     std::cout << "\tcontains " << segmentationHypotheses.size() << " segmentation hypotheses" << std::endl;

     for(int i = 0; i < (int)segmentationHypotheses.size(); i++)
     {
       const Json::Value jsonHyp = segmentationHypotheses[i];

       if(!jsonHyp.isMember(JsonTypeNames[JsonTypes::Id]))
         throw std::runtime_error("Cannot read detection hypothesis without Id!");
       size_t id = jsonHyp[JsonTypeNames[JsonTypes::Id]].asInt();
       assert(id == i);

       if(!jsonHyp.isMember(JsonTypeNames[JsonTypes::Features]))
         throw std::runtime_error("Cannot read detection hypothesis without features!");

       if(jsonHyp.isMember(JsonTypeNames[JsonTypes::Timestep]))
       {
         const Json::Value& timeJson = jsonHyp[JsonTypeNames[JsonTypes::Timestep]];
         if(!timeJson.isArray() || timeJson.size() != 2)
           throw std::runtime_error("Node's Timestep is supposed to be a 2-element array");
         std::pair<int, int> timeRange = std::make_pair(timeJson[0].asInt(), timeJson[1].asInt());
         //graphBuilder_->setNodeTimesteps(id, timeRange);
       }

       auto detCosts = weightedSumOfFeatures(extractFeatures(jsonHyp, JsonTypes::Features), weights, detWeightOffset, statesShareWeights);
       auto detCostDeltas = costsToScoreDeltas(detCosts);
       FeatureVector appearanceCost;
       if(jsonHyp.isMember(JsonTypeNames[JsonTypes::AppearanceFeatures])) {
         appearanceCost = weightedSumOfFeatures(extractFeatures(jsonHyp, JsonTypes::AppearanceFeatures), weights, appWeightOffset, statesShareWeights);
       }

       FeatureVector disappearanceCost;
       if(jsonHyp.isMember(JsonTypeNames[JsonTypes::DisappearanceFeatures])) {
         disappearanceCost = weightedSumOfFeatures(extractFeatures(jsonHyp, JsonTypes::DisappearanceFeatures), weights, disWeightOffset, statesShareWeights);
       }

       size_t targetIdx = 0;
       if(jsonHyp.isMember(JsonTypeNames[JsonTypes::DisappearanceTarget]))
       {
         targetIdx = jsonHyp[JsonTypeNames[JsonTypes::DisappearanceTarget]].asUInt();
       }

       //graphBuilder_->addNode(id, detCosts, detCostDeltas, appearanceCostDeltas, disappearanceCostDeltas, targetIdx);
     }

     // read linking hypotheses
     std::cout << "\tcontains " << linkingHypotheses.size() << " linking hypotheses" << std::endl;
     for(int i = 0; i < (int)linkingHypotheses.size(); i++)
     {
       const Json::Value jsonHyp = linkingHypotheses[i];

       const size_t srcId = jsonHyp[JsonTypeNames[JsonTypes::SrcId]].asInt();
       const size_t destId = jsonHyp[JsonTypeNames[JsonTypes::DestId]].asInt();
       // transform hypothesis to detection factor number
       auto linkingCosts = weightedSumOfFeatures(extractFeatures(jsonHyp, JsonTypes::Features), weights, linkWeightOffset, statesShareWeights);
       detection_factor_stat_[destId][0] += linkingCosts.size();
       detection_factor_stat_[srcId][1] += linkingCosts.size();
       //graphBuilder_->addArc(srcId, destId, costsToScoreDeltas(weightedSumOfFeatures(extractFeatures(jsonHyp, JsonTypes::Features), weights, linkWeightOffset, statesShareWeights)));
     }

     // read divisions
     for(int i = 0; i < (int)segmentationHypotheses.size(); i++)
     {
       const Json::Value jsonHyp = segmentationHypotheses[i];
       size_t id = jsonHyp[JsonTypeNames[JsonTypes::Id]].asInt();

       if(jsonHyp.isMember(JsonTypeNames[JsonTypes::DivisionFeatures]))
       {
         auto divisionCosts = weightedSumOfFeatures(extractFeatures(jsonHyp, JsonTypes::DivisionFeatures), weights, divWeightOffset, statesShareWeights);
         //graphBuilder_->allowMitosis(id, costsToScoreDelta(weightedSumOfFeatures(extractFeatures(jsonHyp, JsonTypes::DivisionFeatures), weights, divWeightOffset, statesShareWeights)));
       }
     }

     // read exclusion constraints between detections
     // const Json::Value exclusions = root[JsonTypeNames[JsonTypes::Exclusions]];
     // std::cout << "\tcontains " << exclusions.size() << " exclusions" << std::endl;
     // for(int i = 0; i < (int)exclusions.size(); i++)
     // {
     //  const Json::Value jsonExc = exclusions[i];
     //  // TODO: implement some API for those?!
     // }
     if(root.isMember(JsonTypeNames[JsonTypes::Exclusions]) && root[JsonTypeNames[JsonTypes::Exclusions]].size() > 0)
       throw std::runtime_error("FlowSolver cannot deal with exclusion constraints yet!");


     // now actually construct factors
     assert(detection_factor_stat_.size() == detection_costs_.size());
     detections_.resize(detection_costs_.size(), nullptr);
     for(INDEX i=0; i<detection_costs_.size(); ++i) {
       if(detection_costs_[i].size() > 0) {
         auto* f = new DETECTION_FACTOR_CONTAINER(detection_costs_[i].size(), detection_factor_stat_[0].size(), detection_factor_stat_[1].size());
         for(INDEX j=0; j<detection_costs_[i].size(); ++j) {
           (*f->GetFactor())[j] = detection_costs_[i][j];
         }

         lp_->AddFactor(f);
         detections_[i] = f;
       }
     }
     /*

    std::vector<std::array<INDEX,2>> current_arc_no(detection_factor_stat_.size(),{0,0});
    for(auto a : arcs_) {
      INDEX& outgoing_edge_index = current_arc_no[a.tail][1];
      INDEX& incoming_edge_index = current_arc_no[a.head][0];
      auto* m = new TRANSITION_MESSAGE_CONTAINER(detections_[a.tail], detections_[a.head], false, outgoing_edge_index, incoming_edge_index); 
      lp_->AddMessage(m);

      detections_[a.tail]->GetFactor()->outgoing(outgoing_edge_index) = a.cost;

      ++outgoing_edge_index; assert(current_arc_no[a.tail][1] == outgoing_edge_index);
      ++incoming_edge_index; assert(current_arc_no[a.head][0] == incoming_edge_index);
    }

    // clear temporary data from reading in
    detection_factor_stat_.clear(); detection_factor_stat_.shrink_to_fit();
    detection_costs_.clear(); detection_costs_.shrink_to_fit();
    appearance_costs_.clear(); appearance_costs_.shrink_to_fit();
    disappearance_costs_.clear(); disappearance_costs_.shrink_to_fit();
    arcs_.clear(); arcs_.shrink_to_fit();
    */

  }

private:
  LP* lp_;
  std::vector<DETECTION_FACTOR_CONTAINER*> detections_;

  // temporary data for holding costs before building factors
  //std::vector<std::array<std::vector<std::vector<REAL>>,2>> detection_factor_stat_;
  std::vector<std::array<INDEX,2>> detection_factor_stat_; // no incoming edges, no outgoing edges
  std::vector<std::vector<REAL>> detection_costs_;
  std::vector<std::vector<REAL>> appearance_costs_;
  std::vector<std::vector<REAL>> disappearance_costs_;
  struct arc {INDEX tail, head; REAL cost;};
  std::vector<arc> arcs_; 
};

namespace conservation_tracking_parser {

   // we assume problem with single constructor
   template<typename SOLVER>
   bool ParseProblem(const std::string& filename, SOLVER& s)
   {
     auto& conservation_tracking_constructor = s.template GetProblemConstructor<0>(); // assume it is the first constructor

     // call Carsten Haubolds json graph reader to parse problem on conservation_tracking_constructor. The constructor implements the graphBuilder interface
     std::cout << "reading " << filename << "\n";
     //dpct::JsonGraphReader reader(filename, filename, &conservation_tracking_constructor);




     //reader.createGraphFromJson();
     //conservation_tracking_constructor.construct();

     return true;
   }

} // end namespace conservation_tracking_parser

} // end namespace LP_MP

#endif // LP_MP_CONSERVATION_TRACKING_CONSTRUCTOR

