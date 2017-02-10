#ifndef LP_MP_CONSERVATION_TRACKING_CONSTRUCTOR_HXX
#define LP_MP_CONSERVATION_TRACKING_CONSTRUCTOR_HXX

#include "json_graph_reader.h" // do zrobienia: put all code from Haubold into own directory

namespace LP_MP {

template<typename DETECTION_FACTOR_CONTAINER, typename TRANSITION_MESSAGE_CONTAINER> 
class conservation_tracking_constructor : public dpct::GraphBuilder {
public:
  using CONSTRUCTOR = conservation_tracking_constructor<DETECTION_FACTOR_CONTAINER, TRANSITION_MESSAGE_CONTAINER>;
  template<typename FMC>
  conservation_tracking_constructor(Solver<FMC>& solver) : lp_(&solver.GetLP()) {} 



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
    /*
    INDEX outgoing_edge_index_start = 0;
    for(auto& c : detection_factor_stat_[srcId][1]) { outgoing_edge_index_start += c.size(); }
    detection_factor_stat_[srcId][1].push_back(costDeltas);

    INDEX incoming_edge_index_start = 0;
    for(auto& c : detection_factor_stat_[destId][0]) { incoming_edge_index_start += c.size(); }
    detection_factor_stat_[destId][0].push_back(std::vector<REAL>(costDeltas.size(),0.0));
    */

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

  void construct()
  {
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
   template<typename FMC>
   bool ParseProblem(const std::string& filename, Solver<FMC>& s)
   {
     auto& conservation_tracking_constructor = s.template GetProblemConstructor<0>(); // assume it is the first constructor

     // call Carsten Haubolds json graph reader to parse problem on conservatoin_tracking_constructor. The constructor implements the graphBuilder interface
     std::cout << "reading " << filename << "\n";
     dpct::JsonGraphReader reader(filename, filename, &conservation_tracking_constructor);
     reader.createGraphFromJson();
     conservation_tracking_constructor.construct();
     
     return true;
   }

} // end namespace conservation_tracking_parser

} // end namespace LP_MP

#endif // LP_MP_CONSERVATION_TRACKING_CONSTRUCTOR

