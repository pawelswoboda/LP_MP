#ifndef LP_MP_MLT_CONSTRUCTOR_HXX
#define LP_MP_MLT_CONSTRUCTOR_HXX

// #include "LP_MP.h"

#include "solvers/multicut/multicut_constructor.hxx"

#include <map>
#include <array>

#include "tclap/CmdLine.h"

namespace LP_MP {

template<typename LIFTED_MULTICUT_CONSTRUCTOR,
INDEX MLT_TRIPLET_FACTOR_NO, INDEX MLT_BIFURCATION_FACTOR_NO,
INDEX EDGE_MLT_TRIPLET_MESSAGE_0_NO, INDEX EDGE_MLT_TRIPLET_MESSAGE_1_NO, INDEX EDGE_MLT_TRIPLET_MESSAGE_2_NO, 
INDEX MLT_TRIPLET_BIFURCATION_MESSAGE_0_NO, INDEX MLT_TRIPLET_BIFURCATION_MESSAGE_1_NO, INDEX MLT_TRIPLET_BIFURCATION_MESSAGE_2_NO,
INDEX MC_TRIPLET_BIFURCATION_MESSAGE_NO
>
class MLT_constructor {
public:
	using FMC = typename LIFTED_MULTICUT_CONSTRUCTOR::FMC;

	using mlt_triplet_factor_container = meta::at_c<typename FMC::FactorList, MLT_TRIPLET_FACTOR_NO>;
	using mlt_bifurcation_factor_container = meta::at_c<typename FMC::FactorList, MLT_BIFURCATION_FACTOR_NO>;
	
	using edge_mlt_triplet_message_0_container = meta::at_c<typename FMC::MessageList, EDGE_MLT_TRIPLET_MESSAGE_0_NO>;
	using edge_mlt_triplet_message_1_container = meta::at_c<typename FMC::MessageList, EDGE_MLT_TRIPLET_MESSAGE_1_NO>;
	using edge_mlt_triplet_message_2_container = meta::at_c<typename FMC::MessageList, EDGE_MLT_TRIPLET_MESSAGE_2_NO>;

	using mlt_triplet_bifurcation_message_0_container = meta::at_c<typename FMC::MessageList, MLT_TRIPLET_BIFURCATION_MESSAGE_0_NO>;
	using mlt_triplet_bifurcation_message_1_container = meta::at_c<typename FMC::MessageList, MLT_TRIPLET_BIFURCATION_MESSAGE_1_NO>;
	using mlt_triplet_bifurcation_message_2_container = meta::at_c<typename FMC::MessageList, MLT_TRIPLET_BIFURCATION_MESSAGE_2_NO>;

	using mc_triplet_bifurcation_message_container = meta::at_c<typename FMC::MessageList, MC_TRIPLET_BIFURCATION_MESSAGE_NO>;

	template<typename SOLVER>
    MLT_constructor(SOLVER& s) 
    : lp_(&s.GetLP()),
    lifted_mc_constructor_(s),
    termination_cost_arg_("T", "terminationCost", "cost for node termination", true, 5.0, "real value", s.get_cmd()),
    birth_cost_arg_("B", "birthCost", "cost for node birth", true, 5.0, "real value", s.get_cmd()),
    node_file_arg_("n", "nodeFile", "input file for the nodes", true, "", "file name", s.get_cmd()),
    edge_file_arg_("e", "edgeFile", "input file for the edges", true, "", "file name", s.get_cmd()),
    node_to_timestep_(std::vector<INDEX>())
	{}

	const std::string& get_node_file()
	{
		return node_file_arg_.getValue();
	}

	const std::string& get_edge_file()
	{
		return edge_file_arg_.getValue();
	}

	void SetNumOfNodes(const INDEX num_of_nodes)
	{
		node_to_timestep_.resize(num_of_nodes);
	}

	void AddNode(const INDEX timestep, const INDEX node)
	{
		node_to_timestep_[node] = timestep;
	}

	void AddBaseEdge(const INDEX node_1, const INDEX node_2, const REAL cost)
	{
		auto* f = lifted_mc_constructor_.AddUnaryFactor(node_1, node_2, cost);
	}

	void AddLiftedEdge(const INDEX node_1, const INDEX node_2, const REAL cost)
	{
		auto* f = lifted_mc_constructor_.AddLiftedUnaryFactor(node_1, node_2, cost);
	}

	bool HasEdge(const INDEX node_1, const INDEX node_2) 
	{
		return lifted_mc_constructor_.HasUnaryFactor(node_1, node_2);
	}

	void AddMcTriplet(const INDEX node_1, const INDEX node_2, const INDEX node_3)
	{
		auto* f = lifted_mc_constructor_.AddTripletFactor(node_1, node_2, node_3);
	}

	bool HasMcTriplet(const INDEX node_1, const INDEX node_2, const INDEX node_3)
	{
		return lifted_mc_constructor_.HasTripletFactor(node_1, node_2, node_3);
	}

	void AddMltTriplet(const INDEX node_1, const INDEX node_2, const INDEX node_3)
	{
		assert(node_1 < node_2 && node_2 < node_3);
		assert(!HasMltTriplet(node_1, node_2, node3));

		auto* f = new mlt_triplet_factor_container();
		lp_->AddFactor(f);
		mlt_triplet_factors_.emplace(std::array<INDEX,3>{node_1, node_2, node_3}, f);

		auto* e1 = lifted_mc_constructor_.GetUnaryFactor(node_1, node_2);
		auto* m1 = new edge_mlt_triplet_message_0_container(e1, f);
		lp_->AddMessage(m1);

		auto* e2 = lifted_mc_constructor_.GetUnaryFactor(node_1, node_3);
		auto* m2 = new edge_mlt_triplet_message_1_container(e2, f);
		lp_->AddMessage(m2);		

		auto* e3 = lifted_mc_constructor_.GetUnaryFactor(node_2, node_3);
		auto* m3 = new edge_mlt_triplet_message_2_container(e3, f);
		lp_->AddMessage(m3);
	}

	bool HasMltTriplet(const INDEX node_1, const INDEX node_2, const INDEX node_3)
	{
		assert(node_1 < node_2 && node_2 < node_3);
		return mlt_triplet_factors_.find(std::array<INDEX,3>{node_1, node_2, node_3}) != mlt_triplet_factors_.end();
	}

	bool HasTriplet(const INDEX node_1, const INDEX node_2, const INDEX node_3)
	{
		return HasMcTriplet(node_1, node_2, node_3) || HasMltTriplet(node_1, node_2, node_3);
	}

	void AddBifurcation(const INDEX node_1, const INDEX node_2, const INDEX node_3, const INDEX node_4, const INDEX node_5, const INDEX node_6)
	{
		assert(node_1 < node_2 && node_2 < node_3 && node_3 < node_4 && node_4 < node_5 && node_5 < node_6);
		assert(!HasBifurcation(node_1, node_2, node_3, node_4, node_5, node_6));

		auto* f = new mlt_bifurcation_factor_container();
		lp_->AddFactor(f);
		mlt_bifurcation_factors_.emplace(std::array<INDEX,6>{node_1,node_2,node_3,node_4,node_5,node_6}, f);

		auto* t1 = mlt_triplet_factors_.find(std::array<INDEX,3>{node_1, node_2, node_3})->second;
		auto* m1 = new mlt_triplet_bifurcation_message_0_container(t1, f);
		lp_->AddMessage(m1);

		auto* t2 = mlt_triplet_factors_.find(std::array<INDEX,3>{node_2, node_3, node_5})->second;
		auto* m2 = new mlt_triplet_bifurcation_message_1_container(t2, f);
		lp_->AddMessage(m2);

		auto* t3 = mlt_triplet_factors_.find(std::array<INDEX,3>{node_1, node_3, node_6})->second;
		auto* m3 = new mlt_triplet_bifurcation_message_2_container(t3, f);
		lp_->AddMessage(m3);

		auto* t4 = lifted_mc_constructor_.GetTripletFactor(node_4, node_5, node_6);
		auto* m4 = new mc_triplet_bifurcation_message_container(t4, f);
		lp_->AddMessage(m4);
	}

	bool HasBifurcation(const INDEX node_1, const INDEX node_2, const INDEX node_3, const INDEX node_4, const INDEX node_5, const INDEX node_6)
	{
		assert(node_1 < node_2 && node_2 < node_3 && node_3 < node_4 && node_4 < node_5 && node_5 < node_6);
		return mlt_bifurcation_factors_.find(std::array<INDEX,6>{node_1, node_2, node_3, node_4, node_5, node_6}) != mlt_bifurcation_factors_.end();
	}

	INDEX SearchViolatedMulticutCycles(const INDEX maxTripletsToAdd)
	{
		std::vector<INDEX> number_outgoing_arcs(lifted_mc_constructor_.num_nodes(),0); // number of arcs outgoing arcs of each node
		std::vector<std::tuple<INDEX,INDEX,REAL>> negative_edges;
		INDEX number_arcs_total = 0;
		for(INDEX e=0; e<lifted_mc_constructor_.NumOfUnaryFactors(); ++e) {
			const REAL v = lifted_mc_constructor_.GetEdgeCost(e);
			const INDEX i = std::get<0>(lifted_mc_constructor_.GetEdge(e));
			const INDEX j = std::get<1>(lifted_mc_constructor_.GetEdge(e));

			if(v > 0) {
				number_outgoing_arcs[i]++;
				number_outgoing_arcs[j]++;
			} 
			else if(v < 0) {
				negative_edges.push_back(std::make_tuple(i,j,v));
			}
			number_arcs_total += 2;
		}

		Graph g(number_outgoing_arcs.size(), number_arcs_total, number_outgoing_arcs); // graph consisting of positive edges
		for(INDEX e=0; e<lifted_mc_constructor_.NumOfUnaryFactors(); ++e) {
			const REAL v = lifted_mc_constructor_.GetEdgeCost(e);
			const INDEX i = std::get<0>(lifted_mc_constructor_.GetEdge(e));
			const INDEX j = std::get<1>(lifted_mc_constructor_.GetEdge(e));

			if(v > 0) {
				g.add_edge(i,j,v);
			}
		}      
		g.sort();

		std::sort(negative_edges.begin(), negative_edges.end(), [](auto& a, auto& b) { return std::get<2>(a) < std::get<2>(b); });

		BfsData bfs(g);

		using CycleType = std::tuple<REAL, std::vector<INDEX>>;
		std::vector<CycleType> found_cycles;

		for(auto& it : negative_edges) {
			const INDEX i = std::get<0>(it);
			const INDEX j = std::get<1>(it);
			const REAL v = std::get<2>(it);
			const INDEX timestep = std::min(node_to_timestep_[i], node_to_timestep_[j]);
			// subgraph mask for frame t and t+1
			auto mask_op = [this, timestep] (const INDEX tail, const INDEX head, const REAL weight) {
				return node_to_timestep_[tail] >= timestep && node_to_timestep_[tail] <= timestep + 1 &&
					node_to_timestep_[head] >= timestep && node_to_timestep_[head] <= timestep + 1;
			};

			auto cycle = bfs.FindPath(i,j,g,0.0,mask_op);

			assert(std::get<1>(cycle).size() > 0);
			if(std::get<1>(cycle).size() > 0) {
				const REAL dualIncrease = std::min(-v, std::get<0>(cycle));
				found_cycles.push_back( std::make_tuple(dualIncrease, std::move(std::get<1>(cycle))) );
			}
		}

		// sort by guaranteed increase in decreasing order
		std::sort(found_cycles.begin(), found_cycles.end(), [](const CycleType& i, const CycleType& j) { return std::get<0>(i) > std::get<0>(j); });
		INDEX tripletsAdded = 0;
		for (auto& cycle : found_cycles) {
			if (std::get<1>(cycle).size() > 2) {
				tripletsAdded += AddCycle(std::move(std::get<1>(cycle)));
				if (tripletsAdded > maxTripletsToAdd) {
					return tripletsAdded;
				}
			}
		}

	}

	INDEX AddCycle(std::vector<INDEX> cycle)
	{
		assert(cycle.size() >= 3);
		std::rotate(cycle.begin(), std::min_element(cycle.begin(), cycle.end()), cycle.end());
		const INDEX firstNode = cycle[0];
		const INDEX t = node_to_timestep_[firstNode];

		// triangulate by adding edges from firstNode to all other nodes
		INDEX noTripletsAdded = 0;
		for (INDEX i = 2; i < cycle.size(); ++i)
		{
			if (!HasEdge(firstNode, cycle[i]))
// TODO: shouldn't we add a lifted edge?
				AddBaseEdge(firstNode, cycle[i], 0.0);

			const INDEX secondNode = std::min(cycle[i], cycle[i-1]);
			const INDEX thirdNode = std::max(cycle[i], cycle[i-1]);
			const INDEX t1 = node_to_timestep_[secondNode];
			const INDEX t2 = node_to_timestep_[thirdNode];
			if (!HasTriplet(firstNode, secondNode, thirdNode))
			{
				// firstNode has the minimal index, thus minimal timestep among the triplet
				// if the other timesteps are both higher, then we have an MLT triplet
				if (t < t1 && t < t2)
					AddMltTriplet(firstNode, secondNode, thirdNode);
				else
					AddMcTriplet(firstNode, secondNode, thirdNode);
				noTripletsAdded++;
			}
		}
		return noTripletsAdded;
	}

   	INDEX Tighten(const INDEX maxCuttingPlanesToAdd)
	{}

   	void ComputePrimal()
	{
	}

private:
	LP* lp_;
	LIFTED_MULTICUT_CONSTRUCTOR lifted_mc_constructor_;
	std::map<std::array<INDEX,3>, mlt_triplet_factor_container*> mlt_triplet_factors_;
	std::map<std::array<INDEX,6>, mlt_bifurcation_factor_container*> mlt_bifurcation_factors_;
	std::vector<INDEX> node_to_timestep_;

	TCLAP::ValueArg<REAL> termination_cost_arg_;
	TCLAP::ValueArg<REAL> birth_cost_arg_;
	TCLAP::ValueArg<std::string> node_file_arg_;
	TCLAP::ValueArg<std::string> edge_file_arg_;
};
}

#endif 
