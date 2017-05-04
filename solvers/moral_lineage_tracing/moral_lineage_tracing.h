#ifndef LP_MP_MORAL_LINEAGE_TRACING_H
#define LP_MP_MORAL_LINEAGE_TRACING_H

#include "LP_MP.h"
#include "factors_messages.hxx"

#include "solvers/multicut/multicut_factors_messages.hxx"
#include "solvers/multicut/multicut_constructor.hxx"
#include "solvers/multicut/lifted_multicut_factors_messages.hxx"
#include "factors/constant_factor.hxx"
#include "solvers/multicut/multicut_constructor.hxx"
#include "solvers/multicut/multicut.h"

#include "mlt_factors_messages.hxx"
#include "mlt_constructor.hxx"

#include "andres/functional.hxx"

#include <vector>
#include <fstream>

namespace LP_MP {
	
struct FMC_MLT {
	constexpr static char* name = "moral lineage tracing";

	using mc_edge_factor_container = FactorContainer<multicut_edge_factor, FMC_MLT, 0>;
	using mc_triplet_factor_container = FactorContainer<multicut_triplet_factor, FMC_MLT, 1>;
	using mlt_triplet_factor_container = FactorContainer<mlt_triplet_factor, FMC_MLT, 2>;
	using mlt_bifurcation_factor_container = FactorContainer<mlt_bifurcation_factor, FMC_MLT, 3>;
	using mc_lifted_cut_factor_container = FactorContainer<LiftedMulticutCutFactor, FMC_MLT, 4>;

	using constant_factor_container = FactorContainer<ConstantFactor, FMC_MLT, 5>;

	using mc_edge_triplet_message_0_container = MessageContainer<multicut_edge_triplet_message_0, 0, 1, variableMessageNumber, 1, FMC_MLT, 0>;
	using mc_edge_triplet_message_1_container = MessageContainer<multicut_edge_triplet_message_1, 0, 1, variableMessageNumber, 1, FMC_MLT, 1>;
	using mc_edge_triplet_message_2_container = MessageContainer<multicut_edge_triplet_message_2, 0, 1, variableMessageNumber, 1, FMC_MLT, 2>;

	using edge_mlt_triplet_message_0_container = MessageContainer<edge_mlt_triplet_message_0, 0, 2, variableMessageNumber, 1, FMC_MLT, 3>;
	using edge_mlt_triplet_message_1_container = MessageContainer<edge_mlt_triplet_message_1, 0, 2, variableMessageNumber, 1, FMC_MLT, 4>;
	using edge_mlt_triplet_message_2_container = MessageContainer<edge_mlt_triplet_message_2, 0, 2, variableMessageNumber, 1, FMC_MLT, 5>;

	using mlt_triplet_bifurcation_message_0_container = MessageContainer<mlt_triplet_bifurcation_message_0, 2, 3, variableMessageNumber, 1, FMC_MLT, 6>;
	using mlt_triplet_bifurcation_message_1_container = MessageContainer<mlt_triplet_bifurcation_message_1, 2, 3, variableMessageNumber, 1, FMC_MLT, 7>;
	using mlt_triplet_bifurcation_message_2_container = MessageContainer<mlt_triplet_bifurcation_message_2, 2, 3, variableMessageNumber, 1, FMC_MLT, 8>;

	using mc_triplet_bifurcation_message_container = MessageContainer<mc_triplet_bifurcation_message, 1, 3, variableMessageNumber, 1, FMC_MLT, 9>;

	using mc_base_edge_lifted_cut_factor_message_container = MessageContainer<CutEdgeLiftedMulticutFactorMessage, 0, 4, variableMessageNumber, variableMessageNumber, FMC_MLT, 10>;
	using mc_lifted_edge_lifted_cut_factor_message_container = MessageContainer<LiftedEdgeLiftedMulticutFactorMessage, 0, 4, variableMessageNumber, variableMessageNumber, FMC_MLT, 11>;


	using FactorList = meta::list<
	mc_edge_factor_container, 
	mc_triplet_factor_container, 
	mlt_triplet_factor_container,
	mlt_bifurcation_factor_container,
	mc_lifted_cut_factor_container,
	constant_factor_container
	>;

	using MessageList = meta::list<
	mc_edge_triplet_message_0_container, mc_edge_triplet_message_1_container, mc_edge_triplet_message_2_container,
	edge_mlt_triplet_message_0_container, edge_mlt_triplet_message_1_container, edge_mlt_triplet_message_2_container,
	mlt_triplet_bifurcation_message_0_container, mlt_triplet_bifurcation_message_1_container, mlt_triplet_bifurcation_message_2_container,
	mc_triplet_bifurcation_message_container,
	mc_base_edge_lifted_cut_factor_message_container, mc_lifted_edge_lifted_cut_factor_message_container
	>;

	// input FMC struct, factor numbers, message numbers, constant factor
	using mc_constructor = MulticutConstructor<FMC_MLT, 0, 1, 0, 1, 2, 5>;
	using lifted_mc_constructor = LiftedMulticutConstructor<mc_constructor, 4, 10, 11>;
	using mlt_constructor = MLT_constructor<lifted_mc_constructor, 2, 3, 3, 4, 5, 6, 7, 8, 9>;
	using ProblemDecompositionList = meta::list<mlt_constructor>;
};

namespace mlt_input {

struct Node
{
    size_t t;
    size_t id;
    size_t cx, cy;
    double probability_birth_termination;
};

struct Edge
{
    size_t t0, v0;
    size_t t1, v1;
    double weight;
};

inline size_t NodeIndex(size_t node_num, size_t frame_num, std::vector<size_t>& num_of_nodes_per_frame)
{
	return num_of_nodes_per_frame[frame_num] + node_num;
}

template<typename SOLVER>
bool ParseProblem(const std::string filename, SOLVER& s)
{
	auto& mlt_constr = s.template GetProblemConstructor<0>();
	const auto& node_file_name = mlt_constr.get_node_file();
	const auto& edge_file_name = mlt_constr.get_edge_file();

	std::vector<Node> nodes;
	std::vector<Edge> edges;
	std::vector<size_t> num_of_nodes_per_frame;

	Node node;
	std::ifstream node_file(node_file_name);
	size_t counter = 0;
    while (node_file >> node.t >> node.id >> node.cx >> node.cy >> node.probability_birth_termination)
    {
    	if (node.t == num_of_nodes_per_frame.size())
    		num_of_nodes_per_frame.push_back(counter);

        nodes.push_back(node);
        counter++;
    }
    node_file.close();

    mlt_constr.SetNumOfNodes(counter);

    Edge edge;
    std::ifstream edge_file(edge_file_name);
    while (edge_file >> edge.t0 >> edge.v0 >> edge.t1 >> edge.v1 >> edge.weight)
    {
        edges.push_back(edge);
    }
	edge_file.close();

	// add node frame numbers to MLT constructor
	for (Node node : nodes)
	{
		mlt_constr.AddNode(node.t, NodeIndex(node.id, node.t, num_of_nodes_per_frame));
	}

	// add edge factors
	andres::NegativeLogProbabilityRatio<double,double> func;
	size_t node_0;
	size_t node_1;
	for (Edge edge : edges)
	{
		node_0 = NodeIndex(edge.v0,edge.t0,num_of_nodes_per_frame);
		node_1 = NodeIndex(edge.v1,edge.t1,num_of_nodes_per_frame);
		if (node_0 < node_1)
			mlt_constr.AddBaseEdge(node_0,node_1,func(edge.weight));
		else
			mlt_constr.AddBaseEdge(node_1,node_0,func(edge.weight));
	}
}

} // end namespace mlt_input

} // end namespace LP_MP

#endif
