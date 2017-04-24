#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>

#include <vector>

#include "visitors/standard_visitor.hxx"
#include "solvers/multicut.h"

namespace py = pybind11;

namespace LP_MP {
    using FMC = FMC_ODD_WHEEL_MULTICUT<MessageSendingType::SRMP>;
    using SolverType = ProblemConstructorRoundingSolver<Solver<FMC,LP,StandardTighteningVisitor>>;

    // TODO
    struct MulticutOptions {
    //std::vector<std::string> options = {
    //  "-i", " ",
    //  "--primalComputationInterval", "100",
    //  "--standardReparametrization", "anisotropic",
    //  "--roundingReparametrization", "damped_uniform",
    //  "--tightenReparametrization", "damped_uniform",
    //  "--tighten",
    //  "--tightenInterval", "100",
    //  "--tightenIteration", "2",
    //  "--tightenSlope", "0.05",
    //  "--tightenConstraintsPercentage", "0.1"
    //};
    };

    template<typename LABEL_TYPE, typename WEIGHT_TYPE, typename MC_CONSTRUCTOR_TYPE>
    inline void constructMcProblem(
        const std::vector<std::pair<LABEL_TYPE,LABEL_TYPE>> & uvIds,
        const std::vector<WEIGHT_TYPE> & weights,
        MC_CONSTRUCTOR_TYPE & constructor
    ) {
        for(size_t i = 0; i < uvIds.size(); ++i) {
            const auto & uv = uvIds[i];
            constructor.AddUnaryFactor(uv.first, uv.second, weights[i]);
        }
    }
    
    template<typename LABEL_TYPE, typename MC_CONSTRUCTOR_TYPE>
    inline void getEdgeLabeling(
        const std::vector<std::pair<LABEL_TYPE,LABEL_TYPE>> & uvIds,
        const MC_CONSTRUCTOR_TYPE & constructor,
        std::vector<bool> & edgeLabeling
    ) {
        for(size_t i = 0; i < uvIds.size(); ++i) {
            const auto & uv = uvIds[i];
            edgeLabeling[i] = constructor.get_edge_label(uv.first, uv.second);
        }
    }

    template<typename LABEL_TYPE, typename WEIGHT_TYPE>
    void exportMulticutT(py::module & pyModule) {
        
        // TODO nested vector is ugly, use proper multiarray instead
        pyModule.def("multicut", 
            [](
                const std::vector<std::pair<LabelType,LabelType>> & uvIds,
                const std::vector<WeightType> & weights,
                const MulticutOptions & mcOpts
            ) {
                
                {
                    py::gil_scoped_release allowThreads;
                    SolverType solver(mcOpts.toString());
                    auto & multicut_constructor = solver.template GetProblemConstructor<0>();
                    std::vector<bool> edgeLabels(uvIds.size());
                }
                return edgeLabels;
            },
            py::arg("uvIds"), py::arg("weights")
        );
    }

    void exportMulticut(py::module & pyModule) {
        exportMulticutT<uint32_t, float>(pyModule);
    }
}
