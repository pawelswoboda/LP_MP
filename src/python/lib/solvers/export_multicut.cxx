#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>

#include <vector>

#include "visitors/standard_visitor.hxx"
#include "solvers/multicut/multicut.h"

namespace py = pybind11;

namespace LP_MP {
    using FMC = FMC_ODD_WHEEL_MULTICUT<MessageSendingType::SRMP>;
    using SolverType = ProblemConstructorRoundingSolver<Solver<FMC,LP,StandardTighteningVisitor>>;

    struct MulticutOptions {
        MulticutOptions(
            const size_t primalComputationInterval = 100,
            const std::string & standardReparametrization = "anisotropic",
            const std::string & roundingReparametrization = "damped_uniform",
            const std::string & tightenReparametrization  = "damped_uniform",
            const bool tighten = true,
            const size_t tightenInterval = 100,
            const size_t tightenIteration = 2,
            const double tightenSlope = 0.05,
            const double tightenConstraintsPercentage = 0.1
        ) :
            primalComputationInterval_(primalComputationInterval),
            standardReparametrization_(standardReparametrization),
            roundingReparametrization_(roundingReparametrization),
            tightenReparametrization_(tightenReparametrization),
            tighten_(tighten),
            tightenInterval_(tightenInterval),
            tightenIteration_(tightenIteration),
            tightenSlope_(tightenSlope),
            tightenConstraintsPercentage_(tightenConstraintsPercentage)
        {}

        std::vector<std::string> toOptionsVector() const {
            std::vector<std::string> options = {
              "-i", " ", // empty input file
              "--primalComputationInterval", std::to_string(primalComputationInterval_),
              "--standardReparametrization", standardReparametrization_,
              "--roundingReparametrization", roundingReparametrization_,
              "--tightenReparametrization",  tightenReparametrization_,
              "--tightenInterval",           std::to_string(tightenInterval_),
              "--tightenIteration",          std::to_string(tightenIteration_),
              "--tightenSlope",              std::to_string(tightenSlope_),
              "--tightenConstraintsPercentage", std::to_string(tightenConstraintsPercentage_)
            };
            if(tighten_)
                options.push_back("--tighten");
            return options;
        }
        
    private:
        const size_t primalComputationInterval_;
        const std::string standardReparametrization_;
        const std::string roundingReparametrization_;
        const std::string tightenReparametrization_;
        const bool tighten_;
        const size_t tightenInterval_;
        const size_t tightenIteration_;
        const double tightenSlope_;
        const double tightenConstraintsPercentage_;
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
        edgeLabeling.resize(uvIds.size());
        for(size_t i = 0; i < uvIds.size(); ++i) {
            const auto & uv = uvIds[i];
            edgeLabeling[i] = constructor.get_edge_label(uv.first, uv.second);
        }
    }

    void exportMulticutOptions(py::module pyModule) {
        py::class_<MulticutOptions>(pyModule, "MulticutOptions")
            .def(
                py::init<const size_t, const std::string &, const std::string &, const std::string &,
                         const bool, const size_t, const size_t, const double, const double>(),
                py::arg("primalComputationInterval"), py::arg("standardReparametrization"), py::arg("roundingReparametrization"),
                py::arg("tightenReparametrization"), py::arg("tighten"), py::arg("tightenInterval"),
                py::arg("tightenIteration"), py::arg("tightenSlope"), py::arg("tightenConstraintsPercentage") 
            )
        ;
    }

    template<typename LABEL_TYPE, typename WEIGHT_TYPE>
    void exportMulticutT(py::module & pyModule) {

        typedef LABEL_TYPE LabelType;
        typedef WEIGHT_TYPE WeightType;
        
        // TODO nested vector is ugly, use proper multiarray instead
        pyModule.def("multicut", 
            [](
                const std::vector<std::pair<LabelType,LabelType>> & uvIds,
                const std::vector<WeightType> & weights,
                const MulticutOptions & multicutOptions
            ) {
                
                std::vector<bool> edgeLabels(uvIds.size());
                {
                    py::gil_scoped_release allowThreads;
                    SolverType solver(multicutOptions.toOptionsVector());
                    auto & multicutConstructor = solver.template GetProblemConstructor<0>();
                    getEdgeLabeling(uvIds, multicutConstructor, edgeLabels);
                }
                return edgeLabels;
            },
            py::arg("uvIds"), py::arg("weights"), py::arg("multicutOptions")
        );
    }

    void exportMulticut(py::module & pyModule) {
        exportMulticutOptions(pyModule);
        exportMulticutT<uint32_t, float>(pyModule);
    }
}
