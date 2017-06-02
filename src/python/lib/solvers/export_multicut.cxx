#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>

#include <vector>

#include "visitors/standard_visitor.hxx"
#include "solvers/multicut/multicut.h"

namespace py = pybind11;

namespace LP_MP {
    using Rounder = KlRounder;
    
    //using FMC = FMC_MULTICUT<MessageSendingType::SRMP,Rounder>;
    using FMC = FMC_ODD_WHEEL_MULTICUT<MessageSendingType::SRMP,Rounder>;
    
    using SolverBase = Solver<FMC,LP,StandardTighteningVisitor,Rounder>;
    using SolverType = ProblemConstructorRoundingSolver<SolverBase>;
    
    //*****************************
    // multicut options
    //*****************************

    struct MulticutOptions {
        
        size_t primalComputationInterval;
        std::string standardReparametrization;
        std::string roundingReparametrization;
        std::string tightenReparametrization;
        bool tighten;
        size_t tightenInterval;
        size_t tightenIteration;
        double tightenSlope;
        double tightenConstraintsPercentage;
        size_t maxIter;
        double minDualImprovement;
        size_t minDualImprovementInterval;
        size_t timeout;
        size_t nThreads;
   
        std::vector<std::string> toOptionsVector() const {
            
            std::vector<std::string> options = {
              "export_multicut",
              "-i", " ", // empty input file
              "--primalComputationInterval", std::to_string(primalComputationInterval),
              "--standardReparametrization", standardReparametrization,
              "--roundingReparametrization", roundingReparametrization,
              "--tightenReparametrization",  tightenReparametrization,
              "--tightenInterval",           std::to_string(tightenInterval),
              "--tightenIteration",          std::to_string(tightenIteration),
              "--tightenSlope",              std::to_string(tightenSlope),
              "--tightenConstraintsPercentage", std::to_string(tightenConstraintsPercentage),
              "--maxIter", std::to_string(maxIter),
              #ifdef WITH_OPENMP
              ,"--numLpThreads", std::to_string(nThreads)
              #endif
            };
            if(tighten)
                options.push_back("--tighten");
            if(minDualImprovement > 0) {
                options.push_back("--minDualImprovement");
                options.push_back(std::to_string(minDualImprovement));
            }
            if(minDualImprovementInterval > 0) {
                options.push_back("--minDualImprovementInterval");
                options.push_back(std::to_string(minDualImprovementInterval));
            }
            if(timeout > 0) {
                options.push_back("--timeout");
                options.push_back(std::to_string(timeout));
            }
            return options;
        }
    };
    
    void exportMulticutOptions(py::module pyModule) {
        py::class_<MulticutOptions>(pyModule, "MulticutOptions")
            .def(py::init<>())
            .def_readwrite("primalComputationInterval",&MulticutOptions::primalComputationInterval)
            .def_readwrite("standardReparametrization",&MulticutOptions::standardReparametrization)
            .def_readwrite("roundingReparametrization",&MulticutOptions::roundingReparametrization)
            .def_readwrite("tightenReparametrization",&MulticutOptions::tightenReparametrization)
            .def_readwrite("tighten",&MulticutOptions::tighten)
            .def_readwrite("tightenInterval",&MulticutOptions::tightenInterval)
            .def_readwrite("tightenIteration",&MulticutOptions::tightenIteration)
            .def_readwrite("tightenSlope",&MulticutOptions::tightenSlope)
            .def_readwrite("tightenConstraintsPercentage",&MulticutOptions::tightenConstraintsPercentage)
            .def_readwrite("maxIter",&MulticutOptions::maxIter)
            .def_readwrite("minDualImprovement",&MulticutOptions::minDualImprovement)
            .def_readwrite("minDualImprovementInterval",&MulticutOptions::minDualImprovementInterval)
            .def_readwrite("timeout",&MulticutOptions::timeout)
            .def_readwrite("nThreads",&MulticutOptions::nThreads)
        ;
    }

    //*****************************
    // multicut solver
    //*****************************
    
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
       
       

    template<typename LABEL_TYPE, typename WEIGHT_TYPE>
    void exportMulticutT(py::module & pyModule) {

        typedef LABEL_TYPE LabelType;
        typedef WEIGHT_TYPE WeightType;
        
        // TODO nested vector is ugly, use proper multiarray instead
        pyModule.def("multicutImpl", 
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
                    constructMcProblem(uvIds, weights, multicutConstructor);
                    solver.Solve();
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
