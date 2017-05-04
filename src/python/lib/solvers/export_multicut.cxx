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
    
    //*****************************
    // multicut options
    //*****************************

    struct MulticutOptions {
        MulticutOptions(
            const size_t primalComputationInterval,
            const std::string & standardReparametrization,
            const std::string & roundingReparametrization,
            const std::string & tightenReparametrization,
            const bool tighten,
            const size_t tightenInterval,
            const size_t tightenIteration,
            const double tightenSlope,
            const double tightenConstraintsPercentage,
            const size_t maxIter,                    // default : 1000
            const double minDualImprovement,         // default 0 -> not used
            const size_t minDualImprovementInterval, // default 0 -> not used
            const size_t timeout,                    // default 0 -> not used
            const size_t nThreads
        ) :
            primalComputationInterval_(primalComputationInterval),
            standardReparametrization_(standardReparametrization),
            roundingReparametrization_(roundingReparametrization),
            tightenReparametrization_(tightenReparametrization),
            tighten_(tighten),
            tightenInterval_(tightenInterval),
            tightenIteration_(tightenIteration),
            tightenSlope_(tightenSlope),
            tightenConstraintsPercentage_(tightenConstraintsPercentage),
            maxIter_(maxIter),
            minDualImprovement_(minDualImprovement),
            minDualImprovementInterval_(minDualImprovementInterval),
            timeout_(timeout),
            nThreads_(nThreads)
        {}

        std::vector<std::string> toOptionsVector() const {
            std::vector<std::string> options = {
              "export_multicut",
              "-i", " ", // empty input file
              "--primalComputationInterval", std::to_string(primalComputationInterval_),
              "--standardReparametrization", standardReparametrization_,
              "--roundingReparametrization", roundingReparametrization_,
              "--tightenReparametrization",  tightenReparametrization_,
              "--tightenInterval",           std::to_string(tightenInterval_),
              "--tightenIteration",          std::to_string(tightenIteration_),
              "--tightenSlope",              std::to_string(tightenSlope_),
              "--tightenConstraintsPercentage", std::to_string(tightenConstraintsPercentage_),
              "--maxIter", std::to_string(maxIter_),
              "--numLpThreads", std::to_string(nThreads_)
            };
            if(tighten_)
                options.push_back("--tighten");
            if(minDualImprovement_ > 0) {
                options.push_back("--minDualImprovement");
                options.push_back(std::to_string(minDualImprovement_));
            }
            if(minDualImprovementInterval_ > 0) {
                options.push_back("--minDualImprovementInterval");
                options.push_back(std::to_string(minDualImprovementInterval_));
            }
            if(timeout_ > 0) {
                options.push_back("--timeout");
                options.push_back(std::to_string(timeout_));
            }
            return options;
        }

        size_t primalComputationInterval() const {return primalComputationInterval_;}
        const std::string & standardReparametrization() const {return standardReparametrization_;}
        const std::string & roundingReparametrization() const {return roundingReparametrization_;}
        const std::string & tightenReparametrization() const {return tightenReparametrization_;}
        bool tighten() const {return tighten_;}
        size_t tightenInterval() const {return tightenInterval_;}
        size_t tightenIteration() const {return tightenIteration_;}
        double tightenSlope() const {return tightenSlope_;}
        double tightenConstraintsPercentage() const {return tightenConstraintsPercentage_;}
        size_t maxIter() const {return maxIter_;}
        double minDualImprovement() const {return minDualImprovement_;}
        size_t minDualImprovementInterval() const {return minDualImprovementInterval_;}
        size_t timeout() const {return timeout_;}
        size_t nThreads() const {return nThreads_;}
        
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
        const size_t maxIter_;
        const double minDualImprovement_;
        const size_t minDualImprovementInterval_;
        const size_t timeout_;
        const size_t nThreads_;
    };
    
    //
    // TODO these are all parameters, some of them should be exposed as well
    // (max runtime, max memory )
    //
    //[--tightenSlope <positive real number smaller 1>] // check
    //[--tightenMinDualImprovementInterval <positive integer>]
    //[--tightenMinDualImprovement <positive real>]
    //[--tightenConstraintsPercentage <positive real>]
    //[--tightenConstraintsMax <positive integer>]
    //[--tightenInterval <positive integer>] 
    //[--tightenIteration <positive integer>]
    //[--tightenReparametrization <(uniform|anisotropic)>]
    //[--tighten]
    //[--roundingReparametrization <{anisotropic|uniform}>]
    //[--standardReparametrization <{anisotropic|uniform}>]
    //[--minDualImprovementInterval <strictly positive integer>]
    //[--minDualImprovement <positive real number>]
    //[--lowerBoundComputationInterval <strictly positive integer>]
    //[--primalComputationInterval <strictly positive integer>]
    //[--timeout <strictly positive integer>] 
    //[--maxMemory <positive integer>]
    //[--maxIter <strictly positive integer>]
    //[-o <file name>] 
    //[-i <file name>] [--] [--version] [-h]

    void exportMulticutOptions(py::module pyModule) {
        py::class_<MulticutOptions>(pyModule, "MulticutOptions")
            .def(
                py::init<
                    size_t,
                    std::string,
                    std::string,
                    std::string,
                    bool,
                    size_t,
                    size_t,
                    double,
                    double,
                    size_t,
                    double,
                    size_t,
                    size_t,
                    size_t
                >(),
                py::arg_t<size_t>("primalComputationInterval", 100), // check
                py::arg_t<std::string>("standardReparametrization", "anisotropic"), // check
                py::arg_t<std::string>("roundingReparametrization", "damped_uniform"), // damped_uniform ?!
                py::arg_t<std::string>("tightenReparametrization",  "damped_uniform"), // dampeld_uniform ?!
                py::arg_t<bool>("tighten", true),
                py::arg_t<size_t>("tightenInterval", 100),
                py::arg_t<size_t>("tightenIteration", 2),
                py::arg_t<double>("tightenSlope", 0.05),
                py::arg_t<double>("tightenConstraintsPercentage", 0.1),
                py::arg_t<size_t>("maxIter", 1000), 
                py::arg_t<double>("minDualImprovement", 0.),
                py::arg_t<size_t>("minDualImprovementInterval", 0),
                py::arg_t<size_t>("timeout", 0),
                py::arg_t<size_t>("nThreads", 1)
            )
            .def_property_readonly("primalComputationInterval",&MulticutOptions::primalComputationInterval)
            .def_property_readonly("standardReparametrization",&MulticutOptions::standardReparametrization)
            .def_property_readonly("roundingReparametrization",&MulticutOptions::roundingReparametrization)
            .def_property_readonly("tightenReparametrization",&MulticutOptions::tightenReparametrization)
            .def_property_readonly("tighten",&MulticutOptions::tighten)
            .def_property_readonly("tightenInterval",&MulticutOptions::tightenInterval)
            .def_property_readonly("tightenIteration",&MulticutOptions::tightenIteration)
            .def_property_readonly("tightenSlope",&MulticutOptions::tightenSlope)
            .def_property_readonly("tightenConstraintsPercentage",&MulticutOptions::tightenConstraintsPercentage)
            .def_property_readonly("maxIter",&MulticutOptions::maxIter)
            .def_property_readonly("minDualImprovement",&MulticutOptions::minDualImprovement)
            .def_property_readonly("timeout",&MulticutOptions::timeout)
            .def_property_readonly("nThreads",&MulticutOptions::nThreads)
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
