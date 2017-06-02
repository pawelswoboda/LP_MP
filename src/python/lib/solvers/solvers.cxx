#include <pybind11/pybind11.h>

namespace py = pybind11;

namespace LP_MP {
    void exportMulticut(py::module &);
}

PYBIND11_PLUGIN(_solvers) {
    py::module pyModule("_solvers", "LP_MP solvers module");

    LP_MP::exportMulticut(pyModule);
        
    return pyModule.ptr();
}

