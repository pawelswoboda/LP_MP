#include <pybind11/pybind11.h>

namespace py = pybind11;

namespace LP_MP {
    void exportMulticut(py::module &);
}

PYBIND11_PLUGIN(_lp_mp) {
    py::module pyModule("_lp_mp", "LP_MP python bindings");

    LP_MP::exportMulticut(pyModule);
        
    return pyModule.ptr();
}

