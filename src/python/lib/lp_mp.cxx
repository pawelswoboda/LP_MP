#include <pybind11/pybind11.h>

namespace py = pybind11;

PYBIND11_PLUGIN(_lp_mp) {
    py::module pyModule("_lp_mp", "LP_MP python bindings");
    return pyModule.ptr();
}

