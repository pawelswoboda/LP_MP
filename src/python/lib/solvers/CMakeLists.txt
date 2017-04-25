addPythonModule(
    NESTED_NAME
        lp_mp/solvers 
    SOURCES
        solvers.cxx
        export_multicut.cxx
    LIBRRARIES
        m stdc++ pthread ${HDF5_LIBRARIES} ${HDF5_CXX_LIBRARIES})
