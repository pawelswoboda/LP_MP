add_subdirectory(solvers) # TODO flag for building solvers ?
if(BUILD_TESTS)
    add_subdirectory(test)
endif()
if(BUILD_PYTHON)
    add_subdirectory(python)
endif()
#add_subdirectory(external_lib)
