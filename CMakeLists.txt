project(LP_MP)

cmake_minimum_required(VERSION 2.8.12)

set(LP_MP_VERSION_MAJOR 0)
set(LP_MP_VERSION_MINOR 1)

# C++11
add_compile_options(-std=c++14)

# compiler options
add_definitions(-march=native)

set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

# Options
option(BUILD_GRAPH_MATCHING "Build graph matching" ON)
option(BUILD_GRAPH_MATCHING_EVALUATION "Build graph matching evaluation" OFF)
option(BUILD_MULTICUT "Build multicut" ON)
option(BUILD_MULTICUT_EVALUATION "Build multicut evaluation" OFF)
option(BUILD_MAX_CUT "Build max-cut" ON)
option(BUILD_GRAPHICAL_MODEL "Build graphical model" ON)
option(BUILD_DISCRETE_TOMOGRAPHY "Build discrete tomography solver" ON)
option(BUILD_CELL_TRACKING "Build cell tracking solver" OFF)
option(BUILD_MORAL_LINEAGE_TRACING "Build moral lineage tracing solver" OFF)
option(BUILD_CANDIDATE_MULTICUT "Build candidate multicut solver" OFF)
option(BUILD_DISCRETE_TOMOGRAPHY_EVALUATION "Build discrete tomography evaluation" OFF)
option(BUILD_TESTS "Build tests" ON)
option(WITH_GUROBI "LP interface to gurobi" OFF)
option(WITH_CPLEX "LP interface to Cplex" OFF)
option(WITH_SAT_BASED_ROUNDING "Use the lingeling SAT solver to decode a primal solution based on reparametrization" OFF)
option(PARALLEL_OPTIMIZATION "Enable parallel optimization" OFF)

include_directories("${CMAKE_CURRENT_SOURCE_DIR}/external/MQLib/include")
find_library(MAX_CUT_LIB NAMES MQLib.a PATHS "${CMAKE_CURRENT_SOURCE_DIR}/external/MQLib/bin")
#link_directories("${CMAKE_CURRENT_SOURCE_DIR}/external/MQLib/bin/")

include_directories("${CMAKE_CURRENT_SOURCE_DIR}/external/libsimdpp")
include_directories("${CMAKE_CURRENT_SOURCE_DIR}/external/arboricity")
add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/external/arboricity")
include_directories("${CMAKE_CURRENT_SOURCE_DIR}/external/andres_graph/include")
include_directories("${CMAKE_CURRENT_SOURCE_DIR}/external/Catch/include")
include_directories("${CMAKE_CURRENT_SOURCE_DIR}/external/meta/include")
include_directories("${CMAKE_CURRENT_SOURCE_DIR}/external/PEGTL")
include_directories("${CMAKE_CURRENT_SOURCE_DIR}/external/opengm/include")
include_directories("${CMAKE_CURRENT_SOURCE_DIR}/external/cpp-sort/include")
include_directories("${CMAKE_CURRENT_SOURCE_DIR}/external/TCLAP/include")
#find_library(ARBORICITY NAMES libarboricity.a PATHS "${CMAKE_CURRENT_SOURCE_DIR}/external/arboricity")

include(ExternalProject)
list(APPEND DEPENDENCIES lingeling_project)
ExternalProject_ADD(lingeling_project
  SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/external/lingeling"
  CONFIGURE_COMMAND "${CMAKE_CURRENT_SOURCE_DIR}/external/lingeling/configure"
  BUILD_COMMAND make
  BUILD_IN_SOURCE 1
  INSTALL_COMMAND ""
  )
ExternalProject_Get_Property(lingeling_project install_dir)
include_directories("${CMAKE_CURRENT_SOURCE_DIR}/external/lingeling")
link_directories("${CMAKE_CURRENT_SOURCE_DIR}/external/lingeling")

if(WITH_GUROBI)
  add_definitions(-DIL_STD) # legacy setting for CPLEX
endif()


# HDF5 for reading OpenGM and Andres models
# set (HDF5_USE_STATIC_LIBRARIES ON)
if(BUILD_MULTICUT OR BUILD_MULTICUT_EVALUATION OR BUILD_GRAPHICAL_MODEL)
   find_package(HDF5 1.8.15 REQUIRED)
   include_directories (${HDF5_INCLUDE_DIR})
   add_definitions(${HDF5_DEFINITIONS})
   message(STATUS ${HDF5_LIBRARIES})
   message(STATUS ${HDF5_INCLUDE_DIR})
endif()

# GUROBI
if(WITH_GUROBI)
   find_package(Gurobi REQUIRED)
   add_definitions(-DWITH_GUROBI)
endif(WITH_GUROBI)

# CPLEX
if(WITH_CPLEX)
   find_package(Cplex REQUIRED)
   add_definitions(-DWITH_CPLEX)
endif(WITH_CPLEX)

if(WITH_SAT_BASED_ROUNDING)
   add_definitions(-DWITH_SAT)
endif()

# Parallelisation support
if(PARALLEL_OPTIMIZATION)

  add_definitions(-DLP_MP_PARALLEL) 

  FIND_PACKAGE(OpenMP REQUIRED)
  if(OPENMP_FOUND)
     message("OPENMP FOUND")
     set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
     set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
     set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
  endif()

endif(PARALLEL_OPTIMIZATION)

IF(UNIX AND NOT APPLE)
   find_library(TR rt)
   set(LINK_RT true)
   message(STATUS "Linking to RT is enabled")
else()
   set(LINK_RT false)
   message(STATUS "Linking to RT is disabled")
endif()

file(GLOB_RECURSE headers include/*.hxx)
include_directories(include)
include_directories(lib)
include_directories(.)
add_subdirectory(solvers)
add_subdirectory(lib)
add_subdirectory(test)

