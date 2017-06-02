cmake_minimum_required(VERSION 2.8.12)
project(LP_MP)

set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake/modules)

set(LP_MP_VERSION_MAJOR 0)
set(LP_MP_VERSION_MINOR 1)

# Options
OPTION(BUILD_GRAPH_MATCHING "Build graph matching" ON)
OPTION(BUILD_GRAPH_MATCHING_EVALUATION "Build graph matching evaluation" OFF)
OPTION(BUILD_MULTICUT "Build multicut" ON)
OPTION(BUILD_MULTICUT_EVALUATION "Build multicut evaluation" OFF)
OPTION(BUILD_GRAPHICAL_MODEL "Build graphical model" ON)
OPTION(BUILD_DISCRETE_TOMOGRAPHY "Build discrete tomography solver" ON)
OPTION(BUILD_CELL_TRACKING "Build cell tracking solver" OFF)
OPTION(BUILD_MORAL_LINEAGE_TRACING "Build moral lineage tracing solver" OFF)
OPTION(BUILD_DISCRETE_TOMOGRAPHY_EVALUATION "Build discrete tomography evaluation" OFF)
OPTION(BUILD_TESTS "Build tests" ON)
OPTION(WITH_GUROBI "LP interface to gurobi" OFF)
OPTION(WITH_CPLEX "LP interface to Cplex" OFF)
OPTION(WITH_SAT_BASED_ROUNDING "Use the glucose SAT solver to decode a primal solution based on reparametrization" OFF)
OPTION(PARALLEL_OPTIMIZATION "Enable parallel optimization" OFF)

# C++14
add_compile_options(-std=c++14)

# compiler options
add_definitions(-DIL_STD) # legacy setting for CPLEX
if(CMAKE_BUILD_TYPE STREQUAL "Release")
   #add_definitions(-ffast-math -fno-finite-math-only) # adding only -ffast-math will result in infinity and nan not being checked (but e.g. graph matching and discrete tomography have infinite costs)
   add_definitions(-march=native)
endif()

option(BUILD_PYTHON OFF)

# automatically downloaded repositories
# can this possibly be done in one place only, i.e. in the superbuild?
include_directories("${CMAKE_CURRENT_SOURCE_DIR}/external/meta/include")
include_directories("${CMAKE_CURRENT_SOURCE_DIR}/external/Catch/include")
include_directories("${CMAKE_CURRENT_SOURCE_DIR}/external/cpp_sort/include")
include_directories("${CMAKE_CURRENT_SOURCE_DIR}/external/opengm/include")
include_directories("${CMAKE_CURRENT_SOURCE_DIR}/external/PEGTL")
include_directories("${CMAKE_CURRENT_SOURCE_DIR}/external/graph/include")
include_directories("${CMAKE_CURRENT_SOURCE_DIR}/external/cereal/include")
include_directories("${CMAKE_CURRENT_SOURCE_DIR}/external/tclap/include")
#include_directories("${CMAKE_CURRENT_SOURCE_DIR}/external/lingeling-ala")

# TODO this seems to be a pretty old version, if the license permits it make a new upload ?! however, this seems to be compatible with the one that was downloaded before
include(ExternalProject)
ExternalProject_ADD(
   lingeling
   URL "http://fmv.jku.at/lingeling/lingeling-bal-2293bef-151109.tar.gz"
   BUILD_IN_SOURCE 1
   CONFIGURE_COMMAND "./configure.sh"
   BUILD_COMMAND "make"
   INSTALL_COMMAND ""
   )
include_directories("${CMAKE_CURRENT_BINARY_DIR}/lingeling-prefix/src/lingeling") 
link_directories("${CMAKE_CURRENT_BINARY_DIR}/lingeling-prefix/src/lingeling")

#include_directories("${CMAKE_CURRENT_BINARY_DIR}/Dependencies/Build/CryptoMiniSat_Project/include")
#link_directories("${CMAKE_CURRENT_BINARY_DIR}/Dependencies/Build/CryptoMiniSat_Project/lib")

#add_subdirectory("${CMAKE_CURRENT_BINARY_DIR}/Dependencies/Source/LEMON_Project")
#set(LEMON_INCLUDE_DIRS
#   ${LEMON_SOURCE_ROOT_DIR}
#   ${CMAKE_BINARY_DIR}/deps/lemon
#   )
#include_directories(${LEMON_INCLUDE_DIRS})
#include_directories("${CMAKE_CURRENT_BINARY_DIR}/Dependencies/Source/LEMON_Project/lemon")

#include_directories("${CMAKE_CURRENT_BINARY_DIR}/Dependencies/Source/Hana_Project/include")
#include_directories("${CMAKE_CURRENT_BINARY_DIR}/Dependencies/Source/CS2_CPP_Project")

# manually downloaded repositories of Kolmogorov's code. How to automate?
#add_subdirectory(lib/MinCost)

#-----------
# pybindings
#-----------
include_directories( "${CMAKE_CURRENT_SOURCE_DIR}/external/pybind11/include")
add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/external/pybind11")


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

# TODO refactor this too !
include_directories(lib)
add_subdirectory(lib)

add_subdirectory(src)
