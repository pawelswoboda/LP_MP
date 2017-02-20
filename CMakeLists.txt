cmake_minimum_required(VERSION 2.8)
set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake/modules)

option(DOWNLOAD_DEPENDENCIES "Download and compile dependencies automatically" ON)
# Options
OPTION(BUILD_GRAPH_MATCHING "Build graph matching" ON)
OPTION(BUILD_GRAPH_MATCHING_EVALUATION "Build graph matching evaluation" OFF)
OPTION(BUILD_MULTICUT "Build multicut" ON)
OPTION(BUILD_MULTICUT_EVALUATION "Build multicut evaluation" OFF)
OPTION(BUILD_GRAPHICAL_MODEL "Build graphical model" ON)
OPTION(BUILD_DISCRETE_TOMOGRAPHY "Build discrete tomography solver" ON)
OPTION(BUILD_CELL_TRACKING "Build cell tracking solver" OFF)
OPTION(BUILD_DISCRETE_TOMOGRAPHY_EVALUATION "Build discrete tomography evaluation" OFF)
OPTION(BUILD_TESTS "Build tests" ON)
OPTION(WITH_GUROBI "LP interface to gurobi" OFF)
OPTION(WITH_CPLEX "LP interface to Cplex" OFF)
<<<<<<< HEAD
OPTION(WITH_SAT_BASED_ROUNDING "Use the glucose SAT solver to decode a primal solution based on reparametrization" OFF)
=======
>>>>>>> 50f351e47fc311f373f1af914ee3938caee18993
OPTION(PARALLEL_OPTIMIZATION "Enable parallel optimization" OFF)

if(DOWNLOAD_DEPENDENCIES)
   # download external projects here 
   project(DOWNLOAD_DEPENDENCIES_PROJECT NONE)
   include(cmake/SuperBuild.cmake)
   return()
else()
   project(LP_MP)
   include(cmake/LP_MP.cmake)
endif()
