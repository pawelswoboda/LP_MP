cmake_minimum_required(VERSION 2.8)
set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake/modules)

# rename to DOWNLOAD_DEPENDENCIES
option(DOWNLOAD_DEPENDENCIES "Download and compile dependencies automatically" ON)
# Options
OPTION(BUILD_GRAPH_MATCHING "Build graph matching" ON)
OPTION(BUILD_MULTICUT "Build multicut" ON)
OPTION(BUILD_GRAPHICAL_MODEL "Build graphical model" ON)
OPTION(BUILD_DISCRETE_TOMOGRAPHY "Build discrete tomography solver" ON)
OPTION(BUILD_EVALUATION "Build evaluation" ON)
OPTION(BUILD_TESTS "Build tests" ON)
OPTION(BUILD_GUROBI_INTERFACE "LP interface to gurobi" OFF)

if(DOWNLOAD_DEPENDENCIES)
   # download external projects here
   project(DOWNLOAD_DEPENDENCIES_PROJECT NONE)
   include(cmake/SuperBuild.cmake)
   return()
else()
   project(LP_MP)
   include(cmake/LP_MP.cmake)
endif()
