cmake_minimum_required(VERSION 2.8)
set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake/modules)

option(USE_SUPERBUILD "Download and compile dependencies automatically" ON)
# Options
OPTION(BUILD_GRAPH_MATCHING "Build graph matching" ON)
OPTION(BUILD_MULTICUT "Build multicut" ON)
OPTION(BUILD_GRAPHICAL_MODEL "Build graphical model" ON)
OPTION(BUILD_DISCRETE_TOMOGRAPHY "Build discrete tomography solver" ON)
OPTION(BUILD_EVALUATION "Build evaluation" ON)
OPTION(BUILD_TESTS "Build tests" ON)

OPTION(MULTICUT_DATASETS "Download multicut datasets" ON)
if(MULTICUT_DATASETS)

endif()

if(USE_SUPERBUILD)
   # download external projects here
   project (SUPERBUILD NONE)
   include(cmake/SuperBuild.cmake)
   return()
else()
   project(LP_MP)
   include(cmake/LP_MP.cmake)
endif()


#include_directories(lib/PEGTL-master)
#include_directories(lib/opengm/include)
#include_directories(lib/spdlog/include)
#include_directories(lib/hana-0.7.0/include)
