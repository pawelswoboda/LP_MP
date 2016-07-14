cmake_minimum_required(VERSION 2.6)
set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake/modules)
set(LP_MP_VERSION_MAJOR 0)
set(LP_MP_VERSION_MINOR 1)

# global headers
file(GLOB_RECURSE headers include/*.hxx)
file(GLOB_RECURSE sources src/*.cpp)
include_directories(include)
include_directories(lib)
include_directories(.)

# C++11
SET(C++_STD_FLAG "c++14")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=${C++_STD_FLAG}")
MARK_AS_ADVANCED(C++_STD_FLAG)

# compiler options
add_definitions(-DIL_STD)
add_definitions(-ffast-math)
add_definitions(-march=native)

# packages

# GUROBI
include_directories("/home/pswoboda/gurobi/linux64/include")
# more robust: use add_library with IMPORTED
LINK_DIRECTORIES("/home/pswoboda/gurobi/linux64/lib")
#LINK_DIRECTORIES("/home/pswoboda/gurobi/linux64/src/build")
SET(GUROBI_LIB "/home/pswoboda/gurobi/linux64/lib/libgurobi_c++.a /home/pswoboda/gurobi/linux64/lib/libgurobi65.so")

# Vc SIMD library
set(Vc_DIR "/home/pswoboda/repo/LP_MP/lib/Vc-1.2.0/Release/cmake")
set(Vc_LIBRARIES "/home/pswoboda/repo/LP_MP/lib/Vc-Installation/lib/libVc.a") # needed?
find_package(Vc 1.0.0 REQUIRED)
#if(Vc_FOUND)
#   message(STATUS "Vc found")
   include_directories(${Vc_INCLUDE_DIR}) 
   add_definitions(${Vc_DEFINITIONS})
   link_directories(${Vc_LIB_DIR})
   #add_library(Vc STATIC Vc/Vc-0.7.5-Installation/lib/libVc.a)
#endif()

# HDF5 for reading OpenGM models
find_package(HDF5)
#find_package(HDF5 NAMES hdf5 )
#find_package (HDF5 NAMES hdf5 COMPONENTS C static)
INCLUDE_DIRECTORIES (${HDF5_INCLUDE_DIR})
set (LINK_LIBS ${LINK_LIBS} ${HDF5_C_STATIC_LIBRARY})
message(STATUS ${HDF5_INCLUDE_DIR})

# Options
OPTION(BUILD_GRAPH_MATCHING "Build graph matching" OFF)
OPTION(BUILD_MULTICUT "Build multicut" OFF)
OPTION(BUILD_GRAPHICAL_MODEL "Build graphical model" OFF)
OPTION(BUILD_DISCRETE_TOMOGRAPHY "Build discrete tomography solver" ON)
OPTION(BUILD_EVALUATION "Build evaluation" OFF)


IF(UNIX AND NOT APPLE)
   find_library(TR rt)
   set(LINK_RT true)
   message(STATUS "Linking to RT is enabled")
else()
   set(LINK_RT false)
   message(STATUS "Linking to RT is disabled")
endif()

#add_subdirectory(lib/cpp-sort)
include_directories(lib/cpp-sort/include)
#add_subdirectory(lib/lemon-1.3.1)
include_directories(lib/lemon-1.3.1)
add_subdirectory(lib/MinCost)
include_directories(lib/PEGTL-master)
include_directories(lib/opengm/include)
include_directories(lib/spdlog/include)
include_directories(lib/andres/include)
#include_directories(lib/hana-0.7.0/include)
add_subdirectory(solvers)
