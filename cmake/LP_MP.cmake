project(LP_MP)

cmake_minimum_required(VERSION 2.8.12)

set(LP_MP_VERSION_MAJOR 0)
set(LP_MP_VERSION_MINOR 1)

# C++11
add_compile_options(-std=c++14)

# compiler options
add_definitions(-DIL_STD)
add_definitions(-ffast-math)
add_definitions(-march=native)

# Vc for SIMD
find_package(Vc 1.2.0 REQUIRED PATHS "${CMAKE_CURRENT_BINARY_DIR}/Dependencies/Build/Vc_Project/cmake" "${CMAKE_CURRENT_BINARY_DIR}/Dependencies/Source/Vc_Project/cmake")
include_directories(${Vc_INCLUDE_DIR}) 
add_definitions(${Vc_DEFINITIONS})
link_directories(${Vc_LIB_DIR})

# automatically downloaded repositories
# can this possibly be done in one place only, i.e. in the superbuild?
include_directories("${CMAKE_CURRENT_BINARY_DIR}/Dependencies/Source/meta_Project/include")
include_directories("${CMAKE_CURRENT_BINARY_DIR}/Dependencies/Source/Catch_Project/include")
include_directories("${CMAKE_CURRENT_BINARY_DIR}/Dependencies/Source/cpp_sort_Project/include")
include_directories("${CMAKE_CURRENT_BINARY_DIR}/Dependencies/Source/OpenGM_Project/include")
include_directories("${CMAKE_CURRENT_BINARY_DIR}/Dependencies/Source/PEGTL_Project")
include_directories("${CMAKE_CURRENT_BINARY_DIR}/Dependencies/Source/Andres_Project/include")
include_directories("${CMAKE_CURRENT_BINARY_DIR}/Dependencies/Source/TCLAP_Project/include")

#add_subdirectory("${CMAKE_CURRENT_BINARY_DIR}/Dependencies/Source/LEMON_Project")
#set(LEMON_INCLUDE_DIRS
#   ${LEMON_SOURCE_ROOT_DIR}
#   ${CMAKE_BINARY_DIR}/deps/lemon
#   )
#include_directories(${LEMON_INCLUDE_DIRS})
#include_directories("${CMAKE_CURRENT_BINARY_DIR}/Dependencies/Source/LEMON_Project/lemon")

include_directories("${CMAKE_CURRENT_BINARY_DIR}/Dependencies/Source/Hana_Project/include")
include_directories("${CMAKE_CURRENT_BINARY_DIR}/Dependencies/Source/CS2_CPP_Project")

# manually downloaded repositories of Kolmogorov's code. How to automate?
#add_subdirectory(lib/MinCost)

# HDF5 for reading OpenGM and Andres models 
find_package(HDF5 1.8.15 REQUIRED)
#find_package(HDF5 NAMES hdf5 )
#find_package (HDF5 NAMES hdf5 COMPONENTS C static)
INCLUDE_DIRECTORIES (${HDF5_INCLUDE_DIR})
set (LINK_LIBS ${LINK_LIBS} ${HDF5_C_STATIC_LIBRARY})
message(STATUS ${HDF5_INCLUDE_DIR})

# GUROBI
if(WITH_GUROBI)
   find_package(Gurobi REQUIRED)
endif(WITH_GUROBI)

# CPLEX
OPTION(WITH_CPLEX "Activate CPLEX-Code" OFF)
if(WITH_CPLEX)
  find_package(Cplex)
endif(WITH_CPLEX)

# LOCALSOLVER
OPTION(WITH_LOCALSOLVER "Activate LocalSolver-Code" OFF)
if(WITH_LOCALSOLVER)
  find_package(LocalSolver)
endif(WITH_LOCALSOLVER)

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
add_subdirectory(test)


message("build solvers")

