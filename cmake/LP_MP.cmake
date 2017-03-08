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

#set(CMAKE_FIND_LIBRARY_SUFFIXES ".a")
#set(BUILD_SHARED_LIBRARIES OFF)
#set(CMAKE_EXE_LINKER_FLAG "-static")

# Vc for SIMD
#find_package(Vc 1.2.0 REQUIRED PATHS "${CMAKE_CURRENT_BINARY_DIR}/Dependencies/Build/Vc_Project/cmake" "${CMAKE_CURRENT_BINARY_DIR}/Dependencies/Source/Vc_Project/cmake")
#include_directories(${Vc_INCLUDE_DIR}) 
#add_definitions(${Vc_DEFINITIONS})
#link_directories(${Vc_LIB_DIR})

# automatically downloaded repositories
# can this possibly be done in one place only, i.e. in the superbuild?
include_directories("${CMAKE_CURRENT_BINARY_DIR}/Dependencies/Source/meta_Project/include")
include_directories("${CMAKE_CURRENT_BINARY_DIR}/Dependencies/Source/Catch_Project/include")
include_directories("${CMAKE_CURRENT_BINARY_DIR}/Dependencies/Source/cpp_sort_Project/include")
include_directories("${CMAKE_CURRENT_BINARY_DIR}/Dependencies/Source/OpenGM_Project/include")
include_directories("${CMAKE_CURRENT_BINARY_DIR}/Dependencies/Source/PEGTL_Project")
include_directories("${CMAKE_CURRENT_BINARY_DIR}/Dependencies/Source/Andres_Project/include")
include_directories("${CMAKE_CURRENT_BINARY_DIR}/Dependencies/Source/Cereal_Project/include")
include_directories("${CMAKE_CURRENT_BINARY_DIR}/Dependencies/Source/TCLAP_Project/include")
#include_directories("${CMAKE_CURRENT_BINARY_DIR}/Dependencies/Build/CryptoMiniSat_Project/include")
#link_directories("${CMAKE_CURRENT_BINARY_DIR}/Dependencies/Build/CryptoMiniSat_Project/lib")
include_directories("${CMAKE_CURRENT_BINARY_DIR}/Dependencies/Source/Lingeling_Project")
link_directories("${CMAKE_CURRENT_BINARY_DIR}/Dependencies/Source/Lingeling_Project")

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


# HDF5 for reading OpenGM and Andres models
# set (HDF5_USE_STATIC_LIBRARIES ON)
find_package(HDF5 1.8.15 REQUIRED)
INCLUDE_DIRECTORIES (${HDF5_INCLUDE_DIR})
message(STATUS ${HDF5_LIBRARIES})
message(STATUS ${HDF5_INCLUDE_DIR})

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

