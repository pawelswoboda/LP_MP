project(LP_MP)

cmake_minimum_required(VERSION 2.8.12)

set(LP_MP_VERSION_MAJOR 0)
set(LP_MP_VERSION_MINOR 1)

# C++11
add_compile_options(-std=c++14)

# compiler options
add_definitions(-march=native)

option(PARALLEL_OPTIMIZATION "Enable parallel optimization" OFF)
option(WITH_SAT_BASED_ROUNDING "SAT based rounding" ON)
option(LP_MP_TEST "build tests" OFF)

include_directories("${CMAKE_CURRENT_SOURCE_DIR}/external/libsimdpp")
include_directories("${CMAKE_CURRENT_SOURCE_DIR}/external/Catch/include")
include_directories("${CMAKE_CURRENT_SOURCE_DIR}/external/meta/include")
include_directories("${CMAKE_CURRENT_SOURCE_DIR}/external/cpp-sort/include")
include_directories("${CMAKE_CURRENT_SOURCE_DIR}/external/TCLAP/include")
include_directories("${CMAKE_CURRENT_SOURCE_DIR}/external/BCFW-Bundle")
add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/external/BCFW-Bundle")

include_directories("${CMAKE_CURRENT_SOURCE_DIR}/lib/ConicBundle/include")
link_directories("${CMAKE_CURRENT_SOURCE_DIR}/lib/ConicBundle/lib")

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

include_directories(include)

if(LP_MP_TEST)
  add_subdirectory(test)
endif()

