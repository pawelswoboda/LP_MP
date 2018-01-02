project(LP_MP)

cmake_minimum_required(VERSION 2.8.12)

set(LP_MP_VERSION_MAJOR 0)
set(LP_MP_VERSION_MINOR 1)

# C++11
add_compile_options(-std=c++17)

# compiler options
add_definitions(-march=native)

option(PARALLEL_OPTIMIZATION "Enable parallel optimization" OFF)

include(ExternalProject)
externalproject_add( conicBundle_Project
  URL http://www-user.tu-chemnitz.de/~helmberg/ConicBundle/CB_v0.3.11.tgz
  URL_HASH SHA256=f5a6774c24b2f71c31b5bdcadc0ac310f3862275ea4097a4ca3d6731876701c8
  PATCH_COMMAND patch -s -d ${PROJECT_BINARY_DIR}/external/ConicBundle -p1 < ${PROJECT_SOURCE_DIR}/external/ConicBundle-v0.3.11.patch -N -r -
  SOURCE_DIR "${PROJECT_BINARY_DIR}/external/ConicBundle"
  CONFIGURE_COMMAND ""
  BUILD_COMMAND "make"
  BUILD_IN_SOURCE 1
  INSTALL_COMMAND ""
  )

add_library(CONIC_BUNDLE STATIC IMPORTED GLOBAL)
set_target_properties(CONIC_BUNDLE PROPERTIES INCLUDE_DIRECTORIES "${PROJECT_BINARY_DIR}/external/ConicBundle/include")
set_target_properties(CONIC_BUNDLE PROPERTIES IMPORTED_LOCATION "${PROJECT_BINARY_DIR}/external/ConicBundle/lib/libcb.a")

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

#IF(UNIX AND NOT APPLE)
#   find_library(TR rt)
#   set(LINK_RT true)
#   message(STATUS "Linking to RT is enabled")
#else()
#   set(LINK_RT false)
#   message(STATUS "Linking to RT is disabled")
#endif()


add_library(LP_MP INTERFACE)
target_include_directories(LP_MP INTERFACE include)
target_include_directories(LP_MP INTERFACE external/TCLAP/include)

target_include_directories(LP_MP INTERFACE "external/libsimdpp")
target_include_directories(LP_MP INTERFACE "external/meta/include")
target_include_directories(LP_MP INTERFACE "external/cpp-sort/include")
target_include_directories(LP_MP INTERFACE "external/TCLAP/include")
target_include_directories(LP_MP INTERFACE "external/BCFW-Bundle/include")
target_include_directories(LP_MP INTERFACE "${PROJECT_BINARY_DIR}/external/ConicBundle/include")

add_subdirectory("external/BCFW-Bundle") 
add_subdirectory("external/DD_ILP")
target_link_libraries(LP_MP INTERFACE DD_ILP)

enable_testing()
add_subdirectory(test) 
