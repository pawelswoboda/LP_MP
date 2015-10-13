cmake_minimum_required(VERSION 2.6)
set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake/modules)
set(LP_MP_VERSION_MAJOR 0)
set(LP_MP_VERSION_MINOR 1)

# global headers
file(GLOB_RECURSE headers include/*.hxx)
file(GLOB_RECURSE sources src/*.cpp)
include_directories(include)
include_directories(lib)

# C++11
SET(C++_STD_FLAG "c++11" CACHE STRING "e.g. c++11,c++14")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=${C++_STD_FLAG}")
MARK_AS_ADVANCED(C++_STD_FLAG)

# compiler options
add_definitions(-DIL_STD)
add_definitions(-ffast-math)


# packages
set(Vc_DIR "/home/pswoboda/svn/praca doktorska/code/LP_MP/lib/Vc/Vc-0.7.5-Installation/lib/cmake/Vc")
set(Vc_LIBRARIES "/home/pswoboda/svn/praca doktorska/code/LP_MP/lib/Vc/Vc-0.7.5-Installation/lib/libVc.a") # needed?
find_package(Vc 0.7.5 REQUIRED)
if(Vc_FOUND)
   message(STATUS "Vc found")
   include_directories(${Vc_INCLUDE_DIR}) 
   add_definitions(${Vc_DEFINITIONS})
   link_directories(${Vc_LIB_DIR})
   #add_library(Vc STATIC Vc/Vc-0.7.5-Installation/lib/libVc.a)
endif()

# Options
OPTION(BUILD_GRAPH_MATCHING "Build graph matching" ON)
OPTION(BUILD_MULTICUT "Build multicut" OFF)
OPTION(BUILD_GRAPHICAL_MODEL "Build graphical model" OFF)


IF(UNIX AND NOT APPLE)
   find_library(TR rt)
   set(LINK_RT true)
   message(STATUS "Linking to RT is enabled")
else()
   set(LINK_RT false)
   message(STATUS "Linking to RT is disabled")
endif()

add_subdirectory(solvers)
add_subdirectory(lib)
