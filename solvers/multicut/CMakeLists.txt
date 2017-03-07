#FIND_PACKAGE(OpenMP REQUIRED)
#if(OPENMP_FOUND)
#   message("OPENMP FOUND")
#   set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
#   set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
#   set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
#endif()


SET(SOURCE_FILES 
   multicut_srmp.cpp 
   multicut_opengm_srmp_cycle.cpp 
   multicut_opengm_srmp_cycle_odd_wheel.cpp
   lifted_multicut_text.cpp
   lifted_multicut_h5.cpp
   lifted_multicut_h5_grid.cpp
   )


# note: HDF5 is only required for the opengm versions. Do not require HDF5 per default
foreach( source_file ${SOURCE_FILES} )
   string( REPLACE ".cpp" "" executable_file ${source_file} )
   add_executable( ${executable_file} ${source_file} ${headers} ${sources})
   target_link_libraries( ${executable_file} m stdc++ pthread ${HDF5_LIBRARIES})
endforeach( source_file ${SOURCE_FILES} )

