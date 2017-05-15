SET(SOURCE_FILES 
   max_cut_cycle.cpp
   max_cut_cycle_odd_bicycle_wheel.cpp
   QUBO_cycle.cpp
   QUBO_cycle_odd_bicycle_wheel.cpp
   )


# note: HDF5 is only required for the opengm versions. Do not require HDF5 per default
foreach( source_file ${SOURCE_FILES} )
   string( REPLACE ".cpp" "" executable_file ${source_file} )
   add_executable( ${executable_file} ${source_file} ${headers} ${sources})
   target_link_libraries( ${executable_file} m stdc++ pthread ${HDF5_LIBRARIES} ${HDF5_CXX_LIBRARIES} lgl)
endforeach( source_file ${SOURCE_FILES} )

