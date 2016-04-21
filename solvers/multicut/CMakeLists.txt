SET(SOURCE_FILES 
   multicut_srmp.cpp multicut_mplp.cpp
   multicut_opengm_srmp.cpp multicut_opengm_mplp.cpp
   )

# note: HDF5 is only required for the opengm versions. Do not require HDF5 per default
foreach( source_file ${SOURCE_FILES} )
   string( REPLACE ".cpp" "" executable_file ${source_file} )
   add_executable( ${executable_file} ${source_file} ${headers} ${sources})
   target_link_libraries( ${executable_file} Vc m stdc++ pthread ${HDF5_LIBRARIES})
endforeach( source_file ${SOURCE_FILES} )

