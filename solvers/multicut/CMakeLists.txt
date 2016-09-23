SET(SOURCE_FILES 
   multicut_srmp.cpp 
   multicut_opengm_srmp_cycle.cpp 
   multicut_opengm_srmp_cycle_odd_wheel.cpp
   lifted_multicut_text.cpp
   lifted_multicut_h5.cpp
   )


# note: HDF5 is only required for the opengm versions. Do not require HDF5 per default
foreach( source_file ${SOURCE_FILES} )
   string( REPLACE ".cpp" "" executable_file ${source_file} )
   add_executable( ${executable_file} ${source_file} ${headers} ${sources})
   add_dependencies(${executable_file} Andres)
   add_dependencies(${executable_file} Vc)
   target_link_libraries( ${executable_file} Vc m stdc++ pthread ${HDF5_LIBRARIES} ncurses)
endforeach( source_file ${SOURCE_FILES} )

