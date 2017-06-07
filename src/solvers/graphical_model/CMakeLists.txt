SET(SOURCE_FILES
   mplp_uai.cpp
   mplp_opengm.cpp
   srmp_uai.cpp
   srmp_opengm.cpp
   srmp_tightening_opengm.cpp
   )

foreach( source_file ${SOURCE_FILES} )
   string( REPLACE ".cpp" "" executable_file ${source_file} )
   add_executable( ${executable_file} ${source_file} ${headers} ${sources})
   add_dependencies( ${executable_file} all_solvers)
   target_link_libraries( ${executable_file} m stdc++ pthread ${HDF5_LIBRARIES} lgl)
endforeach( source_file ${SOURCE_FILES} )

