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
   target_link_libraries( ${executable_file} Vc m stdc++ pthread MinCost ${HDF5_LIBRARIES} )
endforeach( source_file ${SOURCE_FILES} )

