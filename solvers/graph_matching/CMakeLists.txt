SET(SOURCE_FILES
   graph_matching_via_mp_left.cpp graph_matching_via_mp_right.cpp graph_matching_via_mp_both_sides.cpp 
   graph_matching_via_mcf_left.cpp graph_matching_via_mcf_right.cpp graph_matching_via_mcf_both_sides.cpp 
   graph_matching_via_gm_left.cpp graph_matching_via_gm_right.cpp 
   graph_matching_via_mp_uai.cpp graph_matching_via_mcf_uai.cpp graph_matching_via_gm_uai.cpp
   graph_matching_via_mp_tightening_uai.cpp graph_matching_via_mcf_tightening_uai.cpp graph_matching_via_gm_tightening_uai.cpp
   )

foreach( source_file ${SOURCE_FILES} )
   string( REPLACE ".cpp" "" executable_file ${source_file} )
   add_executable( ${executable_file} ${source_file} ${headers} ${sources})
   target_link_libraries( ${executable_file} Vc m stdc++ pthread)
endforeach( source_file ${SOURCE_FILES} )

