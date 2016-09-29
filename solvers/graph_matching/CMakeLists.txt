SET(SOURCE_FILES
   graph_matching_via_mp_left.cpp graph_matching_via_mp_right.cpp graph_matching_via_mp_both_sides.cpp 
   graph_matching_via_mcf_left.cpp graph_matching_via_mcf_right.cpp graph_matching_via_mcf_both_sides.cpp 
   graph_matching_via_gm_left.cpp graph_matching_via_gm_right.cpp 
   graph_matching_via_mp_uai.cpp graph_matching_via_mcf_uai.cpp graph_matching_via_gm_uai.cpp
   hungarian_bp_left.cpp hungarian_bp_right.cpp hungarian_bp_both_sides.cpp hungarian_bp_uai.cpp

   # with frustrated cycle tightening
   graph_matching_via_mp_left_tightening.cpp graph_matching_via_mp_right_tightening.cpp graph_matching_via_mp_both_sides_tightening.cpp 
   graph_matching_via_mcf_left_tightening.cpp graph_matching_via_mcf_right_tightening.cpp graph_matching_via_mcf_both_sides_tightening.cpp 
   graph_matching_via_gm_left_tightening.cpp graph_matching_via_gm_right_tightening.cpp 
   graph_matching_via_mp_uai_tightening.cpp graph_matching_via_mcf_uai_tightening.cpp graph_matching_via_gm_uai_tightening.cpp
   hungarian_bp_left_tightening.cpp hungarian_bp_right_tightening.cpp hungarian_bp_both_sides_tightening.cpp hungarian_bp_uai_tightening.cpp

   # lp solver
   )

foreach( source_file ${SOURCE_FILES} )
   string( REPLACE ".cpp" "" executable_file ${source_file} )
   add_executable( ${executable_file} ${source_file} ${headers} ${sources})
   target_link_libraries( ${executable_file} Vc m stdc++ pthread)
endforeach( source_file ${SOURCE_FILES} )

