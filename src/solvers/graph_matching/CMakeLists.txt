SET(SOURCE_FILES
   graph_matching_via_mp_left.cpp graph_matching_via_mp_right.cpp graph_matching_via_mp_both_sides.cpp 
   graph_matching_via_mcf_left.cpp graph_matching_via_mcf_right.cpp graph_matching_via_mcf_both_sides.cpp 
   graph_matching_via_gm_left.cpp graph_matching_via_gm_right.cpp 
   #graph_matching_via_mp_uai.cpp graph_matching_via_mcf_uai.cpp graph_matching_via_gm_uai.cpp
   hungarian_bp_left.cpp hungarian_bp_right.cpp hungarian_bp_both_sides.cpp 
   #hungarian_bp_uai.cpp

   # with frustrated cycle tightening
   graph_matching_via_mp_left_tightening.cpp graph_matching_via_mp_right_tightening.cpp graph_matching_via_mp_both_sides_tightening.cpp 
   graph_matching_via_mcf_left_tightening.cpp graph_matching_via_mcf_right_tightening.cpp graph_matching_via_mcf_both_sides_tightening.cpp 
   graph_matching_via_gm_left_tightening.cpp graph_matching_via_gm_right_tightening.cpp 
   #graph_matching_via_mp_uai_tightening.cpp graph_matching_via_mcf_uai_tightening.cpp graph_matching_via_gm_uai_tightening.cpp
   hungarian_bp_left_tightening.cpp hungarian_bp_right_tightening.cpp hungarian_bp_both_sides_tightening.cpp 
   #hungarian_bp_uai_tightening.cpp
   )

add_executable(convert_to_hdf5 convert_to_hdf5.cpp ${headers} ${sources})
target_link_libraries(convert_to_hdf5 m stdc++ pthread ${HDF5_LIBRARIES})

foreach( source_file ${SOURCE_FILES} )
   string( REPLACE ".cpp" "" executable_file ${source_file} )
   add_executable( ${executable_file} ${source_file} ${headers} ${sources})
   add_dependencies( ${executable_file} all_solvers)
   target_link_libraries( ${executable_file} m stdc++ pthread lgl)
endforeach( source_file ${SOURCE_FILES} )

if(WITH_GUROBI)
   SET(GUROBI_FILES
      # lp solver
      graph_matching_via_mp_left_tightening_gurobi.cpp graph_matching_via_mp_right_tightening_gurobi.cpp graph_matching_via_mp_both_sides_tightening_gurobi.cpp 
      graph_matching_via_mcf_left_tightening_gurobi.cpp graph_matching_via_mcf_right_tightening_gurobi.cpp graph_matching_via_mcf_both_sides_tightening_gurobi.cpp 
      graph_matching_via_gm_left_tightening_gurobi.cpp graph_matching_via_gm_right_tightening_gurobi.cpp 
      #graph_matching_via_mp_uai_tightening_gurobi.cpp graph_matching_via_mcf_uai_tightening_gurobi.cpp graph_matching_via_gm_uai_tightening_gurobi.cpp
      )

   foreach( source_file ${GUROBI_FILES} )
      string( REPLACE ".cpp" "" executable_file ${source_file} )
      add_executable( ${executable_file} ${source_file} ${headers} ${sources})
      add_dependencies( ${executable_file} all_solvers)
      target_link_libraries( ${executable_file} m stdc++ pthread ${GUROBI_LIBRARIES})

      target_compile_definitions(${executable_file} PUBLIC -DUSE_GUROBI)
      target_include_directories(${executable_file} PUBLIC ${GUROBI_INCLUDE_DIR})
   endforeach( source_file ${GUROBI_FILES} )
endif()

if(WITH_CPLEX)
        SET(CPLEX_FILES
      # lp solver
      graph_matching_via_mp_left_tightening_cplex.cpp graph_matching_via_mp_right_tightening_cplex.cpp graph_matching_via_mp_both_sides_tightening_cplex.cpp 
      graph_matching_via_mcf_left_tightening_cplex.cpp graph_matching_via_mcf_right_tightening_cplex.cpp graph_matching_via_mcf_both_sides_tightening_cplex.cpp 
      graph_matching_via_gm_left_tightening_cplex.cpp graph_matching_via_gm_right_tightening_cplex.cpp 
      #graph_matching_via_mp_uai_tightening_cplex.cpp graph_matching_via_mcf_uai_tightening_cplex.cpp graph_matching_via_gm_uai_tightening_cplex.cpp
      )

   foreach( source_file ${CPLEX_FILES} )
      string( REPLACE ".cpp" "" executable_file ${source_file} )
      add_executable( ${executable_file} ${source_file} ${headers} ${sources})
      add_dependencies( ${executable_file} all_solvers)
      target_link_libraries( ${executable_file} m stdc++ pthread ${CPLEX_LIBRARIES})

      target_compile_definitions(${executable_file} PUBLIC -DUSE_CPLEX)
      target_include_directories(${executable_file} PUBLIC ${CPLEX_INCLUDE_DIRS})
   endforeach( source_file ${CPLEX_FILES} )
endif()



