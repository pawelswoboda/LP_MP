
add_executable(discrete_tomography discrete_tomography.cpp  ${headers} ${sources})
target_link_libraries(discrete_tomography Vc m stdc++ pthread)

if(WITH_GUROBI)
  add_executable(discrete_tomography_gurobi discrete_tomography_gurobi.cpp  ${headers} ${sources})
  target_link_libraries(discrete_tomography_gurobi Vc m stdc++ pthread ${GUROBI_LIBRARIES})

  target_compile_definitions(discrete_tomography_gurobi PUBLIC -DUSE_GUROBI)
  target_include_directories(discrete_tomography_gurobi PUBLIC ${GUROBI_INCLUDE_DIR})


  add_executable(discrete_tomography_standard_relax_gurobi discrete_tomography_standard_relax_gurobi.cpp  ${headers} ${sources})
  target_link_libraries(discrete_tomography_standard_relax_gurobi Vc m stdc++ pthread ${GUROBI_LIBRARIES})

  target_compile_definitions(discrete_tomography_standard_relax_gurobi PUBLIC -DUSE_GUROBI)
  target_include_directories(discrete_tomography_standard_relax_gurobi PUBLIC ${GUROBI_INCLUDE_DIR})
endif()

if(WITH_CPLEX)
  add_executable(discrete_tomography_cplex discrete_tomography_cplex.cpp  ${headers} ${sources})
  target_link_libraries(discrete_tomography_cplex Vc m stdc++ pthread ${CPLEX_LIBRARIES})

  target_compile_definitions(discrete_tomography_cplex PUBLIC -DUSE_CPLEX)
  target_include_directories(discrete_tomography_cplex PUBLIC ${CPLEX_INCLUDE_DIRS})

  INSTALL(TARGETS discrete_tomography_cplex
    ARCHIVE DESTINATION ${PROJECT_BINARY_DIR}/install/discrete_tomography_cplex/lib
    LIBRARY DESTINATION ${PROJECT_BINARY_DIR}/install/discrete_tomography_cplex/lib
    RUNTIME DESTINATION ${PROJECT_BINARY_DIR}/install/discrete_tomography_cplex/bin)  
  
  SET(_files ${Vc_LIBRARIES} ${CPLEX_LIBRARIES})
  SET(_resolvedFiles "")
  foreach (_file ${_files})
    get_filename_component(_resolvedFile "${_file}" REALPATH)
    list (APPEND _resolvedFiles "${_resolvedFile}")
  endforeach()

  INSTALL(FILES ${_resolvedFiles} DESTINATION ${PROJECT_BINARY_DIR}/install/discrete_tomography_cplex/lib)
endif()


