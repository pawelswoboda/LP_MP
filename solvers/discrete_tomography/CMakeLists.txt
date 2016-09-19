
add_executable(discrete_tomography discrete_tomography.cpp  ${headers} ${sources})
target_link_libraries(discrete_tomography Vc m stdc++ pthread)

if(WITH_GUROBI)
  add_executable(discrete_tomography_gurobi discrete_tomography_gurobi.cpp  ${headers} ${sources})
  target_link_libraries(discrete_tomography_gurobi Vc m stdc++ pthread ${GUROBI_LIBRARIES})

  target_compile_definitions(discrete_tomography_gurobi PUBLIC -DUSE_GUROBI)
  target_include_directories(discrete_tomography_gurobi PUBLIC ${GUROBI_INCLUDE_DIR})
endif()

if(WITH_CPLEX)
  add_executable(discrete_tomography_cplex discrete_tomography_cplex.cpp  ${headers} ${sources})
  target_link_libraries(discrete_tomography_cplex Vc m stdc++ pthread ${CPLEX_LIBRARIES})

  target_compile_definitions(discrete_tomography_cplex PUBLIC -DUSE_CPLEX)
  target_include_directories(discrete_tomography_cplex PUBLIC ${CPLEX_INCLUDE_DIRS})
endif()
