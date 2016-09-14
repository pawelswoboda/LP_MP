
add_executable(discrete_tomography discrete_tomography.cpp  ${headers} ${sources})
target_link_libraries(discrete_tomography Vc m stdc++ pthread)

if(WITH_GUROBI)
   add_executable(discrete_tomography_lp discrete_tomography_lp.cpp  ${headers} ${sources})
   target_link_libraries(discrete_tomography_lp Vc m stdc++ pthread ${GUROBI_LIBRARIES})
   target_compile_definitions(discrete_tomography_lp PUBLIC -DUSE_GUROBI)
   target_include_directories(discrete_tomography_lp PUBLIC ${GUROBI_INCLUDE_DIR})
endif()

