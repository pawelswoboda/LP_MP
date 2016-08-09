
add_executable(discrete_tomography discrete_tomography.cpp  ${headers} ${sources})
target_link_libraries(discrete_tomography Vc m stdc++ pthread)

add_definitions(-DUSE_GUROBI)
include_directories(${GUROBI_INCLUDE_DIR})

add_executable(discrete_tomography_lp discrete_tomography_lp.cpp  ${headers} ${sources})
target_link_libraries(discrete_tomography_lp Vc m stdc++ pthread ${GUROBI_LIBRARIES})


