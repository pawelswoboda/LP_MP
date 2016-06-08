add_executable(discrete_tomography discrete_tomography.cpp  ${headers} ${sources})
target_link_libraries(discrete_tomography Vc m stdc++ pthread)