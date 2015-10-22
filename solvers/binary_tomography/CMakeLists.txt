
add_executable(binary_tomography binary_tomography.cpp ${headers} ${sources})
target_link_libraries(binary_tomography stdc++ Vc)