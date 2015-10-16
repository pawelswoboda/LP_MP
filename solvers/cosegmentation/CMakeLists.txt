file(GLOB cosegmentation cosegmentation.cpp)
add_executable(exec_cosegmentation exec_cosegmentation.cpp ${headers} ${sources} ${cosegmentation})
target_link_libraries(exec_cosegmentation Vc m stdc++ pthread MinCost)

