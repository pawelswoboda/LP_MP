file(GLOB graph_matching graph_matching.cpp)
add_executable(exec_graph_matching exec_graph_matching.cpp ${headers} ${sources} ${graph_matching})
target_link_libraries(exec_graph_matching Vc m stdc++ pthread MinCost)
