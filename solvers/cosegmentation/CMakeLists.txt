#file(GLOB cosegmentation cosegmentation.cpp)
#add_executable(exec_cosegmentation exec_cosegmentation.cpp ${headers} ${sources} ${cosegmentation})
#target_link_libraries(exec_cosegmentation Vc m stdc++ pthread MinCost)

add_executable(cosegmentation cosegmentation.cpp ${headers} ${sources})
target_link_libraries(cosegmentation Vc m stdc++ pthread MinCost)

add_executable(cosegmentation_bins cosegmentation_bins.cpp ${headers} ${sources})
target_link_libraries(cosegmentation_bins Vc m stdc++ pthread MinCost)

