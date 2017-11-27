add_executable(vector vector.cpp ${headers})
target_link_libraries( vector LP_MP m stdc++ pthread )
add_executable(serialization serialization.cpp ${headers})
target_link_libraries( serialization LP_MP m stdc++ pthread )
