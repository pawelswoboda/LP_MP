add_executable(test_two_dimensional_variable_array test_two_dimensional_variable_array.cpp)
target_link_libraries(test_two_dimensional_variable_array  LP_MP m stdc++)
add_test(test_two_dimensional_variable_array  test_two_dimensional_variable_array) 

add_executable(test_three_dimensional_variable_array test_three_dimensional_variable_array.cpp)
target_link_libraries(test_three_dimensional_variable_array  LP_MP m stdc++)
add_test(test_three_dimensional_variable_array  test_three_dimensional_variable_array) 

add_executable(vector vector.cpp)
target_link_libraries( vector LP_MP m stdc++ pthread )
add_test( vector vector )

add_executable(serialization serialization.cpp)
target_link_libraries( serialization LP_MP m stdc++ pthread )
add_test( serialization serialization )

add_executable(test_model test_model.cpp)
target_link_libraries(test_model LP_MP DD_ILP lingeling)
add_test( test_model test_model )

add_executable(test_FWMAP test_FWMAP.cpp)
target_link_libraries(test_FWMAP LP_MP FW-MAP lingeling)
add_test(test_FWMAP test_FWMAP)

add_executable(test_conic_bundle test_conic_bundle.cpp)
target_link_libraries(test_conic_bundle cb LP_MP lingeling)
add_test(test_conic_bundle test_conic_bundle)

add_executable(graph_test graph_test.cpp)
target_link_libraries(graph_test LP_MP)
add_test(graph_test graph_test)
