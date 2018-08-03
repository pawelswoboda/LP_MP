
add_executable(vector vector.cpp ${headers})
target_link_libraries( vector LP_MP m stdc++ pthread )
add_test( vector vector )

add_executable(serialization serialization.cpp ${headers})
target_link_libraries( serialization LP_MP m stdc++ pthread )
add_test( serialization serialization )

add_executable(test_model test_model.cpp ${headers})
target_link_libraries(test_model LP_MP DD_ILP lingeling)
add_test( test_model test_model )

add_executable(test_FWMAP test_FWMAP.cpp)
target_link_libraries(test_FWMAP LP_MP FW-MAP lingeling)
add_test(test_FWMAP test_FWMAP)

add_executable(test_conic_bundle test_conic_bundle.cpp)
target_link_libraries(test_conic_bundle cb LP_MP lingeling)
add_test(test_conic_bundle test_conic_bundle)

#add_executable(two_dimensional_variable_array two_dimensional_variable_array.cpp)
#target_link_libraries(two_dimensional_variable_array LP_MP)
#add_test(two_dimensional_variable_array two_dimensional_variable_array) 
