
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
target_link_libraries(test_conic_bundle LP_MP CONIC_BUNDLE lingeling)
add_test(test_conic_bundle test_conic_bundle)
