if(BUILD_TESTS)
   SET(TEST_FILES 
      simplex.cpp
      simplex_marginalization.cpp
      multicut.cpp
      )
   enable_testing()


   add_executable(test_main test_main.cpp ${headers} ${sources})
   target_link_libraries(test_main Vc m stdc++ pthread )

   foreach( source_file ${TEST_FILES} )
      string( REPLACE "_unit_test.cpp" "" object ${source_file} )
      add_library(${object} SHARED ${source_file}) # for some reason only shared library is linked to test_main
      target_link_libraries(test_main ${object})
   endforeach( source_file ${TEST_FILES} )

   add_test(NAME LP_MP_TEST COMMAND test_main)
endif()

