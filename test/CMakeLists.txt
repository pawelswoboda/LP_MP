if(BUILD_TESTS)
   SET(TEST_FILES 
      #simplex.cpp
      #simplex_marginalization.cpp
      #multicut.cpp
      #discrete_tomography.cpp
      #graphical_model.cpp
      graph_matching.cpp
      )
   enable_testing()

   add_executable(test_main test_main.cpp ${headers} ${sources} ${TEST_FILES})
   target_link_libraries(test_main Vc m stdc++ pthread )

   add_test(NAME LP_MP_TEST COMMAND test_main)

   #add_subdirectory(discrete_tomography)
endif()

