if(BUILD_TESTS)
   # tests that are compiled unconditionally
   SET(BASE_TEST_FILES 
      simplex.cpp
      simplex_marginalization.cpp
      min_cost_flow_cs2.cpp
      min_cost_flow.cpp
      )
   enable_testing()

   add_executable(test_main test_main.cpp ${headers} ${sources})
   target_link_libraries(test_main Vc m stdc++ pthread ${GUROBI_LIBRARIES}) # do zrobienia: include GUROBI in general?

   target_sources(test_main PUBLIC ${BASE_TEST_FILES})

   if(BUILD_GRAPHICAL_MODEL)
      target_sources(test_main PUBLIC graphical_model.cpp)
   endif()

   if(BUILD_GRAPH_MATCHING)
      target_sources(test_main PUBLIC graph_matching.cpp)
   endif()
   
   if(BUILD_MULTICUT)
      target_sources(test_main PUBLIC multicut.cpp)
   endif()

   if(BUILD_DISCRETE_TOMOGRAPHY)
      target_sources(test_main PUBLIC discrete_tomography.cpp)
      target_sources(test_main PUBLIC discrete_tomography_file_input.cpp)
      if(WITH_GUROBI)
         target_compile_definitions(test_main PUBLIC -DUSE_GUROBI)
         target_include_directories(test_main PUBLIC ${GUROBI_INCLUDE_DIR})
         target_sources(test_main PUBLIC discrete_tomography_gurobi.cpp)
      endif()
   endif()

   add_test(NAME LP_MP_TEST COMMAND test_main)

endif()

