if(BUILD_TESTS)
   # tests that are compiled unconditionally
   SET(BASE_TEST_FILES 
      simplex.cpp
      potts_factor.cpp
      simplex_marginalization.cpp
      serialization.cpp
      vector.cpp
      labeling_factor_messages.cpp
      #min_cost_flow.cpp
      #min_conv.cpp
      shortest_path.cpp
      cycle_inequalities.cpp
      graphical_model.cpp
      #graph_matching.cpp
      multicut.cpp
      discrete_tomography_chain.cpp
      sat_solver.cpp
      )
   enable_testing()

   add_library(Catch INTERFACE)
   add_executable(test_main test_main.cpp ${headers} ${sources} ${BASE_TEST_FILES})
   target_link_libraries(test_main m stdc++ pthread lgl ${HDF5_LIBRARIES}) # do zrobienia: include GUROBI in general?
   target_link_libraries(test_main Catch)


   #add_executable(test_main test_main.cpp ${headers} ${sources})

   #target_sources(test_main PUBLIC ${BASE_TEST_FILES})

   if(BUILD_GRAPHICAL_MODEL)
      #target_sources(test_main PUBLIC graphical_model.cpp)
   endif()

   if(BUILD_GRAPH_MATCHING)
      #target_sources(test_main PUBLIC graph_matching.cpp)
   endif()
   
   if(BUILD_MULTICUT)
      #target_sources(test_main PUBLIC multicut.cpp)
   endif()

   if(BUILD_DISCRETE_TOMOGRAPHY)
     #target_sources(test_main PUBLIC discrete_tomography.cpp)
      #target_sources(test_main PUBLIC discrete_tomography_file_input.cpp)
      if(WITH_GUROBI)
      #   target_sources(test_main PUBLIC discrete_tomography_gurobi.cpp)
      endif()
   endif()

   if(WITH_GUROBI)
      #target_compile_definitions(test_main PUBLIC -DUSE_GUROBI)
      #target_include_directories(test_main PUBLIC ${GUROBI_INCLUDE_DIR})
   endif()

   add_test(NAME LP_MP_TEST COMMAND test_main)

endif()

