if(BUILD_GRAPH_MATCHING)
   add_subdirectory(graph_matching)
endif()

if(BUILD_MULTICUT)
   add_subdirectory(multicut)
endif()

if(BUILD_GRAPHICAL_MODEL)
   add_subdirectory(graphical_model)
endif()

if(BUILD_DISCRETE_TOMOGRAPHY)
   add_subdirectory(discrete_tomography)
endif()

add_subdirectory(evaluation)

