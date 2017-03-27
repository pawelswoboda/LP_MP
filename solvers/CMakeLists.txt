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

if(BUILD_CELL_TRACKING)
  add_subdirectory(cell_tracking)
endif()

if(BUILD_MORAL_LINEAGE_TRACING)
  add_subdirectory(moral_lineage_tracing)
endif()

add_subdirectory(evaluation)

