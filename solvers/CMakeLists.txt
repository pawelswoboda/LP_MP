if(BUILD_GRAPH_MATCHING)
   add_subdirectory(graph_matching)
endif()

if(BUILD_MULTICUT)
   add_subdirectory(multicut)
endif()

if(BUILD_GRAPHICAL_MODEL)
   add_subdirectory(graphical_model)
endif()

if(BUILD_BINARY_TOMOGRAPHY)
   add_subdirectory(binary_tomography)
 endif()

if(BUILD_TOMOGRAPHY)
   add_subdirectory(tomography)
endif()

