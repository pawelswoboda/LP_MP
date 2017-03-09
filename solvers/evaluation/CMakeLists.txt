function(DOWNLOAD_AND_UNZIP URL DOWNLOAD_TO EXTRACT_TO)
   # possibly download only to source dir. Then it need not be downloaded multiple times
   if(NOT EXISTS "${CMAKE_CURRENT_BINARY_DIR}/${DOWNLOAD_TO}")
      message("Download ${URL} to ${CMAKE_CURRENT_BINARY_DIR}/${DOWNLOAD_TO}")
      file(DOWNLOAD ${URL} "${CMAKE_CURRENT_BINARY_DIR}/${DOWNLOAD_TO}"
         )
      message(unzip "${CMAKE_CURRENT_BINARY_DIR}/${DOWNLOAD_TO}" -d "${CMAKE_CURRENT_BINARY_DIR}/${EXTRACT_TO}")
      execute_process(
         COMMAND unzip "${CMAKE_CURRENT_BINARY_DIR}/${DOWNLOAD_TO}" -d "${CMAKE_CURRENT_BINARY_DIR}/${EXTRACT_TO}"
         )
   endif()
endfunction(DOWNLOAD_AND_UNZIP)

if(BUILD_MULTICUT_EVALUATION)
   add_executable(multicut multicut.cpp ${headers} ${sources})
   target_link_libraries(multicut m stdc++ pthread sqlite3 ${HDF5_LIBRARIES})

   message("Download multicut datasets")
   set(MC_DIRECTORY "multicut_datasets")

   file(MAKE_DIRECTORY ${MC_DIRECTORY})
   DOWNLOAD_AND_UNZIP("http://hciweb2.iwr.uni-heidelberg.de/opengm/benchmarks/models/zip/knott-3d-150.zip" "${MC_DIRECTORY}/knott-3d-150.zip" ${MC_DIRECTORY})
   DOWNLOAD_AND_UNZIP("http://hciweb2.iwr.uni-heidelberg.de/opengm/benchmarks/models/zip/knott-3d-300.zip" "${MC_DIRECTORY}/knott-3d-300.zip" ${MC_DIRECTORY})
   DOWNLOAD_AND_UNZIP("http://hciweb2.iwr.uni-heidelberg.de/opengm/benchmarks/models/zip/knott-3d-450.zip" "${MC_DIRECTORY}/knott-3d-450.zip" ${MC_DIRECTORY})
   DOWNLOAD_AND_UNZIP("http://hciweb2.iwr.uni-heidelberg.de/opengm/benchmarks/models/zip/knott-3d-550.zip" "${MC_DIRECTORY}/knott-3d-550.zip" ${MC_DIRECTORY}/knott-3d-550/)
   DOWNLOAD_AND_UNZIP("http://hciweb2.iwr.uni-heidelberg.de/opengm/benchmarks/models/zip/modularity-clustering.zip" "${MC_DIRECTORY}/modularity-clustering.zip" "${MC_DIRECTORY}")
endif()

if(BUILD_GRAPH_MATCHING_EVALUATION)
   add_executable(graph_matching graph_matching.cpp ${headers} ${sources})
   target_link_libraries(graph_matching m stdc++ pthread sqlite3 ${HDF5_LIBRARIES})

   message("Download graph matching datasets")
   set(GM_DIRECTORY "graph_matching_datasets")

   file(MAKE_DIRECTORY "${GM_DIRECTORY}")
   DOWNLOAD_AND_UNZIP("http://www.cs.dartmouth.edu/~lorenzo/Papers/tkr_pami13_data.zip" "graph_matching_datasets/tkr_pami13_data.zip" "${GM_DIRECTORY}/tkr_pami13_data/")
endif()

if(BUILD_DISCRETE_TOMOGRAPHY_EVALUATION)
   add_executable(discrete_tomography_eval discrete_tomography_eval.cpp ${headers} ${sources})
   target_link_libraries(discrete_tomography_eval m stdc++ pthread pq pqxx ${HDF5_LIBRARIES})# ${GUROBI_LIBRARIES})

   #target_compile_definitions(discrete_tomography_eval PUBLIC -DUSE_GUROBI)
   #target_include_directories(discrete_tomography_eval PUBLIC ${GUROBI_INCLUDE_DIR})

   message("Download discrete tomography datasets")
   set(DT_DIRECTORY "discrete_tomography_datasets")

   DOWNLOAD_AND_UNZIP("https://datarep.app.ist.ac.at/46/1/discrete_tomography_synthetic.zip" "${DT_DIRECTORY}/discrete_tomography_synthetic.zip" "${DT_DIRECTORY}/")
endif()
