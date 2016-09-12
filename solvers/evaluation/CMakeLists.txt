SET(SOURCES graph_matching.cpp multicut.cpp)

function(DOWNLOAD_AND_UNZIP URL DOWNLOAD_TO EXTRACT_TO)
   if(NOT EXISTS ${DOWNLOAD_TO})
      message("Download ${URL}")
      file(DOWNLOAD "${URL}" "${DOWNLOAD_TO}")
      execute_process(
         COMMAND unzip  "${DOWNLOAD_TO}" -d "${EXTRACT_TO}"
         WORKING_DIRECTORY "${CMAKE_BINARY_DIR}"
         )
   endif()
endfunction(DOWNLOAD_AND_UNZIP)

# download test datasets
if(BUILD_MULTICUT)
   add_executable(multicut multicut.cpp ${headers} ${sources})
   target_link_libraries(multicut Vc m stdc++ pthread sqlite3 ${HDF5_LIBRARIES})

   message("Download multicut datasets")
   set(MC_DIRECTORY "solvers/evaluation/multicut_datasets")

   file(MAKE_DIRECTORY "${MC_DIRECTORY}")
   DOWNLOAD_AND_UNZIP("http://hciweb2.iwr.uni-heidelberg.de/opengm/benchmarks/models/zip/knott-3d-150.zip" "${MC_DIRECTORY}/knott-3d-150.zip" "${MC_DIRECTORY}")
   DOWNLOAD_AND_UNZIP("http://hciweb2.iwr.uni-heidelberg.de/opengm/benchmarks/models/zip/knott-3d-300.zip" "${MC_DIRECTORY}/knott-3d-300.zip" "${MC_DIRECTORY}")
   DOWNLOAD_AND_UNZIP("http://hciweb2.iwr.uni-heidelberg.de/opengm/benchmarks/models/zip/knott-3d-450.zip" "${MC_DIRECTORY}/knott-3d-450.zip" "${MC_DIRECTORY}")
   DOWNLOAD_AND_UNZIP("http://hciweb2.iwr.uni-heidelberg.de/opengm/benchmarks/models/zip/knott-3d-550.zip" "${MC_DIRECTORY}/knott-3d-550.zip" "${MC_DIRECTORY}/knott-3d-550/")
   DOWNLOAD_AND_UNZIP("http://hciweb2.iwr.uni-heidelberg.de/opengm/benchmarks/models/zip/modularity-clustering.zip" "${MC_DIRECTORY}/modularity-clustering.zip" "${MC_DIRECTORY}")
endif()

if(BUILD_GRAPH_MATCHING)
   add_executable(graph_matching graph_matching.cpp ${headers} ${sources})
   target_link_libraries(graph_matching Vc m stdc++ pthread sqlite3 ${HDF5_LIBRARIES})

   message("Download graph matching datasets")
   set(GM_DIRECTORY "solvers/evaluation/graph_matching_datasets")

   file(MAKE_DIRECTORY "${GM_DIRECTORY}")
   DOWNLOAD_AND_UNZIP("http://www.cs.dartmouth.edu/~lorenzo/Papers/tkr_pami13_data.zip" "${GM_DIRECTORY}/tkr_pami13_data.zip" "${GM_DIRECTORY}/tkr_pami13_data/")
endif()

