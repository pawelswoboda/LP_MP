SET(SOURCES graph_matching.cpp multicut.cpp)

# download test datasets

if(BUILD_MULTICUT)
   message("Download multicut datasets")
   set(URL "http://hciweb2.iwr.uni-heidelberg.de/opengm/benchmarks/models/zip/knott-3d-150.zip")
   set(DIRECTORY "solvers/evaluation/multicut_datasets")
   set(DOWNLOAD_PATH "${DIRECTORY}/knott-3d-150.zip")
   set(EXTRACTED_FILE "${DIRECTORY}/knott-3d-150/")
   file(MAKE_DIRECTORY "${DIRECTORY}")
   if(NOT EXISTS "${DOWNLOAD_PATH}")
      message("actually download multicut datasets")
      file(DOWNLOAD "${URL}" "${DOWNLOAD_PATH}")
      add_custom_command(
         OUTPUT "${EXTRACTED_FILE}"
         COMMAND "unzip ${DOWNLOAD_PATH}"
         DEPENDS "${DOWNLOAD_PATH}"
         COMMENT "unpack multicut datasets"
         )
   endif()
endif()


# do zrobienia: Only build evaluation for the respective solvers if the respective solvers are enabled as well

# write a function for this
foreach( source_file ${SOURCES} )
   string( REPLACE ".cpp" "" executable_file ${source_file} )
   add_executable( ${executable_file} ${source_file} ${headers} ${sources})
   target_link_libraries( ${executable_file} Vc m stdc++ pthread sqlite3 ${HDF5_LIBRARIES})
endforeach( source_file ${SOURCES} )



