
add_executable(discrete_tomography discrete_tomography.cpp  ${headers} ${sources})
target_link_libraries(discrete_tomography m stdc++ pthread)

if(WITH_GUROBI)

        SET(SOURCE_FILES
                discrete_tomography_gurobi.cpp
                discrete_tomography_standard_relax_gurobi.cpp
                discrete_tomography_combined_gurobi.cpp)

        foreach( source_file ${SOURCE_FILES} )
                string( REPLACE ".cpp" "" executable_file ${source_file} )
                add_executable( ${executable_file} ${source_file} ${headers} ${sources})
                target_link_libraries( ${executable_file} m stdc++ pthread ${GUROBI_LIBRARIES})

                target_compile_definitions(${executable_file} PUBLIC -DUSE_GUROBI)
                target_include_directories(${executable_file} PUBLIC ${GUROBI_INCLUDE_DIR})
        endforeach( source_file ${SOURCE_FILES} )


endif()

if(WITH_LOCALSOLVER)
  add_executable(discrete_tomography_localsolver discrete_tomography_localsolver.cpp  ${headers} ${sources})
  target_link_libraries(discrete_tomography_localsolver m stdc++ pthread ${LOCALSOLVER_LIB})

  target_compile_definitions(discrete_tomography_localsolver PUBLIC -DUSE_LOCALSOLVER)
  target_include_directories(discrete_tomography_localsolver PUBLIC ${LOCALSOLVER_INCLUDE_DIR})
endif()

if(WITH_CPLEX)
  add_executable(discrete_tomography_cplex discrete_tomography_cplex.cpp  ${headers} ${sources})
  target_link_libraries(discrete_tomography_cplex m stdc++ pthread ${CPLEX_LIBRARIES})

  target_compile_definitions(discrete_tomography_cplex PUBLIC -DUSE_CPLEX)
  target_include_directories(discrete_tomography_cplex PUBLIC ${CPLEX_INCLUDE_DIRS})

  INSTALL(TARGETS discrete_tomography_cplex
    ARCHIVE DESTINATION ${PROJECT_BINARY_DIR}/install/discrete_tomography_cplex/lib
    LIBRARY DESTINATION ${PROJECT_BINARY_DIR}/install/discrete_tomography_cplex/lib
    RUNTIME DESTINATION ${PROJECT_BINARY_DIR}/install/discrete_tomography_cplex/bin)  
  
  SET(_files ${CPLEX_LIBRARIES}) 
  SET(_resolvedFiles "")
  foreach (_file ${_files})
    get_filename_component(_resolvedFile "${_file}" REALPATH)
    list (APPEND _resolvedFiles "${_resolvedFile}")
  endforeach()

  INSTALL(FILES ${_resolvedFiles} DESTINATION ${PROJECT_BINARY_DIR}/install/discrete_tomography_cplex/lib)
endif()


add_executable(discrete_tomography_fw discrete_tomography_fw.cpp ${headers} ${sources})
target_link_libraries(discrete_tomography_fw m stdc++ pthread SVM)
