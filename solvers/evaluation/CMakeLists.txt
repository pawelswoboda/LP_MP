SET(SOURCES graph_matching.cpp multicut.cpp)
# do zrobienia: Only build evaluation for the respective solvers if the respective solvers are enabled as well

# write a function for this
foreach( source_file ${SOURCES} )
   string( REPLACE ".cpp" "" executable_file ${source_file} )
   add_executable( ${executable_file} ${source_file} ${headers} ${sources})
   target_link_libraries( ${executable_file} Vc m stdc++ pthread MinCost sqlite3 ${HDF5_LIBRARIES})
endforeach( source_file ${SOURCES} )


