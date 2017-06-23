add_executable( cell_tracking_mother_machine cell_tracking_mother_machine.cpp ${headers} ${sources} )
target_link_libraries( cell_tracking_mother_machine m stdc++ pthread lgl )

add_executable( cell_tracking_with_division_distance cell_tracking_with_division_distance.cpp ${headers} ${sources} )
target_link_libraries( cell_tracking_with_division_distance m stdc++ pthread lgl )

add_executable( cell_tracking cell_tracking.cpp ${headers} ${sources} )
target_link_libraries( cell_tracking m stdc++ pthread lgl )

add_executable( cell_tracking_duplicate_edges cell_tracking_duplicate_edges.cpp ${headers} ${sources} )
target_link_libraries( cell_tracking_duplicate_edges m stdc++ pthread lgl )

add_library( json_cpp jsoncpp.cpp )

add_executable( conservation_tracking conservation_tracking.cpp graphreader.cpp jsongraphreader.cpp jsoncpp.cpp ${headers} ${sources} )
target_link_libraries( conservation_tracking m stdc++ pthread lgl)
