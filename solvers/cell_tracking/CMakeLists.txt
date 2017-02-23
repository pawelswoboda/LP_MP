add_executable( cell_tracking cell_tracking.cpp ${headers} ${sources} )
target_link_libraries( cell_tracking Vc m stdc++ pthread cryptominisat5 )
#target_link_libraries( cell_tracking Vc m stdc++ pthread glucose )

add_library( json_cpp jsoncpp.cpp )

add_executable( conservation_tracking conservation_tracking.cpp graphreader.cpp jsongraphreader.cpp jsoncpp.cpp ${headers} ${sources} )
target_link_libraries( conservation_tracking Vc m stdc++ pthread )
