add_executable(discrete_tomography discrete_tomography.cpp  ${headers} ${sources})
target_link_libraries(discrete_tomography Vc m stdc++ pthread)

#add_executable(minConvTest minConvTest.cxx  ${headers} ${sources})
#target_link_libraries(minConvTest Vc m stdc++ pthread)
