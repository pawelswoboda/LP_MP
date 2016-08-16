
add_executable(discrete_tomography_debug discrete_tomography.cpp  ${headers} ${sources})
target_link_libraries(discrete_tomography_debug Vc m stdc++ pthread)

add_custom_command(TARGET discrete_tomography_debug POST_BUILD
  COMMAND ${CMAKE_COMMAND} -E copy
  ${CMAKE_SOURCE_DIR}/test/discrete_tomography/RunScript.sh $<TARGET_FILE_DIR:discrete_tomography_debug>/)
