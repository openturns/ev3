find_package (muParser)

include_directories (${CMAKE_SOURCE_DIR}/lib/src )
add_custom_target (tests)

if (MUPARSER_FOUND)
  include_directories (${MUPARSER_INCLUDE_DIRS})
  add_executable (t_ev3_std EXCLUDE_FROM_ALL t_ev3_std.cxx)
  target_link_libraries (t_ev3_std ev3 ${MUPARSER_LIBRARIES})
  add_dependencies (tests t_ev3_std)

  add_test (NAME t_ev3_std COMMAND t_ev3_std${CMAKE_EXECUTABLE_SUFFIX})

endif ()