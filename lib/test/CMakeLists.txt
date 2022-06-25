find_package (muParser)

include_directories (${CMAKE_SOURCE_DIR}/lib/src )

if (MUPARSER_FOUND AND BUILD_TESTING)
  include_directories (${MUPARSER_INCLUDE_DIRS})
  add_executable (t_ev3_std t_ev3_std.cxx)
  target_link_libraries (t_ev3_std ev3 ${MUPARSER_LIBRARIES})

  add_test (NAME t_ev3_std COMMAND t_ev3_std)

endif ()
