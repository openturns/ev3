find_package (muParser)

include_directories (${CMAKE_SOURCE_DIR}/lib/src)
include_directories (${CMAKE_CURRENT_SOURCE_DIR})

if (BUILD_TESTING)
  add_executable (t_ev3_std t_ev3_std.cxx)
  target_link_libraries (t_ev3_std ev3)
  if (CMAKE_CXX_COMPILER_ID MATCHES "MSVC")
    set_source_files_properties (t_ev3_std.cxx PROPERTIES COMPILE_FLAGS "/bigobj")
  endif ()
  add_test (NAME t_ev3_std COMMAND t_ev3_std)
endif ()
