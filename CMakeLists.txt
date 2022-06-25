cmake_minimum_required (VERSION 3.2)

set (CMAKE_BUILD_TYPE Release CACHE STRING "Build type")

project (ev3 CXX)
set (project_dir ${ev3_SOURCE_DIR})

list (APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)

#set (CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

option (BUILD_SHARED_LIBS "shared/static" ON)
option (BUILD_TESTING "tests" ON)

if (MSVC)
  list (APPEND EV3_DEFINITIONS -D_USE_MATH_DEFINES)
endif ()

add_definitions (${EV3_DEFINITIONS})

enable_testing()
add_subdirectory(lib)

