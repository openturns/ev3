cmake_minimum_required (VERSION 2.8)

set (CMAKE_BUILD_TYPE Release CACHE STRING "Build type")

project (ev3 CXX)
set (project_dir ${ev3_SOURCE_DIR})

list (APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)

if (MSVC)
  if (NOT ${CMAKE_BUILD_TYPE} STREQUAL "")
    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${project_dir}/${CMAKE_BUILD_TYPE}/bin)
    set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${project_dir}/${CMAKE_BUILD_TYPE}/lib)
    set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${project_dir}/${CMAKE_BUILD_TYPE}/lib)
    link_directories(${CMAKE_ARCHIVE_OUTPUT_DIRECTORY})
  else()
    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_DEBUG ${project_dir}/Debug/bin)
    set(CMAKE_LIBRARY_OUTPUT_DIRECTORY_DEBUG ${project_dir}/Debug/lib)
    set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY_DEBUG ${project_dir}/Debug/lib)
    link_directories(${CMAKE_ARCHIVE_OUTPUT_DIRECTORY_DEBUG})

    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_RELEASE ${project_dir}/Release/bin)
    set(CMAKE_LIBRARY_OUTPUT_DIRECTORY_RELEASE ${project_dir}/Release/lib)
    set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY_RELEASE ${project_dir}/Release/lib)
    link_directories(${CMAKE_ARCHIVE_OUTPUT_DIRECTORY_RELEASE})
  endif()
endif ()

option (BUILD_SHARED_LIBS "shared/static" ON)
option (BUILD_TESTING "tests" ON)

if (MSVC)
  find_package (Boost REQUIRED)
  include_directories (${Boost_INCLUDE_DIRS})
  list (APPEND EV3_DEFINITIONS -D_USE_MATH_DEFINES)
endif ()

add_definitions (${EV3_DEFINITIONS})

enable_testing()
add_subdirectory(lib)

