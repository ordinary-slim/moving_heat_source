cmake_minimum_required(VERSION 3.16)
project(MovingHeatSource LANGUAGES CXX)

FILE(GLOB_RECURSE mySources cpp/*.cpp)

add_subdirectory(external/pybind11)

find_package(CGAL COMPONENTS Qt5)
add_definitions(-DCGAL_USE_BASIC_VIEWER)

#set(CMAKE_DEBUG_POSTFIX _DEBUG)
pybind11_add_module(MovingHeatSource ${mySources})
target_link_libraries(MovingHeatSource PUBLIC CGAL::CGAL_Qt5)
#set_target_properties(MovingHeatSource PROPERTIES
#DEBUG_POSTFIX ${CMAKE_DEBUG_POSTFIX})

target_include_directories(MovingHeatSource
  PRIVATE ${CMAKE_SOURCE_DIR}/includes)


# Tests
enable_testing()
add_subdirectory(tests)
