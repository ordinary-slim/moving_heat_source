cmake_minimum_required(VERSION 3.16)
project(MovingHeatSource LANGUAGES CXX)

FILE(GLOB_RECURSE mySources src/*.cpp)

add_subdirectory(external/pybind11)

set(CMAKE_DEBUG_POSTFIX _DEBUG)
pybind11_add_module(MovingHeatSource ${mySources})
set_target_properties(MovingHeatSource PROPERTIES
  DEBUG_POSTFIX ${CMAKE_DEBUG_POSTFIX})

target_include_directories(MovingHeatSource
  PRIVATE ${CMAKE_SOURCE_DIR}/includes)
