cmake_minimum_required(VERSION 3.16)
project(MovingHeatSource LANGUAGES CXX)

FILE(GLOB_RECURSE mySources src/*.cpp)

add_subdirectory(external/pybind11)

pybind11_add_module(MovingHeatSource ${mySources})

target_include_directories(MovingHeatSource
  PRIVATE ${CMAKE_SOURCE_DIR}/includes)
