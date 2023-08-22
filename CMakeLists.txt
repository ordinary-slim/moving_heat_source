cmake_minimum_required(VERSION 3.16)
project(MovingHeatSource LANGUAGES CXX)

FILE(GLOB_RECURSE mySources cpp/*.cpp)

add_subdirectory(external/pybind11)

find_package(CGAL COMPONENTS Qt5)
add_definitions(-DCGAL_USE_BASIC_VIEWER)

pybind11_add_module(cpp ${mySources})

target_include_directories(cpp PUBLIC ${PROJECT_SOURCE_DIR}/external/pybind11/include)
target_include_directories(cpp PUBLIC ${PROJECT_SOURCE_DIR}/external/eigen)
target_link_libraries(cpp PUBLIC CGAL::CGAL_Qt5)
