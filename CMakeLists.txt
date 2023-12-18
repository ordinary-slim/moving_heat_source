cmake_minimum_required(VERSION 3.13)
project(MovingHeatSource LANGUAGES CXX)

# Build compiler_commands.json for LSPs
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
add_custom_target(
    copy-compile-commands ALL
    ${CMAKE_COMMAND} -E copy_if_different
        ${CMAKE_BINARY_DIR}/compile_commands.json
        ${CMAKE_CURRENT_LIST_DIR}
    )

FILE(GLOB_RECURSE mySources cpp/*.cpp cpp/linearAlgebra/*.cpp)

add_subdirectory(external/pybind11)

set(CGAL_DIR ${PROJECT_SOURCE_DIR}/external/cgal)
find_package(CGAL REQUIRED)
add_definitions(-DCGAL_USE_BASIC_VIEWER)

# Profiling with Valgrind
if($ENV{Valgrind})
  set(CMAKE_BUILD_TYPE Debug)
  add_definitions(-DCGAL_DISABLE_ROUNDING_MATH_CHECK)
endif()

pybind11_add_module(cpp ${mySources})

set_property(TARGET cpp PROPERTY CXX_STANDARD 17)
target_include_directories(cpp PUBLIC ${PROJECT_SOURCE_DIR}/external/pybind11/include)
target_include_directories(cpp PUBLIC ${PROJECT_SOURCE_DIR}/external/eigen)
target_link_libraries(cpp PUBLIC CGAL::CGAL)

# MATLAB
if($ENV{hasMATLAB})
  target_compile_definitions(cpp PUBLIC hasMATLAB)
  set(matlabroot "/usr/local/MATLAB/R2023a")
  set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
  target_include_directories(cpp PUBLIC "${matlabroot}/extern/include/")
  target_link_directories(cpp PUBLIC "${matlabroot}/extern/bin/glnxa64/")
  target_link_libraries(cpp PUBLIC MatlabDataArray)
  target_link_libraries(cpp PUBLIC MatlabEngine)
  SET(CMAKE_EXE_LINKER_FLAGS  "${CMAKE_EXE_LINKER_FLAGS} -pthread")
endif()
