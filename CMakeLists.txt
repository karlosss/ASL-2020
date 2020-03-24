cmake_minimum_required(VERSION 3.16)

project(asl_project)
set(CMAKE_CXX_FLAGS "-std=c++11 -Wall -pedantic -Wextra -O3")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY bin)

set(test_base_SOURCES test_base.cpp base.cpp utils.cpp)
set(test_base_HEADERS base.h common.h)

add_executable(test ${test_base_SOURCES} ${test_base_HEADERS})
add_executable(base main.cpp utils.cpp base.cpp measure.h common.h)
