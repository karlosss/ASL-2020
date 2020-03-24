cmake_minimum_required(VERSION 3.16)
project(asl_project)

set(CMAKE_CXX_FLAGS -std=c++11)

add_executable(measure_runtime main.cpp utils.cpp base.cpp measure.h common.h )

project(test_base)
set(test_base_SOURCES test_base.cpp base.cpp utils.cpp)
set(test_base_HEADERS base.h common.h)
add_definitions(-std=c++11 -Wall -Wextra -O3)
add_executable(test_base ${test_base_SOURCES} ${test_base_HEADERS})