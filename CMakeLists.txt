cmake_minimum_required(VERSION 3.16)

project(asl_project)
set(CMAKE_CXX_FLAGS "-std=c++17 -Wall -pedantic -Wextra -O3")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY bin)

set(COMMON common/clock.h common/common.h common/generate.cpp common/main.cpp common/perf.cpp common/utils.cpp common/inline.h)

add_executable(test test/test.cpp ${COMMON})
add_executable(base opts/base.cpp ${COMMON})
