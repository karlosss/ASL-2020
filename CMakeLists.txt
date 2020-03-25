cmake_minimum_required(VERSION 3.16)

project(asl_project)
set(CMAKE_CXX_FLAGS "-std=c++17 -Wall -pedantic -Wextra -O3")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY bin)

set(LIB lib/clock.h lib/common.h lib/utils.cpp lib/main.cpp lib/perf.cpp base/base.cpp lib/modes.cpp lib/modes.h)

add_executable(no_opt opts/no_opt.cpp ${LIB})
add_executable(validate base/validate.cpp ${LIB})
