cmake_minimum_required(VERSION 3.10)

project(asl_project)
set(CMAKE_CXX_FLAGS "-O3 -std=c++17 -Wall -pedantic -Wextra")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY bin)

set(LIB lib/clock.h lib/common.h lib/utils.cpp lib/main.cpp lib/runtime.cpp base/base.cpp lib/modes.cpp lib/modes.h)
add_executable(validate base/validate.cpp ${LIB})

# ======================================================================================================================
# Register your executables here.

add_executable(no_opt opts/no_opt.cpp ${LIB})

