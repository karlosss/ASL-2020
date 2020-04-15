cmake_minimum_required(VERSION 3.10)

project(asl_project)
set(CMAKE_CXX_FLAGS "-O3 -std=c++17 -Wall -pedantic -Wextra")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY bin)
set(CMAKE_MODULE_PATH "${CMAKE_MODULE_PATH};${CMAKE_CURRENT_SOURCE_DIR}")
find_package(PAPI REQUIRED)
include_directories(${PAPI_INCLUDE_DIRS})
include_directories(lib)
add_library(LIBS lib/utils.cpp lib/main.cpp base/base.cpp lib/modes.cpp)
add_executable(validate base/validate.cpp)
target_link_libraries(validate LIBS ${PAPI_LIBRARIES})
# ======================================================================================================================

# Register your executables here.
add_executable(no_opt opts/no_opt.cpp)
target_link_libraries(no_opt LIBS ${PAPI_LIBRARIES})
