cmake_minimum_required(VERSION 3.10)

project(asl_project)
set(CMAKE_CXX_FLAGS "-std=c++17 -Wall -pedantic -Wextra ${FLAGS}")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY bin/${COMPILER}_${UNDERSCORED_FLAGS})
set(CMAKE_CXX_COMPILER ${COMPILER})
set(CMAKE_MODULE_PATH "${CMAKE_MODULE_PATH};${CMAKE_CURRENT_SOURCE_DIR}")
set(PAPI_INCLUDE_DIRS ${CMAKE_MODULE_PATH}/papi/include)
set(PAPI_LIBRARIES ${CMAKE_MODULE_PATH}/papi/lib/libpapi.so.6.0.0.1)
include_directories(${PAPI_INCLUDE_DIRS})
include_directories(lib)
add_library(LIBS lib/utils.cpp lib/main.cpp base/base.cpp lib/modes.cpp)
# ======================================================================================================================
# For each executable, add the following two lines:

# Manually verify the base algorithm
add_executable(base_validate base/base_validate.cpp)
target_link_libraries(base_validate LIBS ${PAPI_LIBRARIES})

# Verify the base algorithm using assertions to check properties
add_executable(validate base/validate.cpp)
target_link_libraries(validate LIBS ${PAPI_LIBRARIES})

add_executable(no_opt opts/no_opt.cpp)
target_link_libraries(no_opt LIBS ${PAPI_LIBRARIES})

add_executable(opt1 opts/opt1.cpp)
target_link_libraries(opt1 LIBS ${PAPI_LIBRARIES})
