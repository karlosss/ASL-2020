cmake_minimum_required(VERSION 3.16)
project(asl_project C)

set(CMAKE_C_STANDARD 11)
set(CMAKE_CXX_FLAGS -std=gnu11)

add_executable(asl_project main.c measure.h common.h)
