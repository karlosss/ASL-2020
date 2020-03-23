#!/bin/bash
number_of_cores=$(grep -c ^processor /proc/cpuinfo)
taskset -c $(( number_of_cores-1 )) ./asl_project
