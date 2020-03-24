#!/bin/bash
if [[ -z "$1" ]]
then
    echo "Usage: $0 <binary>" >&2
    exit 1
fi

if ! [[ -f "bin/$1" ]]
then
    echo "Binary $1 not found!" >&2
    exit 1
fi

number_of_cores=$(grep -c ^processor /proc/cpuinfo)
taskset -c $(( number_of_cores-1 )) "bin/$1"
