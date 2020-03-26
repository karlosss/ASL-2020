#!/bin/bash

VALIDATE_INPUT="100 100 100 5"
number_of_cores=$(grep -c ^processor /proc/cpuinfo)

if [[ "$1" = "h" || "$1" = "help" || "$1" = "-h" || "$1" = "--help" ]]
then
    echo "Usage:" >&2
    echo "$0 validate - run validation" >&2
    echo "$0 - run all binaries (specify input on stdin)" >&2
    echo "$0 BIN - run binary named BIN (specify input on stdin)" >&2
    echo "input is 4 integers - N, M, T, num_iters - separated by a space"
    exit 2
fi

if [[ "$1" = "validate" ]]
then
    if taskset -c $(( number_of_cores-1 )) "bin/$1" <<< "$VALIDATE_INPUT" > /dev/null
    then
        echo "Successfully validated."
        exit 0
    else
        echo "!!! Validation error !!!"
        exit 1
    fi
fi

echo "Awaiting input (N, M, T, num_iters)"
read in
cd bin || exit 3

if [[ -z "$1" ]]
then
    for i in *
    do
        if ! [[ "$i" = ".gitkeep" || "$i" = "validate" ]]
        then
            cyc=$(taskset -c $(( number_of_cores-1 )) "./$i" <<< "$in")
            if [[ $? -eq 2 ]]
            then
                echo "$i: $cyc"
            else
                printf '%.20s %s cycles\n' "$i.................................." "$cyc"
            fi
        fi
    done
else
    if ! [[ -f "$1" ]]
    then
        echo "Binary $1 not found!" >&2
        exit 1
    fi

    cyc=$(taskset -c $(( number_of_cores-1 )) "./$1" <<< "$in")
    if [[ $? -eq 2 ]]
    then
        echo "$1: $cyc"
    else
        printf '%.20s %s cycles\n' "$1.................................." "$cyc"
    fi
fi
