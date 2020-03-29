#!/bin/bash

number_of_cores=$(grep -c ^processor /proc/cpuinfo)
in=""
fmt="%-25s%-16s%-16s\n"

print_usage () {
    echo "Usage:" >&2
    echo "$0 validate - run validation (specify input on stdin)" >&2
    echo "$0 debug BIN - run binary BIN in debug mode (no runtime measurement, specify input on stdin)" >&2
    echo "$0 runtime - run all binaries in runtime mode (specify input on stdin)" >&2
    echo "$0 runtime BIN - run binary BIN in runtime mode (specify input on stdin)" >&2
    echo "input is 4 integers - N, M, T, num_iters - separated by a space"
    exit 2
}

read_input () {
    echo "Awaiting input (N, M, T, num_iters)"
    read -r in
}

print_not_found () {
    echo "Binary $1 not found!" >&2
    exit 1
}

run_binary () {
    cyc=$(taskset -c $(( number_of_cores-1 )) ./"$1" runtime <<< "$in")
    flops=$(taskset -c $(( number_of_cores-1 )) ./"$1" flops <<< "$in")
    flops_cyc=$(bc -l <<< "$flops/$cyc")

    sci_cyc=$(printf "%.5g" $cyc)
    sci_fc=$(printf "%.5g" $flops_cyc)


    if [[ $? -eq 2 ]]
    then
        echo "$i: $cyc"
    else
        printf "$fmt" "$1" "$sci_cyc" "$sci_fc"
    fi
}

if [[ "$1" = "h" || "$1" = "help" || "$1" = "-h" || "$1" = "--help" || -z "$1" ]]
then
    print_usage
fi

cd bin || exit 3

if [[ "$1" = "validate" ]]
then
    read_input
    if taskset -c $(( number_of_cores-1 )) ./validate debug <<< "$in" > /dev/null
    then
        echo "Successfully validated."
        exit 0
    else
        echo "!!! Validation error !!!"
        exit 1
    fi
fi

if [[ "$1" = "debug" ]]
then
    if [[ -z "$2" ]]
    then
        print_usage
    else
        read_input
        if ! [[ -f "$2" ]]
        then
            print_not_found "$2"
        fi
        taskset -c $(( number_of_cores-1 )) ./"$2" <<< "$in"
        exit $?
    fi
fi

if [[ "$1" = "runtime" ]]
then
    read_input
    printf "$fmt" BINARY CYCLES FLOPS/CYCLE
    if [[ -z "$2" ]]
    then
        for i in *
        do
            if ! [[ "$i" = ".gitkeep" || "$i" = "validate" ]]
            then
                run_binary "$i"
            fi
        done
    else
        if ! [[ -f "$2" ]]
        then
            print_not_found "$2"
        fi
        run_binary "$2"
    fi
    exit 0
fi
