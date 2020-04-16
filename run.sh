#!/bin/bash

print_usage () {
    echo "Usage:" >&2
    echo "$0 validate - run validation of base implementation (specify input on stdin)" >&2
    echo "$0 BIN - execute binary BIN in run mode and log measurements(specify input on stdin)" >&2
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

parse_json () {
    log=$(echo logs/papi_hl_output/*)
    # parse here
}

if [[ "$1" = "h" || "$1" = "help" || "$1" = "-h" || "$1" = "--help" || -z "$1" ]]
then
    print_usage
fi

perf_paranoid=$(cat /proc/sys/kernel/perf_event_paranoid)
if [[ perf_paranoid -ne 0 ]]
then
    echo "I need to set /proc/sys/kernel/perf_event_paranoid as root, please give me sudo access:"
    $(sudo sh -c "echo 0 >/proc/sys/kernel/perf_event_paranoid")
fi

echo "!!! Make sure you disabled TurboBoost via setup_scripts/turbo_boost.sh disable, I cannot check that for you !!!"

export PAPI_EVENTS="PAPI_TOT_CYC,PAPI_REF_CYC,FP_ARITH:SCALAR_DOUBLE,FP:ARITH:256B_PACKED_DOUBLE,perf::CYCLES,perf::CPU-CLOCK,perf::TASK-CLOCK"
export PAPI_OUTPUT_DIRECTORY="$PWD/logs"
# Uncomment line below to see print output
#export PAPI_REPORT="1"
log_dir="$PAPI_OUTPUT_DIRECTORY/papi_hl_output"

number_of_cores=$(grep -c ^processor /proc/cpuinfo)
in=""

if [[ "$1" = "validate" ]]
then
    read_input
    if taskset -c $(( number_of_cores-1 )) ./bin/validate debug <<< "$in" > /dev/null
    then
        rm -rf logs/papi_hl_output-*
        echo "Successfully validated."
        exit 0
    else
        rm -rf logs/papi_hl_output-*
        echo "!!! Validation error !!!"
        exit 1
    fi
else
    read_input
    if ! [[ -f "bin/$1" ]]
    then
        print_not_found "$1"
    fi
    taskset -c $(( number_of_cores-1 )) ./bin/"$1" "execute" <<< "$in"
    rm -rf logs/papi_hl_output-*
    echo "Outputs written to" logs/papi_hl_output/*
    parse_json
    exit 0
fi
