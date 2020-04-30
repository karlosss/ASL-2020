#!/bin/bash

print_usage () {
    echo "Usage:" >&2
    echo "$0 validate COMPILER [FLAGS...] - run validation of base implementation (specify input on stdin)" >&2
    echo "$0 BIN COMPILER [FLAGS...] - execute binary BIN in run mode and log measurements(specify input on stdin)" >&2
    echo "input is 4 integers - N, M, T, num_iters - separated by a space"
    exit 2
}

read_input () {
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

export PAPI_EVENTS="FP_ARITH:SCALAR_DOUBLE,FP_ARITH:256B_PACKED_DOUBLE,PAPI_L3_TCM,API_L3_TCM,PAPI_L3_TCA"
export PAPI_OUTPUT_DIRECTORY="$PWD/logs"
# Uncomment line below to see print output
#export PAPI_REPORT="1"
log_dir="$PAPI_OUTPUT_DIRECTORY/papi_hl_output"

number_of_cores=$(grep -c ^processor /proc/cpuinfo)
in=""

binary="$1"

if [[ "$1" = "validate" ]]
then
  run_mode="debug"
else
  run_mode="execute"
fi

shift
compiler="$1"
shift
flags="$*"
underscored_flags=${flags// /_}

binary="bin/${compiler}_$underscored_flags/$binary"

if ! [[ -f "$binary" ]]
then
    print_not_found "$binary"
fi

read_input

if [[ "$run_mode" = "debug" ]]
then
    if taskset -c $(( number_of_cores-1 )) "$binary" "$run_mode" <<< "$in" > /dev/null
    then
        echo "Successfully validated."
        exit 0
    else
        echo "!!! Validation error !!!"
        exit 1
    fi
else
    taskset -c $(( number_of_cores-1 )) "$binary" "$run_mode" <<< "$in"
    rm -rf logs/papi_hl_output-*
    echo "Outputs written to" logs/papi_hl_output/*
    parse_json
    exit 0
fi
