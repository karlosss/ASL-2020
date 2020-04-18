#!/bin/bash

PROJECT_ROOT=$(pwd)

function set_perf_paranoid() {
    perf_paranoid=$(cat /proc/sys/kernel/perf_event_paranoid)
    if [[ perf_paranoid -ne 0 ]]
    then
        echo "I need to set /proc/sys/kernel/perf_event_paranoid as root, please give me sudo access:"
        $(sudo sh -c "echo 0 >/proc/sys/kernel/perf_event_paranoid")
    fi
}

# install PAPI
if ! [[ -d "papi" ]]
then
    echo "PAPI is not installed. Installing PAPI..."

    set_perf_paranoid

    git clone https://bitbucket.org/icl/papi.git /tmp/papi
    mkdir papi
    cd /tmp/papi/src || exit 3
    ./configure --prefix="$PROJECT_ROOT/papi"
    make && make install
    rm -rf /tmp/papi

    echo "PAPI sucessfully installed."
    cd "$PROJECT_ROOT" || exit 3
fi

# compile
set_perf_paranoid
cmake .
make

echo "======================================================================"

export PAPI_OUTPUT_DIRECTORY="/tmp/trash"

executables=$(grep ^add_executable CMakeLists.txt | cut -d" " -f1 | cut -d"(" -f2)

if [[ -z "$executables" ]]
then
    echo "No executables found!" >&2
    exit 2
fi

cd bin || exit 3

# delete old binaries which are no longer mentioned in CMakeLists.txt
for i in *
do
    if ! [[ "$i" = ".gitkeep" ]]
    then
        found=0
        for j in $executables
        do
            if [[ "$i" = "$j" ]]
            then
                found=1
            fi
        done
        if [[ "$found" -eq 0 ]]
        then
            rm "$i"
        fi
    fi
done

# test binaries
wrong=0
correct=0

N=20
M=10
T=100
num_iter=50

for i in $executables
do
    if ! ./"$i" test <<<"$N $M $T $num_iter"
    then
        echo "!!! WARNING: binary \`$i\` produces wrong output !!!" >&2
        ((wrong++))
    else
        ((correct++))
    fi
done

echo "Testing done. Correct: $correct, wrong: $wrong."

rm -rf /tmp/trash

if [[ "$wrong" -gt 0 ]]
then
    exit 1
fi
