#!/bin/bash
cmake .
make

echo "======================================================================"

executables=$(grep ^add_executable CMakeLists.txt | cut -d" " -f1 | cut -d"(" -f2)

wrong=0
correct=0

for i in $executables
do
    if ! bin/"$i" test
    then
        echo "!!! WARNING: binary \`$i\` produces wrong output !!!" >&2
        ((wrong++))
    else
        ((correct++))
    fi
done

echo "Testing done. Correct: $correct, wrong: $wrong."