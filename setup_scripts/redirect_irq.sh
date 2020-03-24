#!/bin/bash

GRUB_PATH="/etc/default/grub"

if [[ $EUID -ne 0 ]]
then
    echo 'This script must be run with root privileges.' >&2
    exit 1
fi

number_of_cores=$(grep -c ^processor /proc/cpuinfo)

echo "$number_of_cores cores detected."

b=0
for i in $(seq 2 1 "$number_of_cores")
do
    b="${b}1"
done

b=$(printf '%x\n' $((2#$b)))

if [[ $(( ${#b} % 2 )) -ne 0 ]]
then
    b="0$b"
fi

written=0
fail=0

for i in /proc/irq/*
do
  if ! [[ -d "$i" ]]
  then
      continue
  fi

  if echo "$b" > "$i/smp_affinity" 2>/dev/null
  then
      ((written++))
  else
      ((fail++))
  fi
done

echo "Written $written IRQs, $fail not writable."
echo "Success."
