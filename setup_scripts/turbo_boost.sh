if [[ $EUID -ne 0 ]]
then
    echo 'This script must be run with root privileges.' >&2
    exit 1
fi

if ! command -v wrmsr &> /dev/null
then
    echo '`wrmsr` seems not to be installed.' >&2
    exit 127
fi

if [[ $# -ne 1 ]]
then
    echo "Usage: $0 PARAM" >&2
    echo "If PARAM = e, then TurboBoost will be enabled, disabled otherwise." >&2
    exit 1
fi

number_of_cores=$(grep -c ^processor /proc/cpuinfo)

echo "$number_of_cores cores detected."

modprobe msr


if [[ "$1" != e ]]
then
    wrmsr -p"$(( $number_of_cores-1 ))" 0x1a0 0x4000850089
else
    wrmsr -p"$(( $number_of_cores-1 ))" 0x1a0 0x850089
fi

if [[ "$1" != e ]]
then
    echo "TurboBoost disabled on core $(( $number_of_cores-1 ))."
else
    echo "TurboBoost enabled on core $(( $number_of_cores-1 ))."
fi
