GRUB_PATH="/etc/default/grub"

if [[ $EUID -ne 0 ]]
then
    echo 'This script must be run with root privileges.' >&2
    exit 1
fi

number_of_cores=$(grep -c ^processor /proc/cpuinfo)

echo "$number_of_cores cores detected."

if [[ $number_of_cores -lt 4 ]]
then
    echo "At least 4 cores are needed. Exiting." >&2
    exit 1
fi

if ! [[ -f "$GRUB_PATH" ]]
then
    echo "GRUB configuration file not found in $GRUB_PATH. Please change the path in the script." >&2
    exit 1
fi

if grep '^GRUB_CMDLINE_LINUX_DEFAULT.*isolcpus' "$GRUB_PATH" &>/dev/null
then
    echo "Already isolated. Nothing to do." >&2
    exit 1
fi

echo "Generating new GRUB config."
sed -Ei "s/^(GRUB_CMDLINE_LINUX_DEFAULT.*)\"/\1 isolcpus=$(( number_of_cores-1 ))\"/g" "$GRUB_PATH"
grub-mkconfig -o "/boot/grub/grub.cfg" &>/dev/null

echo "Done. Please reboot."
