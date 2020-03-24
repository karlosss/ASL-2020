number_of_cores=$(grep -c ^processor /proc/cpuinfo)

cmake .
make
echo =============================================
taskset -c $(( number_of_cores-1 )) ./measure_runtime
