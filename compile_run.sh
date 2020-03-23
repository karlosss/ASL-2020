number_of_cores=$(grep -c ^processor /proc/cpuinfo)

cmake .
make
echo =============================================
taskset -c $(( number_of_cores-1 )) ./asl_project
