# Measurement setup

Run the following scripts (you should check what they do, provided without warranty):
- `setup_scripts/isolate_core.sh` (ensures the scheduler won't schedule anything on last CPU)
- Reboot your computer
- `setup_scripts/redirect_irq.sh` (ensures all possible IRQs will be handled by other CPUs)
- `setup_scripts/turbo_boost.sh` (disables/enables turbo boost, you should have it disabled)

# Register implementations
- Create a new `*.cpp` file in `opts` directory
- Implement functions `size_t flop_count(int N, int M, int T, int n_iter)` and `void baum_welch(double* PI, double* A, double* B, int* O, double* FW, double* BW, double* C, int N, int M, int T, int n_iter)`
- Add a record to `CMakeLists.txt`

# Compile and run

Prerequisities:
- `gcc`
- `cmake >= 3.10`
- `make`
- `git`
- `sudo`

## Compile and run
Run `compile.sh` to compile all the binaries. After compilation, all of them are automatically tested. Old binaries are automatically pruned. 

The `run.sh` has multiple modes:
- `run.sh validate` - run validation of base implementation (specify input on stdin)
- `run.sh BIN` - execute binary BIN in run mode and log measurements(specify input on stdin)
- Input is 4 integers - N, M, T, num_iters - separated by a space

## Running experiments
An 'experiment' consists of:
- A source file containing optimization code to be run.
- A compiler to to compile the source file.
- A list of compiler flags for the specified comiler.
- Three ranges for the three performance plots to be generated:
  * N plot: Varies N while keeping T and M fixed.
  * M plot: Varies M while keeping N and T fixed.
  * T plot: Varies T while keeping N and M fixed.
- A flag 'e' specified whether the ranges are sampled linearly or exponentially (base 2)
For more details run `./experiments.py --help`.

The experiment involves compiling the source code with the specified compiler and flags. The resulting binary is then run over the specified ranges while benchmarking the performances for the different inputs and storing the measured data in a .csv file in 'output_data'. Furthermore, a performance plot is generated automatically after a run of `experiments.py` and stored under `figures/`.

Example:
`./experiments.py -s no_opt -c 'g++' -f 'Ofast' -N='1,4,3' -M='1,6,5' -T='6,10' -p -e `

## Comparison Plots
Run `plot.py` to generate plot comparing the performance measurements already measured from previous runs of `experiments.py`. The options are:
- `-f FILELIST`: Compare the performances of all .csv files in the comma separated list `FILELIST`.
- `-d DIRECTORY`: Compare the performances of all .csv files inside `DIRECTORY`.

If no arguments are given the folder `./output_dir` is used as the directory.