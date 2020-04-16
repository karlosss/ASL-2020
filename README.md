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
