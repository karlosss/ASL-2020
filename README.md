# Measurement setup

Run the following scripts (you should check what they do, provided without warranty):
- `setup_scripts/isolate_core.sh` (ensures the scheduler won't schedule anything on last CPU)
- reboot your computer
- `setup_scripts/redirect_irq.sh` (ensures all possible IRQs will be handled by other CPUs)
- `setup_scripts/turbo_boost.sh` (disables/enables turbo boost, you should have it disabled)

# Compile and run

Prerequisities:
- `gcc`
- `cmake >= 3.16`
- `make`

Run `compile.sh` to compile all the binaries. After compilation, all of them are automatically validated. Old binaries are automatically pruned. 

The `run.sh` has multiple modes:
- `run.sh validate` validates the code using a default input hardcoded in the script.
- `run.sh` runs all binaries and reports their performance. The input is read from stdin.
- `run.sh BIN` runs only binary `BIN` and reports its performance. The input is read from stdin.
