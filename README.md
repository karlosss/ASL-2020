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

Run `compile.sh` to compiler all the binaries.

Run `run.sh <binary>` to run a binary. The script reads the standard input and passes it to the program: 
4 integers are expected and the input is not verified (yet).