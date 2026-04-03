#!/bin/bash

#SBATCH --nodes=1                  # Number of nodes
#SBATCH --time=24:00:00             # Wall clock time
#SBATCH --ntasks-per-node 20
#SBATCH --partition gpu

python3 wrapper.py <scan_inputs.txt

find . -mindepth 1 ! -name '*mm*' ! -name 'slurm*' ! -name 'run_job.sh' -exec rm -rf {} +
