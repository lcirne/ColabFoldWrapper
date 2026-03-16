#!/bin/bash

#SBATCH --nodes=1                  # Number of nodes
#SBATCH --time=24:00:00             # Wall clock time
#SBATCH --ntasks-per-node 20
#SBATCH --partition gpu

python3 wrapper.py <wrapper_input.txt
