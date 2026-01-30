#!/bin/bash
# Grid search experiment constants
input_file="UvrD.fasta"
template_dir="UvrD_PDB"
num_recycles=0
num_seeds=15

# Grid search independant vars
n_values=("30" "50" "70")
max_msa_values=("16" "24" "32")

for n in "${n_values[@]}"; do
  for msa in "${max_msa_values[@]}"; do
    echo "$input_file" >>scan_inputs.txt
    echo "$template_dir" >>scan_inputs.txt
    echo "$num_recycles" >>scan_inputs.txt
    echo "$num_seeds" >>scan_inputs.txt
    echo "$n" >>scan_inputs.txt
    echo "$msa" >>scan_inputs.txt
    python3 wrapper.py <scan_inputs.txt
    >scan_inputs # Truncate all contents of scan_inputs.txt
  done
done
