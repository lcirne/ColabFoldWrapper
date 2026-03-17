#!/bin/bash
# Grid search experiment constants
input_file="UvrD.fasta"
template_dir="UvrD_PDB"
num_recycles=0
num_seeds=15

# Grid search independant vars
n_values=("30" "50" "70")
max_msa_values=("16" "24" "32")

n_len=${#n_values[@]}
msa_len=${#max_msa_values[@]}
num_containers=$((n_len * msa_len))

# Containers
containers=()
for ((i = 0; i < num_containers; i++)); do
  containers+=("container-$i")
done

i=0
for n in "${n_values[@]}"; do
  for msa in "${max_msa_values[@]}"; do
    echo "$input_file" >>scan_inputs.txt
    echo "$template_dir" >>scan_inputs.txt
    echo "$num_recycles" >>scan_inputs.txt
    echo "$num_seeds" >>scan_inputs.txt
    echo "$n" >>scan_inputs.txt
    echo "$msa" >>scan_inputs.txt

    mkdir -p "${containers[$i]}"
    cp scan_inputs.txt "${containers[$i]}/"
    cp "$input_file" "${containers[$i]}/"
    cp wrapper.py "${containers[$i]}/"
    cp data_engine.py "${containers[$i]}/"
    cp -r "$template_dir" "${containers[$i]}/"
    cp -r distance_finder/ "${containers[$i]}/"
    cp queue_wrapper.sh "${containers[$i]}/"
    cd "${containers[$i]}/"
    sbatch --job-name="wrapper-${containers[$i]}" queue_wrapper.sh
    cd ..
    >scan_inputs.txt # Truncate all contents of scan_inputs.txt
    ((i++))
  done
done
