#!/bin/bash
JID=0
num_c=0
seed=1
num_s=25
m_e_msa=32
m_msa=16
inputfile=./UvrD.fasta
outputdir=/home/lcirne/scratch/container-2/lcirne0mm16-container/lcirne0mm16
temp_dir=iterations/iteration5

colabfold_batch --pair-mode unpaired_paired --templates \
--msa-mode mmseqs2_uniref_env \
--custom-template-path $temp_dir \
--max-msa $m_msa:$m_e_msa \
--use-dropout \
--num-seeds $num_s \
--num-recycle $num_c \
$inputfile $outputdir
    