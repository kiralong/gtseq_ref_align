#!/bin/bash
#SBATCH -p long
#SBATCH -J pyra_gstacks
#SBATCH --cpus-per-task=8

# Script for genotyping using gstacks

# Load required modules
module load stacks/2.60

# Declare paths and variables
work_dir=/path/to/working/directory/gstacks_runs
out=$(date +${work_dir}/%y%m%d_gstacks_syfl_pyra) #to create an output directory with the date
mkdir -p $out
popmap=/path/to/popmap/pyra_popmap.tsv
aligned_samples=/path/to/aligned_samples
threads=8

gstacks \
        -I $aligned_samples \
        -O $out \
        -M $popmap \
        -t $threads \
        --details
