#!/bin/bash
#SBATCH -p reg              # Specify parition
#SBATCH --cpus-per-task=1   # Number of threads
#SBATCH -J syfl_bwa_db      # Name of job

# Bash script to make reference genome database on a computing cluster running SLURM

# Load the required modules
module load gcc/7.2.0
module load bwa/0.7.17

# Declare your variables and paths
work=/path/to/working/directory
ref=$work/path/to/reference/genome/genomic.fna.gz
pref=$work/specify/prefix/for/files

# Run the bwa index command
bwa index -p $pref $ref
