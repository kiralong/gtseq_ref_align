#!/bin/bash
#SBATCH -p reg
#SBATCH -J pyra4_process_radtags
#SBATCH --cpus-per-task=1
#SBATCH --mem=64G

# Load required modules
module load stacks/2.60

# Define paths and variables
in_fastq=/path/to/RAW_S0_L001_R1_001.fastq.gz
barcodes=/path/to/barcodes.tsv
out_dir=/path/to/processed_reads

# stacks process_radtags command (@KML not doing clean and quality flags because doing cleaning in fastp)
cmd=(
    process_radtags
    -f $in_fastq #In this case manually specifying fastq file, this might change
    -b $barcodes
    -o $out_dir
    --rescue
    --index-index #for i7 and i5 on fastq header MUST specify this
    --retain-header #This will retain the barcodes in the header for checking later
    --disable-rad-check #since this is not RAD data, specify this
)
"${cmd[@]}"
