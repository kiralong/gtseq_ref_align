#!/bin/bash
#SBATCH -p long
#SBATCH -J pyra3_fastp_populations_R50_maf0.01_hwe
#SBATCH --cpus-per-task=4

# Load required modules
module load stacks/2.60

# Declare paths and variables
thr=4
popmap=/path/to/popmap/pyra_popmap.tsv
gstacks_out=/path/to/gstacks/catalog/gstacks_runs/240305_gstacks_syfl_pyra
populations_out_path=/path/to/directory/for/populations_runs

# Create an output directory for this specific populations run
populations_output=$populations_out_path/populations_pyra3_R50_maf0.01_hwe
mkdir -p $populations_output

# Set filter parameters
R=0.50    # minimum precentage of individuals required to process a locus
maf=0.01  # minor allele frequence of 10% required to process a site at a locus

cmd=(
	populations
	--in-path $gstacks_out
	--out-path $populations_output
	--popmap $popmap
	--threads $thr
	--min-samples-overall $R
	--min-maf $maf
	--hwe
	--ordered-export
	--vcf
)
"${cmd[@]}"

