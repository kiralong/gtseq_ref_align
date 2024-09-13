#!/bin/bash
#SBATCH -p long
#SBATCH -J pyra_fastp
#SBATCH --cpus-per-task=16

#Script to run fastp to trim illumina adapter sequences and polyG tails from reads on a computing cluster with SLURM

#Load the required modules
module load fastp

# Declare variables and paths
sample_names_path=/path/to/file/with/sample/names/pyra_name_map.tsv
in_dir=/path/to/input/directory/raw_data
out_dir=/path/to/output/directory/trimmed_reads_fastp

cd $out_dir

cat $sample_names_path | cut -f 2 |
while read sample; do
	cmd=(
		fastp
		--in1 ${in_dir}/${sample}.1.fastq.gz
		--out1 ${out_dir}/${sample}.1.fq.gz
		--adapter_sequence=AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
		--trim_poly_g
		--length_required 50
		--thread 4
		--json ${out_dir}/${sample}fastp.json
		--html ${out_dir}/${sample}.fastp.html
		)
		echo ${cmd[@]}
done | parallel -j 4
