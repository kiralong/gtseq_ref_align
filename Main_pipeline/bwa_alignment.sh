#!/bin/bash
#SBATCH -p long
#SBATCH -J align_trimmed_pyra_gtseq_data
#SBATCH --cpus-per-task=12

# Script to align your reads to a reference genome with bwa mem 
# and sort the resulting .bam files with samtools

# Load required modules
module load bwa/0.7.17
module load samtools/1.15

# Declare variables and paths
sample_names_path=/path/to/file/with/sample/names/raw_data/info/pyra3_name_map.tsv
reads_path=/path/to/input/reads/trimmed_reads_fastp
database_path=/path/to/reference/genome/database/reference_genomes/mSylFlo1.10/mSylFlo1
out_path=/path/to/output/directory/for/aligned/reads/aligned_samples

# Run command as a while loop to go through sample list 
cat $sample_names_path | cut -f 2 |
while read sample; do
	echo "bwa mem -t 2 $database_path $reads_path/${sample}.1.fq.gz | \
		samtools view -b -h | samtools sort --threads 2 -o $out_path/${sample}.bam";
done | parallel -j 4
