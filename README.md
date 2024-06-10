# gtseq_ref_align

Pipeline and utility scripts for running a [Genotyping-in-Thousands](https://onlinelibrary.wiley.com/doi/10.1111/1755-0998.12357) (GTseq) analysis with reference aligned data

## Overall Pipeline Summary

raw GTseq reads (fastq file) -> demultiplex -> trim with `fastp` -> align with `bwa` -> genotype with *Stacks* `gstacks` -> filter with *Stacks* `populations`

## Pipeline Steps

### Step 1: Demultiplex

The sequencing facility will usually give you a single large fastq file with all of your sequencing reads for all your plates and individuals. Sometimes the sequencing facility will demultiplex the samples for you if you provide them the barcode information ahead of time. In order to demultiplex your samples from the single, raw fastq file, you will need your i7 and i5 barcode information for every sample. You can use the "grep method" which uses the unix command `grep`, or you could use the python script `GTseq_BarcodeSplit_KML.py`, which is a modified script based on Nate Campbell's demultiplexing script found here: <https://github.com/GTseq/GTseek_utils/blob/Main/GTseq_BarcodeSplit_MP.py>. 

`GTseq_BarcodeSplit_KML.py` example usage:
```
GTseq_BarcodeSplit_KML.py path/barcode.csv path/raw_reads.fastq(.gz) path/output/dir

The script uses 3 positional arguments:
     1) Path to the barcodes file (Needs plate and individual barcodes) in csv format (required).
     2) Path to the fastq file with the raw reads. It can be gzipped (required).
     3) Path to the output directory (defaults to current directory)
```

To generate the needed barcodes file you can use the script `prepare_barcodes.py` (see usage below in "Utility Scripts"). Note that you can make a bash script to submit the revised `GTseq_BarcodeSplit_KML.py` script to a queue manager on a computing cluster.

Example script to queue `GTseq_BarcodeSplit_KML.py`:
```sh
#!/bin/bash
#SBATCH -p reg
#SBATCH -J pyra4_demultiplexing
#SBATCH --cpus-per-task=10

# Load required modules
module load python

# Set variables and paths
work_dir=/path/to/working/directory
barcodes=$work_dir/path/to/barcodes/file/barcodes_gtseek.csv
fastq=$work_dir/path/to/raw/fastq/GC3F-SN-8373---7594_S0_L001_R1_001.fastq.gz
out_dir=$work_dir/path/to/individual_fastqs

# Go to the script directory
cd $work_dir/scripts

# Run the demultiplexing script
./GTseq_BarcodeSplit_KML.py $barcodes $fastq $out_dir
```
NOTE: The way barcode split script was originally written to run samples in parallel means you MUST ask for 10 cpus to run this script! The above example asks for this with the `#SBATCH --cpus-per-task=10` in the `SLURM` header. This also means if you are running this script on your own laptop you must have 10 available threads for this script to work. 
TODO: rewrite the parallelization so the user can specify the desired number of threads to use (with default 1). 

You will end up with a separate fastq file for each of your samples.

### Step 2: Trim and remove adapters

Use the script `run_fastp.sh` to trim and remove adapter sequences from your raw, demultiplexed reads. To run this script you will need the software `fastp` and `parallel` along with a list of your sample names.

### Step 3: Align to a reference genome

Use the script `bwa_alignment.sh` to align your trimmed fastq files to a reference assembly. This script requires that you have `bwa`, `samtools`, and `parallel` either installed and in your path on your cluster/computer or the module loaded if it is already available on your computing cluster. Make sure you keep track of which versions of the software you are using as some features change with version updates and you will need to report the versions of the software you used in your methods. This script should export a .bam file for each of your samples.

Note that before you can align to a reference assembly you will need to download and make a database of the assembly. You can use the unix command `wget <url>` to download the reference assembly (a genomic.fna.gz file) from NCBI. I usually use the NCBI FTP directory to copy the link address of the file for `wget`. After you download the reference assembly, make sure you record the NCBI RefSeq accession number (ex. GCF_025592945.1) of the genome you are using as this must be reported in future manuscripts. Then you can run the script `bwa_db.sh` to make your reference genome database. When I specify the prefix (`$pref`) I generally make this a shorthand of the species name of the reference genome I am using. This prefix will be added to all the files generated for the database.

When editing the paths in the `bwa_alignment.sh` script, make sure the `database_path` variable ends with the same prefix you assigned the reference database in `bwa_db.sh`. Also note that in the script I currently use cut to get the sample names from the second column of my sample list. This is because I was reusing a list with sample names in the second column. If your sample list is only one column, you don't need to use `cut -f 2`. Similarly, if you are using a csv or tsv file with sample names in the 3rd, 4th etc columns use `cut -f <integer>` to specify which column you want. 

### Step 4: Genotype

Use the script `run_gstacks.sh` to genotype your aligned reads. You will need *Stacks* version 2.0 or later. You will also need a `popmap` file to run `gstacks`. The `popmap` must be a tsv file, and is 2 columns separated by a tab, the first column is the sample ID and the second column is the population the individual belongs to. If you do not want to separate the samples by populations, just put a placeholder in the second column like `1` to note they are all in one population.

After you get your `gstacks` catalog and log files, you should export and check some baseline sample quality information using the `stacks-dist-extract` utily available in *Stacks*. I recommend always checking how well your samples aligned to the reference genome (`bam_stats_per_sample`) and your coverages (`effective_coverages_per_sample`) from your `gstacks.log.distribs` file. Redirect the output with `> file_name.tsv` to keep this information for future reference, reporting, and/or graphing.

Example usage and output:
```sh
stacks-dist-extract ./stacks/gstacks.log.distribs bam_stats_per_sample
sample	records	primary_kept	kept_frac	primary_kept_read2	primary_disc_mapq	primary_disc_sclip	unmapped	secondary	supplementary
S1_2023.01	2780637	2515438	0.905	1195103	26801	98337	80108	0	59953
S1_2023.07	3156646	2860191	0.906	1359700	27987	110763	89513	0	68192
S2_1999.13	2835542	2574684	0.908	1225169	25379	96962	81343	0	57174
...
```

### Step 5: Filter

Use the script `run_populations.sh` to filter your genotyped data and export desired file formats such as vcf and plink files. I recommend when you are first going through your data to run a "base filtering" run of your data to see the number of loci/SNPs kept and to check the missing data per individual. Rerun `populations` with stricter filtering schemes as needed. 

Remember that your filtering scheme (and hence what counts as a 'base' filtering) will change depending on the number of samples you are running. You will want to run a minor allele count (mac) or minor allele frequency (maf) filter but maf especially, being a frequency, will change drastically depending on the number of samples you are running. As a starting point, I often use -mac 3 or -maf 0.01 as a base filter.

After running `populations` and getting a `populations.log.distribs` file, I recommend checking the missing data per individual (variant_sites_per_sample). You can extract this information using `stacks-dist-extract` just like for the `gstacks.log.distribs' file. 

Example usage and output:
```sh
stacks-dist-extract ./populations.log.distribs variant_sites_per_sample
# Number of variant sites per individual sample (after filtering).
sample	n_sites	present_sites	missing_sites	frequency_missing
S1_2023.01	1256	1042	214	0.1704
S1_2023.07	1256	1205	51	0.0406
S2_1999.13	1256	1178	78	0.062
...
```

Now that you have some quality data for your full run from both `gstacks` and `populations` I would go back and remove any samples with low coverage and high missing data from your popmap, rerun `gstacks` to get a new gstacks catalog, and then rerun `populations` with your final filtering scheme. Removing the "bad apple" samples is highly recommended for RADseq datasets, and I believe will help with GTseq datasets as well. See: <https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/2041-210X.13562>

## Utility scripts

### Create barcode file

Generate a barcode file containing the per-sample plate and well barcode sequences. Requires both plate and well barcode files, as well as a sample plate map.

#### Plate barcodes

Tab-separated file describing the per-plate (`i7`) barcodes.

```sh
#plate_ID<tab>i7_name<tab>barcode_i7
plate_1       GTseq_01    CTTGTA
plate_2       GTseq_02    GGATAA
...
plate_n       GTseq_08    TACTGT
```

#### Well barcodes

Tab-separated file describing the per-well individual (`i5`) barcodes. The `row` and `column` describe the position of the barcode on a 96-well plate.

```sh
#row<tab>column<tab>barcode
A        1          CCGTTT
B        1          AAGAGT
...
H        12         AACGTT
```

#### Sample plate map

Tab-separated file describing the position of samples on a 96-well plate. It follows the following format for a single plate:

```sh
plate_1<tab>1 <tab> 2 <tab> 3 <tab> 4 <tab> 5 <tab> 6 <tab> 7 <tab> 8 <tab> 9 <tab> 10 <tab>11 <tab>12
A           SAM443  SAM451  SAM459  SAM467  SAM475  SAM483  SAM491  SAM499  SAM507  SAM515  SAM523  BE0706
B           SAM444  SAM452  SAM460  SAM468  SAM476  SAM484  SAM492  SAM500  SAM508  SAM516  SAM524  BE0707
C           SAM445  SAM453  SAM461  SAM469  SAM477  SAM485  SAM493  SAM501  SAM509  SAM517  SAM525  BE0708
D           SAM446  SAM454  SAM462  SAM470  SAM478  SAM486  SAM494  SAM502  SAM510  SAM518  SAM527  BE0709
E           SAM447  SAM455  SAM463  SAM471  SAM479  SAM487  SAM495  SAM503  SAM511  SAM519  SAM528  BE0710
F           SAM448  SAM456  SAM464  SAM472  SAM480  SAM488  SAM496  SAM504  SAM512  SAM520  SAM529  BE0711
G           SAM449  SAM457  SAM465  SAM473  SAM481  SAM489  SAM497  SAM505  SAM513  SAM521  SAM530  BE0712
H           SAM450  SAM458  SAM466  SAM474  SAM482  SAM490  SAM498  SAM506  SAM514  SAM522  SAM355  BE0713
```

Multiple 96-well plates can be processed by adding additional rows. The information for each plate must start with a column describing the plate ID (`plate_1` in the example above).

```sh
plate_1<tab>1 <tab> 2 <tab> ... 12
A           SAM443  SAM451  ... SAM706
B           SAM444  SAM452  ... SAM707
...
H           SAM450  SAM458  ... SAM713
plate_2     1       2       ... 12
A           SAM654  SAM655  ... SAM662
B           SAM444  SAM452  ... SAM707
...
H           SAM712  SAM713  ... SAM721
plate_n     1       2       ... 12
A           SAM898  SAM899  ... SAM911
B           SAM444  SAM452  ... SAM707
...
H           SAM987  SAM988  ... SAM998
```

#### Output barcodes file

The script produces two different barcode output files.

The first one, `barcodes_stacks.tsv`, is a tab-separated file describing the plate barcode (`i7`) and well barcode (`i5`) of each sample in the *Stacks* format (section 4.1.2 of thr *Stacks* [manual](https://catchenlab.life.illinois.edu/stacks/manual/#specbc)).

```sh
CTTGTA<tab>CCGTTT<tab>SAM443
CTTGTA     AACGTT     SAM451
CTTGTA     TCAGTT     SAM459
...
CAGTCG     ATTAAA     SAM284
CAGTCG     GACAAA     SAM998
```

The second output, `barcodes_gtseek.csv`, is a comma-separated file following the [GTseek format](https://github.com/GTseq/GTseek_utils/blob/Main/GTseq_BarcodeSplit_MP.py). It specifies the sample info, as well as the plate barcode (`i7_name` and `i7_sequence`) and well barcode (`i5_name` and `i5_sequence`).

```sh
Sample,PlateID,i7_name,i7_sequence,i5_name,i5_sequence
SAM403,plate_1,GTseq_01,CTTGTA,A1,CCGTTT
SAM521,plate_1,GTseq_01,CTTGTA,A2,AACGTT
SAM539,plate_1,GTseq_01,CTTGTA,A3,TCAGTT
...
SAM284,plate_8,GTseq_08,CAGTCG,H11,ATTAAA
SAM998,plate_8,GTseq_08,CAGTCG,H12,GACAAA

```

#### Usage

```sh
prepare_barcodes.py started on 2024-03-22 17:04:42

usage: prepare_barcodes.py [-h] -p PLATE_BARCODES -w WELL_BARCODES -m SAMPLE_MAP [-o OUTDIR]

Prepare the set of per-individual barcodes for GTseq analysis.

options:
  -h, --help            show this help message and exit
  -p PLATE_BARCODES, --plate-barcodes PLATE_BARCODES
                        (str) Path to the per-plate barcodes file.
  -w WELL_BARCODES, --well-barcodes WELL_BARCODES
                        (str) Path to the per-well barcodes file.
  -m SAMPLE_MAP, --sample-map SAMPLE_MAP
                        (str) Path to the multi-plate sample map file.
  -o OUTDIR, --outdir OUTDIR
                        (str) Path to output directory.
```

#### Example

```sh
$ python3 prepare_barcodes.py \
    --plate-barcodes plate_barcodes.tsv \   # Path to plate barcodes file
    --well-barcodes well_barcodes.tsv \     # Path to well barcodes file
    --sample-map sample_map.tsv \           # Path to sample plate map
    --outdir output/                        # Path to output directory
```

### Extract probe reads

Extract the reads matching a given probe sequence from the GTseq reads of a given individual and export them in FASTQ format. They can be used to verify the target specificity of a given probe.

#### Probe sequences

The target probe sequences are used to define each genotyped marker during GTseq analysis. Is marker is designed to match a given probe. For example, these probe sequences are specified to the GTseq pipeline in the comma-separated probe file:

```sh
locus_1,A,G,TCCCAAAATGC,TCCCAGAATGC,AACATCAGAGACGCCAGGTG,0,0
locus_2,G,A,CCAGTGTTCAA,CCAGTATTCAA,CTTGGCATTCGTCAAATCCAGT,0,0
locus_3,G,A,GTGACGGATGG,GTGACAGATGG,GAGGAGCGGCAGGTGTG,0,0
locus_4,C,T,TGGCCCGATAT,TGGCCTGATAT,TGCAGGCAGTGGCTTAAAGA,0,0
```

The probes can be found of the 6th field of this file, e.g., `AACATCAGAGACGCCAGGTG` in the first line.

#### Usage

```sh
extract_probe_reads.py started on 2024-03-22 17:06:36

usage: extract_probe_reads.py [-h] -f FASTQ -t TARGET_SEQ [-o OUTDIR] [-n N_MISMATCH]

Extract the records that match a target probe/barcode sequence in a given FASTQ file.

options:
  -h, --help            show this help message and exit
  -f FASTQ, --fastq FASTQ
                        (str) Path to target FASTQ file.
  -t TARGET_SEQ, --target-seq TARGET_SEQ
                        (str) Target barcode/probe sequence to search in reads.
  -o OUTDIR, --outdir OUTDIR
                        (str) Path to output directory.
  -n N_MISMATCH, --n-mismatch N_MISMATCH
                        (int) Number of allowed sequence mismatched when searching target sequence. [default=1]
```

#### Example

```sh
$ ./extract_probe_reads.py \
    --fastq SAM123.GTseq.fastq.gz \        # Path to sample FASTQ
    --target-seq AACATCAGAGACGCCAGGTG \    # Target probe sequence
    --outdir output/                       # Path to output directory
    --n-mismatch 1                         # Tolerate 1 mistmatch between target sequence and read
```

### Extract the allelic ratios of a target site

Extract the per-sample allelic ratios at a given genotyped site. These are analogous allelic ratio diagnostic plots generated by the GTseek pipeline.

These ratios are extracted from the genotype allelic depth (`AD`) field from a VCF generated from the Stacks' `populations` software. The target site is defined according to the Stacks catalog id format, `<locus_id>_<column>`. 

#### Output table

```sh
#Chr<tab>Pos<tab>ID <tab> Ref<tab>Alt<tab>SampleID<tab>Genotype<tab>TotalDP<tab>RefDP<tab>AltDP
chr2     295723  90_102   G       A       SAM509       1/1          10          0         10
chr2     295723  90_102   G       A       SAM511       1/1          12          0         12
chr2     295723  90_102   G       A       SAM514       0/0          13          13        0
chr2     295723  90_102   G       A       SAM516       1/1          7           0         7
chr2     295723  90_102   G       A       SAM545       1/1          8           0         8
chr2     295723  90_102   G       A       SAM546       0/1          10          3         7
chr2     295723  90_102   G       A       SAM548       0/0          17          16        1
```

The output table described the genomic coordinates of the site (`Chr` and `Pos`), when available, as well as the Stacks' catalog id (`ID`). For each individual, the table reports the individual's genotype at that site (using the VCF format's notation), the total depth (`TotalDP`), and the depth of the two reference (`RefDP`) and alternative (`AltDP`) allele.

#### Usage

```sh
usage: extract_vcf_stats.py [-h] -v VCF [-o OUTDIR] [-t TARGET_SITE]

Extract genotype, coverage, and allele depth data for a target site in a VCF file

options:
  -h, --help            show this help message and exit
  -v VCF, --vcf VCF     (str) Path to the genotype data in VCF format.
  -o OUTDIR, --outdir OUTDIR
                        (str) Output directory.
  -t TARGET_SITE, --target-site TARGET_SITE
                        (str) ID of target site to extract.
````

#### Example

```sh
$ python3 extract_vcf_stats.py \
    --vcf input.vcf \               # Path to input VCF
    --target-site 90_102 \          # ID of target SNP, maching the Stacks catalog ID
    --outdir output/                # Path to output directory
```

### Match *Stacks* catalog and marker set

Filter the *Stacks* SUMSTATS file (`populations.sumstats.tsv`) to obtain the ID of all the SNPs that are within a set of specified genomic intervals. These intervals are often the known genomic coordinates of the original GTseq panel marker set.

The SUMSTATS file follows the standard format specified by *Stacks* (section 6.6.1 of the *Stacks* [manual](https://catchenlab.life.illinois.edu/stacks/manual/#pfiles)) and can be generated after generating a reference-based catalog using `gstacks` and applying some baseline filters with `populations`.

#### Genomic coordinates

The genomic coordinates are specified in the form of a BED file specifying the genomic coordinates (chromosome, start bp, and end bp) and the ID of the markers. These can be obtained, e.g., by aligning the sequence of the panel markers to the reference, processing the alignments, and exporting coordinates in BED format (`#TODO: add example of this`). By default, the script keeps any SNP within the interval specified in the BED (can be extended with available options). Each SNP within this interval is assiged to the corresponding marker.

```sh
#chrom<tab>startBP<tab>endBP<tab>marker_id>
chr1       61724575    61725396  panel_marker_01
chr1       61558608    61558960  panel_marker_02
chr2       50458016    50458730  panel_marker_03
chr3       48307895    48308637  panel_marker_04
chr5       47639353    47639563  panel_marker_05
chr8       21953369    21954103  panel_marker_06
chr8       21954099    21954842  panel_marker_07
chr10      12328678    12329324  panel_marker_08
chr11      12057388    12058136  panel_marker_09
chr12      6025344     6025973   panel_marker_10
```

Note: The BED file might have other columns, but only the first four are used and must follow this convention.

#### Output table

The output of the script (`kept_panel_snps.tsv`) describes the ID and location of the retained SNPs.

```sh
#LocusID<tab>SnpColumn<tab>Chrom<tab>BasePair<tab>PanelMarkerID
10           28            chr1      61724582      panel_marker_01
19           18            chr1      61558631      panel_marker_02
22           34            chr2      50458037      panel_marker_03
24           45            chr4      48307901      panel_marker_04
27           63            chr5      47639373      panel_marker_05
32           49            chr6      21954108      panel_marker_06
33           57            chr6      21954122      panel_marker_07;panel_marker_08
41           20            chr11     12057400      panel_marker_09
41           74            chr11     12057454      panel_marker_09
```

The table species the genomic coordinates of the SNP (`Chrom` and `BasePairs`) and well as the marker ID in the *Stacks* catalog (`LocusID` and `SnpColumn`). These coordinates are matched against the panel markers (`PanelMarkerID`).

In the example above, SNP `10_28` (first row) in the *Stacks* catalog corresponds to panel marker 1. It is possible for one SNP to match than one marker (e.g., if the loci truly overlap in the genome and/or the the user specified distance results in the spans of multiple marker intervals). We see this for SNP `33_57`, which is within the span of panel markers 7 and 8.

This output table can be used to generate a whitelist in the *Stacks* format (section 4.4.4 of the [manual](https://catchenlab.life.illinois.edu/stacks/manual/#wl)), by selecting just the first two columns describing the loci ID and SNP columns. This whitelist can be used to, e.g., run an analysis in `populations` including *just* the set of SNPs matching the GTseq panel.

```sh
$ cat kept_panel_snps.tsv | grep -v '#' | cut -f 1,2 > gtseq_panel.whitelist.tsv
$ head gtseq_panel.whitelist.tsv
10<tab>28
19     18
22     34
24     45
27     63
32     49
33     57
41     20
41     74
```

#### Usage

```sh
usage: extract_loci_in_coords.py [-h] -s SUMSTATS -c COORS
       [-o OUTDIR] [-d DISTANCE]

Extract a whitelist of selected catalog markers in a SUMSTATS 
file from a set of target genomic coordinates.

options:
  -h, --help            show this help message and exit
  -s SUMSTATS, --sumstats SUMSTATS
                        (str) Path to SUMSTATS file
  -c COORS, --coors COORS
                        (str) Path to coordinates BED file
  -o OUTDIR, --outdir OUTDIR
                        (str) Path to output directory
  -d DISTANCE, --distance DISTANCE
                        (int) Distance in bp plus/minus target 
                        coordinate to extract a SNP [default=0]
```

#### Example

```sh
$ python3 extract_loci_in_coords.py \
    --sumstats /path/to/populations.sumstats.tsv \    # Stacks SUMSTATS
    --coors-bed /path/to/coordinates.bed \            # Coordinates BED
    --outdir output/                                  # Outdir
```

## Authors

Kira M. Long  
Angel G. Rivera-Colon
