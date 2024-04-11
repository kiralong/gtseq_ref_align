# gtseq_ref_align
Pipeline and utility scripts for running a [Genotyping-in-Thousands](https://onlinelibrary.wiley.com/doi/10.1111/1755-0998.12357) (GTseq) analysis with reference aligned data

## Overall Pipeline Summary
raw GTseq reads (fastq file) -> demultiplex -> trim with fastp -> align with bwa -> genotype with Stacks

## Pipeline Steps

### Step 1: Demultiplex

### Step 2: Trim and remove adapters

### Step 3: Align to a reference genome

### Step 4: Genotype

## Utility scripts

### Create barcode file

Generate a barcode file containing the per-sample plate and well barcode sequences. Requires both plate and well barcode files, as well as a sample plate map.

#### Plate barcodes

Tab-separated file describing the per-plate (`i7`) barcodes.

```sh
#plate_ID<tab>barcode_i7
plate_1       CTTGTA
plate_2       GGATAA
...
plate_n       TACTGT
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

The output is a tab-separated file describing the plate barcode (`i7`) and well barcode (`i5`) of each sample:

```sh
CTTGTA<tab>CCGTTT<tab>SAM443
CTTGTA     AACGTT     SAM451
CTTGTA     TCAGTT     SAM459
...
CAGTCG     ATTAAA     SAM284
CAGTCG     GACAAA     SAM998
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

### Match *Stacks* catalog and maker set

Filter the *Stacks* SUMSTATS file (`populations.sumstats.tsv`) to obtain the ID of all the SNPs that are within a set of specified genomic intervals. These intervals are often the known genomic coordinates of the original GTseq panel marker set.

The SUMSTATS file follows the standard format specified by *Stacks* (section 6.6.1 of the *Stacks* [manual](https://catchenlab.life.illinois.edu/stacks/manual/#pfiles)) can an be generated after generating a reference-based catalog using `gstacks` and applying some baseline filters with `populations`.

#### Genomic coordinates

The genomic coordinates are specified in the form of a TSV file specifying the genomic coordinates (chromosome and basepair) and the ID of the markers. These can be obtained, e.g., by aligning the sequence of the panel makers to the references and obtaining the coordinates based on the SAM/BAM alignment (`#TODO: add example of this`). By default, the script keeps any SNP within a specicied distance (+-100bp by default, but it should match the average read length) from the specified basepair. Each SNP within this interval is assiged to the corresponding marker.

```sh
#chrom<tab>basepair<tab>marker_id
chr1       107998114    panel_marker_01
chr1       117784811    panel_marker_02
chr2       11971911     panel_marker_03
chr3       14963079     panel_marker_04
chr5       18555876     panel_marker_05
chr8       2063927      panel_marker_06
chr9       21533701     panel_marker_07
chr10      30672183     panel_marker_08
chr11      40580125     panel_marker_09
chr12      42676185     panel_marker_10
```

#### Output table

The output of the script (`kept_panel_snps.tsv`) describes the ID and location of the retained SNPs.

```sh
#LocusID<tab>SnpColumn<tab>Chrom<tab>BasePair<tab>PanelMarkerID
10           28            chr1      2064007      panel_marker_01
21           18            chr2      6128794      panel_marker_02
21           34            chr2      6128810      panel_marker_03
21           45            chr2      6128821      panel_marker_04
21           63            chr2      6128839      panel_marker_05
28           49            chr7      8833311      panel_marker_06
28           57            chr7      8833319      panel_marker_07;panel_marker_08
41           20            chr8      11971977     panel_marker_09
41           74            chr8      11972031     panel_marker_09
```

The table species the genomic coordinates of the SNP (`Chrom` and `BasePairs`) and well as the maker ID in the *Stacks* catalog (`LocusID` and `SnpColumn`). These coordinates are matched against the panel makers (`PanelMarkerID`).

In the example above, SNP `10_28` (first row) in the *Stacks* catalog corresponds to panel marker 1. It is possible for one SNP to match than one marker (e.g., if the loci truly overlap in the genome and/or the the user specified distance results in the spans of multiple marker intervals). We see this for SNP `28_57`, which is within the span of panel markers 7 and 8.

#### Usage

```sh
usage: extract_loci_in_coords.py [-h] -s SUMSTATS -c COORDINATES [-o OUTDIR] [-d DISTANCE]

Extract a whitelist of selected catalog markers in a SUMSTATS file from a set of target genomic coordinates.

options:
  -h, --help            show this help message and exit
  -s SUMSTATS, --sumstats SUMSTATS
                        (str) Path to SUMSTATS file
  -c COORDINATES, --coordinates COORDINATES
                        (str) Path to coordinates TSV file
  -o OUTDIR, --outdir OUTDIR
                        (str) Path to output directory
  -d DISTANCE, --distance DISTANCE
                        (int) Distance in bp plus/minus target coordinate to extract a SNP [default=100]
```

#### Example

```sh
$ python3 extract_loci_in_coords.py \
    --sumstats /path/to/populations.sumstats.tsv \    # Stacks SUMSTATS
    --coordinates /path/to/coordinates.tsv \          # Coordinates
    --outdir output/ \                                # Outdir
    --distance 150                                    # Search 150 bp away
```

## Authors

Kira M. Long  
Angel G. Rivera-Colon
