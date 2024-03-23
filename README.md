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
#plate_ID     barcode_i7
plate_1       CTTGTA
plate_2       GGATAA
...
plate_n       TACTGT
```

#### Well barcodes

Tab-separated file describing the per-well individual (`i5`) barcodes. The `row` and `column` describe the position of the barcode on a 96-well plate.

```sh
#row   column  barcode
A      1       CCGTTT
B      1       AAGAGT
...
H      12      AACGTT
```

#### Sample plate map

Tab-separated file describing the position of samples on a 96-well plate. It follows the following format for a single plate:

```sh
plate_1  1       2       3       4       5       6       7       8       9       10      11      12
A        SAM443  SAM451  SAM459  SAM467  SAM475  SAM483  SAM491  SAM499  SAM507  SAM515  SAM523  BE0706
B        SAM444  SAM452  SAM460  SAM468  SAM476  SAM484  SAM492  SAM500  SAM508  SAM516  SAM524  BE0707
C        SAM445  SAM453  SAM461  SAM469  SAM477  SAM485  SAM493  SAM501  SAM509  SAM517  SAM525  BE0708
D        SAM446  SAM454  SAM462  SAM470  SAM478  SAM486  SAM494  SAM502  SAM510  SAM518  SAM527  BE0709
E        SAM447  SAM455  SAM463  SAM471  SAM479  SAM487  SAM495  SAM503  SAM511  SAM519  SAM528  BE0710
F        SAM448  SAM456  SAM464  SAM472  SAM480  SAM488  SAM496  SAM504  SAM512  SAM520  SAM529  BE0711
G        SAM449  SAM457  SAM465  SAM473  SAM481  SAM489  SAM497  SAM505  SAM513  SAM521  SAM530  BE0712
H        SAM450  SAM458  SAM466  SAM474  SAM482  SAM490  SAM498  SAM506  SAM514  SAM522  SAM355  BE0713
```

Multiple 96-well plates can be processed by adding additional rows. The information for each plate must start with a column describing the plate ID (`plate_1` in the example above).

```sh
plate_1  1       2       ...  12
A        SAM443  SAM451  ...  SAM706
B        SAM444  SAM452  ...  SAM707
...
H        SAM450  SAM458  ...  SAM713
plate_2  1       2       ...  12
A        SAM654  SAM655  ...  SAM662
B        SAM444  SAM452  ...  SAM707
...
H        SAM712  SAM713  ...  SAM721
plate_n  1       2       ...  12
A        SAM898  SAM899  ...  SAM911
B        SAM444  SAM452  ...  SAM707
...
H        SAM987  SAM988  ...  SAM998
```

#### Output barcodes file

The output is a tab-separated file describing the plate barcode (`i7`) and well barcode (`i5`) of each sample:

```sh
CTTGTA  CCGTTT  SAM443
CTTGTA  AACGTT  SAM451
CTTGTA  TCAGTT  SAM459
...
CAGTCG  ATTAAA  SAM284
CAGTCG  GACAAA  SAM998
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
#Chr  Pos     ID      Ref  Alt  SampleID  Genotype  TotalDP  RefDP  AltDP
chr2  295723  90_102  G    A    SAM509    1/1       10       0      10
chr2  295723  90_102  G    A    SAM511    1/1       12       0      12
chr2  295723  90_102  G    A    SAM514    0/0       13       13     0
chr2  295723  90_102  G    A    SAM516    1/1       7        0      7
...
chr2  295723  90_102  G    A    SAM545    1/1       8        0      8
chr2  295723  90_102  G    A    SAM546    0/1       10       3      7
chr2  295723  90_102  G    A    SAM548    0/0       17       16     1
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

## Authors

Kira M. Long  
Angel G. Rivera-Colon
