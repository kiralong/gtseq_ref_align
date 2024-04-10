#!/usr/bin/env python3
import sys, os, gzip, argparse, datetime

# (c) Angel G. Rivera-Colon

PROG = sys.argv[0].split('/')[-1]
DESC = '''Extract a whitelist of selected catalog markers in a SUMSTATS file 
from a set of target genomic coordinates.'''

class GenomicRange:
    def __init__(self, chromosome, basepair, name, distance=100):
        self.chr = chromosome
        self.bp  = basepair
        self.id  = name
        self.min = basepair-distance
        if self.min < 0:
            self.min = 0
        self.max = basepair+distance
        self.interval = range(self.min, self.max+1)
    def __str__(self):
        return f'{self.chr}\t{self.bp}\t{self.id}\t{self.min}\t{self.max}'
    def in_interval(self, value):
        if value in self.interval:
            return True
        else:
            return False

def parse_args():
    p = argparse.ArgumentParser(prog=PROG, description=DESC)
    p.add_argument('-s', '--sumstats', required=True,
                   help='(str) Path to SUMSTATS file')
    p.add_argument('-c', '--coordinates', required=True,
                   help='(str) Path to coordinates TSV file')
    p.add_argument('-o', '--outdir', required=False,
                   default='.', help='(str) Path to output directory')
    p.add_argument('-d', '--distance', required=False, type=int,
                   default=100, help='(int) Distance in bp plus/minus target coordinate to extract a SNP [default=100]')
    args = p.parse_args()
    if not os.path.exists(args.sumstats):
        sys.exit(f'Error: `{args.sumstats}`: SUMSTATS does not exist.')
    if not os.path.exists(args.coordinates):
        sys.exit(f'Error: `{args.coordinates}`: outlier list does not exist.')
    if not os.path.exists(args.outdir):
        sys.exit(f'Error: `{args.outdir}`: output directory does not exist.')
    if args.distance < 0:
        sys.exit(f'Error: distance ({args.distance:,}) must be >= 0.')
    return args

def now():
    '''Print the current date and time.'''
    return f'{datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")}'

def parse_coordinates(coors_f, distance=100):
    assert distance>=0
    '''Parse the coordinates file and generate a set of target regions of the genome.'''
    print(f'Parsing coordinates file...\n    Defining ranges +- {distance:,} bp from the provided coordinates.')
    coordinates = dict()
    with open(coors_f) as fh:
        for i, line in enumerate(fh):
            if line.startswith('#'):
                continue
            fields = line.strip('\n').split('\t')
            assert len(fields) == 3, "Error: Coordinates file is a TSV with three columns: chrosome, name, and marker ID."
            chrom = fields[0]
            basepair = int(fields[1])
            name = fields[2]
            # Set a new genomic range object with the corresponding coordinates
            genomic_range = GenomicRange(chrom, basepair, name, distance)
            coordinates.setdefault(chrom, [])
            coordinates[chrom].append(genomic_range)
    # Tally the retained coordinates
    # This is done now in case any duplicates are present
    n_ranges = 0
    for chrom in coordinates:
        n_ranges += len(coordinates[chrom])
    print(f'    Loaded {n_ranges:,} genomic ranges from the coordinates file.\n')
    return coordinates

def parse_sumstats(sumstats_f, coordinates, outdir='.'):
    '''Parse the SUMSTATS file and keep the loci ID located within the target coordinates.'''
    assert type(coordinates) is dict
    print('Parsing SUMSTATS file...')
    # Prepare output
    outfh = open(f'{outdir}/kept_panel_snps.tsv', 'w')
    outfh.write('#LocusID\tSnpColumn\tChrom\tBasePair\tPanelMarkerID\n')
    n_records = 0
    kept_snps = 0
    with open(sumstats_f) as fh:
        for i, line in enumerate(fh):
            if line.startswith('#'):
                continue
            n_records += 1
            fields = line.strip('\n').split('\t')
            # Set the specific columns of interest
            locus   = int(fields[0])
            chrom   = fields[1]
            bp      = int(fields[2])
            column  = int(fields[3])
            markers = list()
            # Search for them in the coordinates
            chr_coordinates = coordinates.get(chrom, None)
            if chr_coordinates is None:
                # This chromosome does not have any selected coordinates
                continue
            for ranges in chr_coordinates:
                assert isinstance(ranges, GenomicRange)
                if ranges.in_interval(bp):
                    markers.append(ranges.id)
                    kept_snps += 1
            # Save the output when the SNP is within one selected interval
            if len(markers) > 0:
                mkr_str = ';'.join(markers)
                row = f'{locus}\t{column}\t{chrom}\t{bp}\t{mkr_str}\n'
                outfh.write(row)
    print(f'    Read {n_records:,} total records from SUMSTATS.')
    print(f'    Kept {kept_snps:,} total SNPs located in the selected genomic intervals.')

def main():
    # Report start
    print(f'{PROG} started on {now()}\n')
    args = parse_args()
    # Load coordinates
    coordinates = parse_coordinates(args.coordinates, args.distance)
    # Parse sumstats and generate output
    parse_sumstats(args.sumstats, coordinates, args.outdir)
    # End report
    print(f'\n{PROG} finished on {now()}')

# Run Code
if __name__ == '__main__':
    main()
