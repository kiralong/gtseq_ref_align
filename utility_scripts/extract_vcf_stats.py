#!/usr/bin/env python3
import sys, os, argparse, gzip, datetime

# (c) 2024 Angel G. Rivera-Colon

PROG = sys.argv[0].split('/')[-1]
DESC = '''Extract genotype, coverage, and allele depth data for a target
site in a VCF file'''

def parse_args():
    p = argparse.ArgumentParser(prog=PROG, description=DESC)
    p.add_argument('-v', '--vcf', required=True,
                   help='(str) Path to the genotype data in VCF format.')
    p.add_argument('-o', '--outdir', required=False, default='.',
                   help='(str) Output directory.')
    p.add_argument('-t', '--target-site', required=False, default=None,
                   help='(str) ID of target site to extract.')
    # Check input arguments
    args = p.parse_args()
    args.outdir = args.outdir.rstrip('/')
    if not os.path.exists(args.outdir):
        sys.exit(f"Error: '{args.outdir}' not found.")
    if not os.path.exists(args.vcf):
        sys.exit(f"Error: '{args.vcf}' not found.")
    return args

#
# Classes
#

class SnpInfo():
    def __init__(self, chrom, bp, snpid, ref, alt):
        self.chr = chrom
        self.bp  = bp
        self.id  = snpid
        self.ref = ref
        self.alt = alt
    def __str__(self):
        return f'{self.chr} {self.bp} {self.id} {self.ref} {self.alt}'
    
class IndGenotype():
    def __init__(self, sample_id, site_id, genotype, total_depth=0, ref_depth=0, alt_depth=0):
        assert type(total_depth) is int
        assert type(ref_depth) is int
        assert type(ref_depth) is int
        assert (ref_depth+alt_depth) <= total_depth
        # Sample and site ID
        self.sam = sample_id
        self.id  = site_id
        # Process genotypes and depth info
        self.geno  = genotype
        self.totdp = total_depth
        self.altdp = alt_depth
        self.refdp = ref_depth
        # Process missing genotypes
        if genotype in ('./.', '.|.'):
            self.geno  = 'NA'
            self.totdp = 'NA'
            self.altdp = 'NA'
            self.refdp = 'NA'
    def __str__(self):
        return f'{self.sam} {self.id} {self.geno} {self.totdp} {self.refdp} {self.altdp}'

#
# Functions
#

def now():
    '''Print the current date and time.'''
    return f'{datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")}'

def time():
    '''Print today's date'''
    return f'{datetime.datetime.now().strftime("%Y%m%d")}'

def initialize_log(args, log=sys.stdout):
    '''Initialize the log file'''
    print(f'''{PROG} started on {now()}\n
    Extracting data for locus {args.target_site}''', file=log, flush=True)

def parse_vcf(vcf_f, target_site=None, log=sys.stdout):
    '''Parse the VCF to obtain the genotype info.'''
    print('\nParsing VCF...', file=log, flush=True)
    vcf_rec  = 0
    snp_info = dict()
    geno_info = dict()
    vcf_samples = None
    # Parse VCF
    with gzip.open(vcf_f, 'rt') if vcf_f.endswith('.gz') else open(vcf_f) as fh:
        for line in fh:
            # The header lines need to be added to the new output VCF file
            if line.startswith('##'):
                continue
            elif line.startswith('#CHROM'):
                vcf_samples = line.strip('\n').split('\t')[9:]
                print(f'    Processing genotypes for {len(vcf_samples):,} samples.', file=log, flush=True)
                continue
            else:
                # Process all other rows which contain variant sites
                vcf_rec += 1
                # Extract the reference info for each variant site,
                # which is used to reconstruct the final VCF
                fields = line.strip('\n').split('\t')
                chrom = fields[0]
                basepair = fields[1]
                id = fields[2]
                ref = fields[3]
                alt = fields[4]
                # Process the ID to match the Default Stacks `locus_col`
                id = '_'.join(id.split(':')[:2])
                # Filter non-target sites
                if target_site is not None:
                    if not id == target_site:
                        continue
                info = SnpInfo(chrom, basepair, id, ref, alt)
                snp_info[id] = info
                # Process the genotype columns
                geno_info.setdefault(id, [])
                for i, geno in enumerate(fields[9:]):
                    sample = vcf_samples[i]
                    site_geno = None
                    geno_field = geno.split(':')
                    # Process missing sites
                    if len(geno_field) == 1:
                        site_geno = IndGenotype(sample, id, geno_field[0])
                    else:
                        geno = geno_field[0]
                        total_dp = int(geno_field[1])
                        ref_dp = int(geno_field[2].split(',')[0])
                        alt_dp = int(geno_field[2].split(',')[1])
                        site_geno = IndGenotype(sample, id, geno, total_dp, ref_dp, alt_dp)
                    geno_info[id].append(site_geno)
        # Report and Return
        print(f'''    Read {vcf_rec:,} records from input VCF.
    Retained data for {len(snp_info):,} variant sites.''', file=log, flush=True)
        return snp_info, geno_info

def write_output_tbl(snp_info, geno_info, outdir='.'):
    '''Export the data for target site as a tsv.'''
    header = '#Chr\tPos\tID\tRef\tAlt\tSampleID\tGenotype\tTotalDP\tRefDP\tAltDP\n'
    for site in snp_info:
        snp = snp_info[site]
        # Open the output for a given SNP
        fh = open(f'{outdir}/vcf_stats.{snp.id}.tsv', 'w')
        fh.write(header)
        # Process the individual genotypes
        snp_genos = geno_info[site]
        for ind in snp_genos:
            row = f'{snp.chr}\t{snp.bp}\t{snp.id}\t{snp.ref}\t{snp.alt}\t{ind.sam}\t{ind.geno}\t{ind.totdp}\t{ind.refdp}\t{ind.altdp}\n'
            fh.write(row)
        fh.close()

def main():
    args = parse_args()
    # Initialize log file
    log_f = open(f'{args.outdir}/gtseq_vcf_stats.log', 'w')
    initialize_log(args, log_f)

    # Parse VCF
    snp_info, geno_info = parse_vcf(args.vcf, target_site=args.target_site, log=log_f)

    # Make output table
    write_output_tbl(snp_info, geno_info, args.outdir)

    # Make plot ???

    # Close outputs
    log_f.write(f'\n{PROG} finished on {now()}\n')
    log_f.close()

# Run Code
if __name__ == '__main__':
    main()
