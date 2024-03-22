#!/usr/bin/env python3
import sys, os, argparse, datetime, gzip

#
# Globals
#
DATE = datetime.datetime.now().strftime("%Y%m%d")
PROG = sys.argv[0].split('/')[-1]
DESC = 'Extract the records that match a target probe/barcode sequence in a given FASTQ file.'

class FastqRecord:
    '''A sequence record from an FASTQ file'''
    def __init__(self, seq_id, sequence, plus, phred):
        assert len(sequence) == len(phred)
        assert not None in [seq_id, sequence, plus, phred]
        self.id    = seq_id
        self.seq   = sequence
        self.plus  = plus
        self.phred = phred
    def __str__(self):
        row = f'{self.id} {self.seq[:10]}... {self.plus} {self.phred[:10]}...'
        return row
    def print_record(self):
        record_str = f'{self.id}\n{self.seq}\n{self.plus}\n{self.phred}\n'
        return record_str


def parse_args():
    '''Command line arguments'''
    p = argparse.ArgumentParser(prog=PROG, description=DESC)
    p.add_argument('-f', '--fastq', required=True,
                   help='(str) Path to target FASTQ file.')
    p.add_argument('-t', '--target-seq', required=True,
                   help='(str) Target barcode/probe sequence to search in reads.')
    p.add_argument('-o', '--outdir', required=False, default='.',
                   help='(str) Path to output directory.')
    p.add_argument('-n', '--n-mismatch', required=False, default=1, type=int,
                   help='(int) Number of allowed sequence mismatched when searching target sequence.  [default=1]')
    # Check input arguments
    args = p.parse_args()
    args.outdir = args.outdir.rstrip('/')
    if not os.path.exists(args.outdir):
        sys.exit(f"Error: '{args.outdir}' not found.")
    if not os.path.exists(args.fastq):
        sys.exit(f"Error: '{args.fastq}' not found.")
    return args


def now():
    '''Print the current date and time.'''
    return f'{datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")}'


def parse_fastq(fastq_f, target_seq, out_f, n_mismatch=1):
    '''Open and parse the FASTQ file'''
    fh = open(fastq_f, 'r')
    if fastq_f.endswith('.gz'):
        fh = gzip.open(fastq_f, 'rt')
    print(f'Processing FASTQ file...', flush=True)
    records = 0
    kept = 0
    # FASTQ elements
    seq_id = None
    sequence = None
    plus = None
    phred = None
    # Parse file
    for i, line in enumerate(fh):
        line = line.strip('\n')
        if i%4 == 0:
            seq_id = line
            if not line.startswith('@'):
                sys.exit(f'Error: missing FASTQ header line {i}')
        elif i%4 == 1:
            sequence = line
        elif i%4 == 2:
            plus = line
            if not line.startswith('+'):
                sys.exit(f'Error: missing FASTQ plus line {i}')
        else:
            phred = line
            if len(phred) != len(sequence):
                sys.exit(f'Length of PHRED line {i} does not match sequence in line {i-2}')
            # Process the FASTQ record
            fastq_record = FastqRecord(seq_id, sequence, plus, phred)
            records += 1
            # Find matches to target sequence
            seq_match = find_target_sequence_match(fastq_record.seq, target_seq, n_mismatch)
            # Process if there is a match
            if seq_match:
                kept += 1
                out_f.write(fastq_record.print_record())
    print(f'    Read {records:,} FASTQ records.\n    Retained {kept:,} records matching the target sequence.', flush=True)
    out_f.close()
    fh.close()


def rev_comp(sequence):
    '''Reveser compliment a sequence'''
    rev = []
    for nt in sequence.upper():
        if nt == 'A':
            rev.append('T')
        elif nt == 'C':
            rev.append('G')
        elif nt == 'G':
            rev.append('C')
        elif nt == 'T':
            rev.append('A')
        elif nt not in ['A', 'C', 'G', 'T']:
            rev.append('N')
    rev_seq = ''.join(rev[::-1])
    return rev_seq


def set_output_fh(fastq_f, target, outdir='.'):
    '''Prepare file handle for output FASTQ'''
    # TODO: This search likely has to be more robust
    in_file_n = fastq_f.split('/')[-1].split('.')[0]
    out_file_n = f'{outdir}/{in_file_n}.{target}.fq.gz'
    fh = gzip.open(out_file_n, 'wt')
    return fh


def compare_seqs(query_seq, target_seq, n_mismatch=1):
    '''Find if two sequence are under N mismatches appart'''
    assert len(query_seq) == len(target_seq)
    obs_mis = 0
    matching = True
    for i in range(len(query_seq)):
        if query_seq[i] != target_seq[i]:
            obs_mis += 1
        if obs_mis > n_mismatch:
            matching = False
            break
    return matching


def find_target_sequence_match(query_seq, target_seq, n_mismatch=1):
    '''Look for the presence of target sequence in query sequence'''
    assert len(query_seq) >= len(target_seq)
    # find the block of sequence to search in query seq
    query_block = query_seq[0:len(target_seq)]
    # Search target in the forward orientation
    matching = compare_seqs(query_block, target_seq, n_mismatch)
    # Search target in the reverse orientation
    rev_target = rev_comp(target_seq)
    rev_matching = compare_seqs(query_block, rev_target, n_mismatch)
    # Return True if found in either forward or reverse complement
    if matching or rev_matching:
        return True
    else:
        return False





def main():
    print(f'{PROG} started on {now()}\n')
    args = parse_args()

    # Report to log
    print(f'Looking for records containing target sequence \'{args.target_seq}\' in input FASTQ:\n    {args.fastq}\n', flush=True)

    # Create output fh
    out_fh = set_output_fh(args.fastq, args.target_seq, args.outdir)

    # Parse input FASTQ
    parse_fastq(args.fastq, args.target_seq, out_fh, args.n_mismatch)


    # Close outputs
    print(f'\n{PROG} finished on {now()}')


# Run Code
if __name__ == '__main__':
    main()