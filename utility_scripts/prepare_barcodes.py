#!/usr/bin/env python3
import sys, os, argparse, datetime

# (c) 2024, Angel G. Rivera-Colon

#
# Globals
#
DATE = datetime.datetime.now().strftime("%Y%m%d")
PROG = sys.argv[0].split('/')[-1]
DESC = 'Prepare the set of per-individual barcodes for GTseq analysis.'
ROWS = set('ABCDEFGH')   # Rows for a 96-well plate
COLS = set(range(1,13))  # Columns for a 96-well plate
NUCS = set('ACGT')       # Valid nucleotides for barcodes

def parse_args():
    '''Command line arguments'''
    p = argparse.ArgumentParser(prog=PROG, description=DESC)
    p.add_argument('-p', '--plate-barcodes', required=True,
                   help='(str) Path to the per-plate barcodes file.')
    p.add_argument('-w', '--well-barcodes', required=True,
                   help='(str) Path to the per-well barcodes file.')
    p.add_argument('-m', '--sample-map', required=True,
                   help='(str) Path to the multi-plate sample map file.')
    p.add_argument('-o', '--outdir', required=False, default='.',
                   help='(str) Path to output directory.')
    # Check input arguments
    args = p.parse_args()
    args.outdir = args.outdir.rstrip('/')
    if not os.path.exists(args.outdir):
        sys.exit(f"Error: '{args.outdir}' not found.")
    if not os.path.exists(args.plate_barcodes):
        sys.exit(f"Error: '{args.plate_barcodes}' not found.")
    if not os.path.exists(args.well_barcodes):
        sys.exit(f"Error: '{args.well_barcodes}' not found.")
    if not os.path.exists(args.sample_map):
        sys.exit(f"Error: '{args.sample_map}' not found.")
    return args

def now():
    '''Print the current date and time.'''
    return f'{datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")}'

def is_nucleotide(seq):
    '''Check if a string contains legal nucleotides.'''
    return set(seq).issubset(NUCS)

def parse_plate_barcodes(plate_barcodes_f):
    '''Processing the plate barcodes'''
    print('Extracting the plate barcodes...')
    plate_barcodes = dict()
    with open(plate_barcodes_f) as fh:
        for i, line in enumerate(fh):
            if line.startswith('#'):
                continue
            fields = line.strip('\n').split('\t')
            assert len(fields) == 2
            plate_id = fields[0]
            barcode  = fields[1].upper()
            if not is_nucleotide(barcode):
                sys.exit(f"Error: not a legal barcode ({plate_id} {barcode}) in line {i+1}.")
            plate_barcodes[plate_id] = barcode
    print(f'    Extracted {len(plate_barcodes):,} barcodes from the plate file.')
    return plate_barcodes

def parse_well_barcodes(well_barcodes_f):
    '''Extract the well barcodes.'''
    print('\nProcessing the well barcodes...')
    well_barcodes = dict()
    records = 0
    with open(well_barcodes_f) as fh:
        for i, line in enumerate(fh):
            if line.startswith('#'):
                continue
            fields = line.strip('\n').split('\t')
            assert len(fields) == 3
            legal = True
            # Rows go from A-H
            row = fields[0]
            if row not in ROWS:
                legal = False
            # Columns go from 1-12
            column = int(fields[1])
            if column not in COLS:
                legal = False
            # Barcode must be nucleotides
            barcode = fields[2].upper()
            if not is_nucleotide(barcode):
                legal = False
            # Check if legal fields
            if not legal:
                sys.exit(f"Error: not a legal barcode ({row} {column} {barcode}) in line {i+1}.")
            # Populate dict
            cols = [ None for _ in range(12) ]
            well_barcodes.setdefault(row, cols)
            well_barcodes[row][column-1] = barcode
            records += 1
    print(f'    Extracted {records:,} barcodes from the well file.')
    return well_barcodes

def parse_sample_map(sample_map_f, plate_barcodes, well_barcodes, outdir='.'):
    '''Process the sample map and link to the plate and well barcodes.'''
    outfh = open(f'{outdir}/barcodes.tsv', 'w')
    print('\nProcessing the sample map...')
    with open(sample_map_f) as fh:
        curr_plate = None
        total_samples = 0
        plate_samples = 0
        for i, line in enumerate(fh):
            if line.startswith('#'):
                continue
            fields = line.strip('\n').split('\t')
            # Check if a plate header and process
            if fields[0].startswith('plate_'):
                if curr_plate is not None:
                    print(f'    Extracted {plate_samples} samples in {curr_plate}.')
                    plate_samples = 0
                curr_plate = fields[0]
                if curr_plate not in plate_barcodes:
                    sys.exit(f"Error: plate ID '{curr_plate}' not in plate barcodes file (line {i+1}).")
            else:
                # The other lines should have the row IDs
                row = fields[0]
                if row not in ROWS:
                    sys.exit(f"Error: '{row}' not a valid plate row (line {i+1}).")
                # Determine the plate barcode from the current plate
                plate_barcode = plate_barcodes[curr_plate]
                # Determine the well barcode for each sample in the row
                for j, sample in enumerate(fields[1:]):
                    well_barcode = well_barcodes[row][j]
                    # Export to file
                    outfh.write(f'{plate_barcode}\t{well_barcode}\t{sample}\n')
                    total_samples += 1
                    plate_samples += 1
    print(f'    Extracted {plate_samples} samples in {curr_plate}.')
    print(f'\nFormatted barcodes for {total_samples:,} total samples.')
    outfh.close()

def main():
    print(f'{PROG} started on {now()}\n')
    args = parse_args()
    # Extract the plate barcodes
    plate_bars = parse_plate_barcodes(args.plate_barcodes)
    # Extract the well barcodes
    well_bars = parse_well_barcodes(args.well_barcodes)
    # Process the sample maps
    parse_sample_map(args.sample_map, plate_bars, well_bars, args.outdir)
    # Close outputs
    print(f'\n{PROG} finished on {now()}')


# Run Code
if __name__ == '__main__':
    main()