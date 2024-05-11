#!/usr/bin/env python3

import sys, os, gzip, argparse
from multiprocessing import Process

# GTseq_BarcodeSplit_KML.py is Kira Long's edited version of Nate Campbell's barcode splitting script for GTseq data

# Original notes from Nate Campbell
# GTseq_BarcodeSplit_MP.py
# by Nate Campbell
# Use this script to split out individual files by barcode combination rather than using "grep method".
# This version will process up to 5,000 individual samples and simultaneously run up to 10 processing cores for faster parsing of individual sequences.
# Input file 1 is a .csv file containing sample and barcode information [individual sample names, PlateID,i7_name,i7_sequence,i5_name,i5_sequence] for each sample.
# example: Sample,PlateID,i7_name,i7_sequence,i5_name,i5_sequence
#          Sample123,P1234,i001,ACCGTA,25,CCCTAA
#          Sample321,P1234,i001,ACCGTA,26,GGCACA
#          ....
# Note: The header line is ignored while executing the script so always include it to avoid missing data.  Also, note that output files are appended
# and therefore not overwritten if the script is called multiple times.  If for some reason the script needs to be run more than once for the same set of samples, 
# the original output files will need to be deleted to avoid having files with duplicated sequences.
# Input file 2 is the .fastq file containing the sequences.
# Output is a set of fastq files for all individuals included in the library named in the format [i7name]_[i5name]_[PlateID]_[SampleID].fastq 
# This script uses multiple processors (1 for every 500 samples) and was clocked at about 20 minutes to complete barcode splitting for one lane of data.
# Speed is faster than grep method while using fewer compute resources.
# Using this script for barcode splitting and the GTseq_Genotyper_v2.pl script for genotyping allows the GTseq pipeline to be run on a linux
# desktop computer.  (Raw data to genotypes in less than 1 hour).

# KML addition: usage example if command is run without additional arguments/options
PROG = sys.argv[0].split('/')[-1]
usage_message='''{} usage example:

    {} path/barcode.csv path/raw_reads.fastq(.gz) path/output/dir

Script uses 3 positional arguments:
     1) Path to the barcodes file (Needs plate and individual barcodes) in csv format (required).
     2) Path to the fastq file with the raw reads. It can be gzipped (required).
     3) Path to the output directory (defaults to current directory)'''.format(PROG,PROG)
if len(sys.argv)<2:
     print(usage_message)
     sys.exit()

# KML addition: Path to input csv file with barcode information
path1 = sys.argv[1]
if os.path.exists(path1) == False:
     sys.exit('Path to barcodes file not found')

# KML addition: Path to the fastq file of the sequencing data to demultiplex
path2 = sys.argv[2]
if os.path.exists(path2) == False:
     sys.exit('Path to fastq file not found')

# KML addition: Path to an output directory
path3 = '.'
if len(sys.argv) == 4:
     path3 = sys.argv[3]
path3 = path3.rstrip('/')

# To run this script, use cmd path/barcode.csv path/raw_reads_fastq path/output/dir

def split_file(individual_list):
     if individual_list == 'list1':
          individual_list = list1
     elif individual_list == 'list2':
          individual_list = list2
     elif individual_list == 'list3':
          individual_list = list3
     elif individual_list == 'list4':
          individual_list = list4
     elif individual_list == 'list5':
          individual_list = list5
     elif individual_list == 'list6':
          individual_list = list6
     elif individual_list == 'list7':
          individual_list = list7
     elif individual_list == 'list8':
          individual_list = list8
     elif individual_list == 'list9':
          individual_list = list9
     else:
          individual_list = list10

     # KML addition: Open the fastq file and check if it is zipped to open it either way
     fq = None
     if path2.endswith('.gz'):
          fq = gzip.open(path2, 'rt')
     else:
          fq = open(path2, 'r')
     BC_Dict = {}
     File_Dict = {}
     handle_dict = {}

     for lines in individual_list:
          # KML addition: renaming variable stuff to fields
          # fields describes the different "fields" in the csv, example sample_name, barcode1, barcode2 etc.
          fields = lines.split(',')
          # KML addition: naming variables based on what they are in the input barcodes file for readability
          # Input barcode file header: Sample,PlateID,i7_name,i7_sequence,i5_name,i5_sequence
          sample_id = fields[0]
          plate_id = fields[1]
          i7_name = fields[2]
          i7_seq = fields[3]
          i5_name = fields[4]
          i5_seq = fields[5]

          # KML addition: changed variable name to output_fq
          # TODO: option for zipping outputs?
          output_fq = '{}/{}_{}_{}_{}.fastq'.format(path3,i7_name,i5_name,plate_id,sample_id)
          barcode_combo = '{}+{}'.format(i7_seq,i5_seq)
          BC_Dict[output_fq] = barcode_combo
          File_Dict[barcode_combo] = output_fq
          handle_dict[barcode_combo] = open(File_Dict[barcode_combo], 'a')

     lineNo2 = 0
     writelines = 0
     seqID = "NA"

     for line in fq:
          lineNo2 = lineNo2 + 1
          if lineNo2 == 1:
               # KML addition: changed things to fq_line
               fq_line = line.split(':')
               seqID = fq_line[0]
          if lineNo2 > 1:
               break

     for line in fq:
          # Selects the 4 lines in the fastq record
          lineNo2 = lineNo2 + 1
          if writelines < 5 and writelines > 0:
               f_out.write(line)
               writelines = writelines + 1
               if writelines == 4:
                    writelines = 0
          # check if the fastq file has the barcode
          # KML addition: changed info to fq_header
          if seqID in line:
               fq_header = line.split(':')
               BC = fq_header[9]
               if BC in File_Dict:
                    f_out = handle_dict[BC]
                    f_out.write(line)
                    writelines = 1
     fq.close()
     return

def Main():
     if len(list1) > 0:
          P1 = Process(target=split_file, args=('list1',))
          print('Process initiated for 1st set of samples')
          P1.start()
     if len(list2) > 0:
          P2 = Process(target=split_file, args=('list2',))
          print('Process initiated for 2nd set of samples')
          P2.start()
     if len(list3) > 0:
          P3 = Process(target=split_file, args=('list3',))
          print('Process initiated for 3rd set of samples')
          P3.start()
     if len(list4) > 0:
          P4 = Process(target=split_file, args=('list4',))
          print('Process initiated for 4th set of samples')
          P4.start()
     if len(list5) > 0:
          P5 = Process(target=split_file, args=('list5',))
          print('Process initiated for 5th set of samples')
          P5.start()
     if len(list6) > 0:
          P6 = Process(target=split_file, args=('list6',))
          print('Process initiated for 6th set of samples')
          P6.start()
     if len(list7) > 0:
          P7 = Process(target=split_file, args=('list7',))
          print('Process initiated for 7th set of samples')
          P7.start()
     if len(list8) > 0:
          P8 = Process(target=split_file, args=('list8',))
          print('Process initiated for 8th set of samples')
          P8.start()
     if len(list9) > 0:
          P9 = Process(target=split_file, args=('list9',))
          print('Process initiated for 9th set of samples')
          P9.start()
     if len(list10) > 0:
          P10 = Process(target=split_file, args=('list10'))
          print('Process initiated for 10th set of samples')
          P10.start()
     print('Keepin it 100!  Your files will be ready shortly...')


list1 = []
list2 = []
list3 = []
list4 = []
list5 = []
list6 = []
list7 = []
list8 = []
list9 = []
list10 = []

individuals = 0
start = 1  # skip line 1 of input file (header line)...
end = 0
f = open(path1, 'r')
lineNo = 0

# Determine the number of individuals in the library...
for line in f:
     lineNo = lineNo + 1
file_end = lineNo
individuals = lineNo - 1
if individuals > 500 + start:
     end = 500 + start
else:
     end = file_end

sets = 0
Samples = True
starting_ind = 1

while Samples == True:
     lineNo = 0

     f = open(path1, 'r')
     list0 = []     # set empty list...
     for line in f:
          lineNo = lineNo + 1
          if lineNo > start and lineNo <= end: # Populate list in sets of 500 individuals ...
               list0.append(line)
     sets = sets + 1
     print(sets)
     size = len(list0)
     if sets == 1:
          list1 = tuple(list0)
     elif sets == 2:
          list2 = tuple(list0)
     elif sets == 3:
          list3 = tuple(list0)
     elif sets == 4:
          list4 = tuple(list0)
     elif sets == 5:
          list5 = tuple(list0)
     elif sets == 6:
          list6 = tuple(list0)
     elif sets == 7:
          list7 = tuple(list0)
     elif sets == 8:
          list8 = tuple(list0)
     elif sets == 9:
          list9 = tuple(list0)
     else:
          list10 = tuple(list0)

     ending_ind = end - 1

     if end == file_end:
          Samples = False
          break
     elif start > file_end:
          Samples = False
          break

     start = end
     end = start + 500
     starting_ind = starting_ind + 500

     if end > file_end:
          end = file_end
          ending_ind = file_end - 1
     #print('Next set of individuals begins at line %s and ends at %s of input file' % (start, end))

     f.close()

if __name__ == '__main__':
     Main()
