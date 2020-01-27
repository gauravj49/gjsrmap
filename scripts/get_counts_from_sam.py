#!/usr/local/bin/python
"""
***********************************************
- PROGRAM: get_counts_from_sam.py
- CONTACT: Gaurav Jain(gaurav.jain@dzne.edu)
***********************************************
"""
print (__doc__)

# Built in modules
import argparse
import os.path
import sys

# 3rd party modules
import textwrap
import re
import subprocess
import math
from collections import *

################ USER CONFIGURATION ###################
#######################################################

def main():
    # Get input options
    args             = check_options()
    input_file       = args.input_file
    input_fasta_file = args.input_fasta_file
    output_file      = args.output_file
    normoutput_file  = "{0}_normalized.txt".format(get_file_info(output_file)[3])
    
    ##########################################
    ## NOTE: allMirName is actually allNcrnaNames
    # Get the list of all mirnames from reference fasta file
    print("- Get the list of all mirnas from the reference fasta file: {0} ....".format(input_fasta_file))
    allSmncrnaNames = get_all_mir_names(input_fasta_file)
    print("\t- Total number of smallrnas: {0}".format(len(allSmncrnaNames)))

    # Read and parse BWA SAM file
    # @SQ	SN:mmu-miR-6998-5p	LN:22
    # HWI-ST1140:186:C77GGACXX:1:2316:19229:100797	0	mmu-miR-486b-5p	1	0	22M	*	0	0	TCCTGTACTGAGCTGCCCCGAG	CCCFFEFFHHHHHJJJJIJJJJ	XT:A:R	NM:i:0	X0:i:2	X1:i:0	XM:i:0	XO:i:0	XG:i:0	MD:Z:22	XA:Z:mmu-miR-486a-5p,+1,22M,0;

    # Initialize the dictionary
    samSmncrnaNames = defaultdict(int)
    
    # Get the counts
    print("\n- Getting counts from the SAM file (may take some time depending on sam file size)...")
    process = subprocess.Popen("grep -v \"^@\" {0} | cut -f3 | sort | uniq -c | column -t ".format(input_file), shell=True, stdout=subprocess.PIPE, encoding='utf-8')
    # print(process.communicate()[0])
    stdout_list = process.communicate()[0].split('\n')
    stdout_list = [_f for _f in stdout_list if _f] # filter out empty lines
    print("\t- Counting finished")

    # Get the dictionary of counts found in sam file
    for line in stdout_list:
        row = line.split()
        if len(row) < 2:
            continue
        else:
            mirCounts = row[0]
            featureid = row[1]
            samSmncrnaNames[featureid] = int(mirCounts)


    # Get library depth
    library_depth = int(sum(samSmncrnaNames.values()))
    print("\n- Process the input SAM file: {0}".format(input_file))
    print("\t- Library depth: {0}".format(library_depth))

    # Get the list of mirnas not in SAM file
    unexpressed_smncrnaNames = list(set(allSmncrnaNames).difference(set(samSmncrnaNames.keys())))
    no_feature_names     = list(set(unexpressed_smncrnaNames + list(samSmncrnaNames.keys())).difference(set(allSmncrnaNames)))
    
    no_features_count = 0
    for nf in no_feature_names:
        no_features_count += samSmncrnaNames[nf]

    # Save the counts to output file
    outf  = open(output_file , 'w')
    outnf = open(normoutput_file , 'w')

    # Get the combined dictionary of all the smallncrnaNames: expressed_smallncrnas + unexpressed_smallncrnas
    allsmncrnadict = samSmncrnaNames.copy() 
    for k in set(unexpressed_smncrnaNames).difference(set(no_feature_names)):
        allsmncrnadict[k] = 0

    # Save counts in output file sorted by small ncrna names (so that they are same in every samples)
    # This is useful when doing other analysis like DE via DeSEQ2
    i = 1
    for k,v in sorted(allsmncrnadict.items()):
        if k == '*':
            continue
        if k in no_feature_names:
            continue
        outf.write("{0}\t{1}\n".format(k,v))
        outnf.write("{0}\t{1}\n".format(k, int(math.ceil(v*1000000.0/library_depth))))
        i += 1
        
    # Print the last 5 lines to make the counts file in htseq format
    # __no_feature, __ambiguous, __too_low_Qual, __not_aligned, __alignment_not_unique
    outf.write("__no_feature\t{0}\n__ambiguous\t0\n__too_low_Qual\t0\n__not_aligned\t0\n__alignment_not_unique\t0\n".format(no_features_count))
    outnf.write("__no_feature\t{0}\n__ambiguous\t0\n__too_low_Qual\t0\n__not_aligned\t0\n__alignment_not_unique\t0\n".format(no_features_count))
      
    # Close the output file handles
    outf.close()
    outnf.close()

################ USER DEFINED FUNCTIONS ###################
def get_all_mir_names(input_file):
    ''' Read the input fasta file and create a dictionary out of it.'''
 
    # Inititalize the list
    d = list()

    # Get the information in a dictionary
    with open(input_file,'rU') as fg:
        for name, seq in read_fasta(fg):
            d.append(name.lstrip('>'))
    return d

def read_fasta(fp):
    ''' Reads the fasta file and returns the name and seq '''
    name, seq = None, []
    for line in fp:
        line = line.rstrip()
        if line.startswith(">"):
            if name: yield (name, ''.join(seq))
            name, seq = line, []
        else:
            seq.append(line)
    if name: yield (name, ''.join(seq))

def get_file_info(file_name_with_path):
    ''' Get path, basename, ext, path+basename and basename+ext of a file '''

    # get the path of the file
    input_file_path = os.path.dirname(os.path.abspath(file_name_with_path))

    # get the filename without path but with extension
    input_file_name = os.path.basename(file_name_with_path)

    # get basename of the file (no path or extension)
    input_file_basename = os.path.splitext(input_file_name)[0]

    # get extension of the file
    input_file_extension = os.path.splitext(input_file_name)[1]
    # remove "."
    input_file_extension = re.sub(r"^\.","",input_file_extension)

    # get the path+filename
    path_filename = input_file_path + "/" + input_file_basename

    return (input_file_path, input_file_basename, input_file_extension, path_filename, input_file_name)

def create_dir(dirname):
    '''Creates a directory if not exists. Returns relevant status messages on error'''
    try:
        os.makedirs(dirname)
    except OSError:
        if os.path.exists(dirname):
            # We are nearly safe
            pass
        else:
            # There was an error on creation, so make sure we know about it
            raise

def useful_lines(fh):
    '''Skips all the comments (right now #) from the file and returns the generator object of lines without comments.
    It remembers the state where it left. So it will after the last executed line '''
    for line in (l.strip() for l in fh if not l.startswith('#')):
        yield line

def print_initial_arguments(parser):
    ''' Prints all the initial arguments entered '''

    print("\n------Input Arguments ------")
    opts = vars(parser.parse_args())
    maxl = max(len(text) for text in list(opts.keys()))
    for k,v in list(opts.items()):
        print("%-*s = %s" % (maxl, k, v))
    print("-"*29, "\n")

# class for Logging
class Log(object):
    def __init__(self, *files):
        self.files = files

    def write(self, obj):
        for f in self.files:
            f.write(obj)

    def flush(self):
        pass

def check_options():
    ''' Checks the options to the program '''

    # Create parser object
    parser = argparse.ArgumentParser(add_help=True,
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=textwrap.dedent('''\
        ----------------- SAMPLE USAGE ------------------
        - python scripts/get_counts_from_sam.py -if=input/reference_genome/hsa_unique_mirna.fa -is=input/sam/mapped_only_mouse.sam -of=output/mapped_only_mouse_counts.txt
        -------------------------------------------------
        CONTACT: 
        	Gaurav Jain
        	gaurav.jain@dzne.de
        -------------------------------------------------
        '''))

    # Add arguments 
    parser.add_argument("-is", metavar="--insamf", help="*Input SAM file"    , dest="input_file"  , type=str, required=True)
    parser.add_argument("-if", metavar="--infile", help="*Input fasta file"  , dest="input_fasta_file"  , type=str, required=True)
    parser.add_argument("-of", metavar="--opfile", help="*Output file name"  , dest="output_file" , type=str, required=True)
    
    # Print the help message only if no arguments are supplied
    if len(sys.argv)==1:
        parser.print_help()
        sys.exit(1)

    # Save the STDOUT output in a log file
    if parser.parse_args().output_file:
        logdir = "{0}/logs".format(get_file_info(parser.parse_args().output_file)[0])
        create_dir(logdir)
        logfile = "{0}/{1}.log".format(logdir, get_file_info(parser.parse_args().output_file)[1])
    else:
        logdir  = "{0}/logs".format(os.getcwd())
        create_dir(logdir)
        logfile = "{0}/{1}.log".format(logdir,get_file_info(sys.argv[0])[1])

    logf = open(logfile, 'w')
    sys.stdout = Log(logf, sys.stdout)

    # Parse command line with parse_args and store it in an object
    args = parser.parse_args()
    print_initial_arguments(parser)
    return args

# main function
if __name__=="__main__":
      main()
