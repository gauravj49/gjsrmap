#!/usr/local/bin/python
"""
***********************************************
- PROGRAM: process_fasta.py
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
import numpy as np
import scipy as sp
import pandas as pd
import matplotlib as mp
#mp.use('Agg') # to use matplotlib without X11
import matplotlib.pyplot as plt
import subprocess
import binascii as bi
import scipy.stats as stats
from collections import *
from numpy import nanmean

# for looping files in a dir
import glob

# user defined modules
from gjainLIB import *      # import all the functions from the Gaurav`s python library

### for color scale
from  matplotlib import colors
from itertools import cycle, islice # barplot colors

################ USER CONFIGURATION ###################
#######################################################

def main():
    # Get input options
    args = check_options()
    input_fasta_file = args.input_fasta_file
    output_file      = args.output_file
    ncrna_prefix     = args.ncrna_prefix
    ftDNA            = args.ftDNA
    
    # Get the dictionary of fasta file
    print "- Get fasta dictionary for the {0} file....".format(input_fasta_file)
    fasta_dict = get_fasta_dict(input_fasta_file, ncrna_prefix)

    # Save the counts to output file
    outf = open(output_file , 'w')

    # Concatenate the names for which there are duplicate sequences and save them in output file
    print "- Saving the processed fasta file..."
    for k,v in fasta_dict.iteritems():
        # Concatenate the names if there are more than one name for the same sequence
        if len(v) > 1:
            name = "|".join(list(OrderedDict.fromkeys(v))) # this removes duplicate keys from the list (as in pirna)
        else:
            name = v[0]

        # Save it in the fasta format
        if ftDNA:
            outf.write(">{0}\n{1}\n".format(name,k))
        else:
            outf.write(">{0}\n{1}\n".format(name,rna2dna(k)))

    print "- Your output fasta file is: {0}\n".format(output_file)
################ USER DEFINED FUNCTIONS ###################
def get_fasta_dict(input_file, ncrna_prefix):
    ''' Read the input fasta file and create a dictionary out of it.'''
 
    # Inititalize the dictionary
    d = defaultdict(list)

    # Get the information in a dictionary
    with open(input_file,'rU') as fg:
        for name, seq in read_fasta(fg):
            # Save only the first part of mirna name or other smallrna name
            ## 	name before: hsa-miR-518a-5p MIMAT0005457
            ## 	name after : hsa-miR-518a-5p
            ## 	name before: hsa_piR_000001|gb|DQ569913|Homo sapiens:21:43519790:43519815:Plus
            ## 	name after : hsa_piR_000001
            name = re.split(r'\s|\|', name)[0].lstrip('>')  # split on either | or space
            if ncrna_prefix:
                name = "{0}_{1}".format(ncrna_prefix, name)

            # To get the duplicate reads, make reads as keys and append names as values
            d[seq].append(name.lstrip('>'))

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

def rna2dna(rna):
    ''' Convert Us to Ts '''
    mapping = string.maketrans('uU', 'tT')
    return rna.translate(mapping)

def check_options():
    ''' Checks the options to the program '''

    # Create parser object
    parser = argparse.ArgumentParser(add_help=True,
        formatter_class=argparse.RawTextHelpFormatter,
        epilog=textwrap.dedent('''\
        ----------------- SAMPLE USAGE ------------------
        - python scripts/process_fasta.py -if=input/reference_genome/rna/hsa_rna_mirna.fa -of=input/reference_genome/dna/hsa_dna_mirna_unique.fa
        - python scripts/process_fasta.py -if=input/reference_genome/dna/hsa_dna_mirna.fa -of=input/reference_genome/dna/hsa_dna_mirna_unique.fa -dn
        -------------------------------------------------
        CONTACT: 
        	Gaurav Jain
        	gaurav.jain@dzne.de
        -------------------------------------------------
        '''))

    # Add arguments 
    parser.add_argument("-if", metavar="--infile", help="*Input fasta file", dest="input_fasta_file", type=str, required=True)
    parser.add_argument("-of", metavar="--opfile", help="*Output file name", dest="output_file"     , type=str, required=True)
    parser.add_argument("-px", metavar="--prefix", help=" Prefix to add before a name.\n Example, add _premirna to mmu-let-7j to make premirna_mmu-let-7j", dest="ncrna_prefix", type=str, default="")
    parser.add_argument("-dn",         "--fstDNA", help=" If set, treat fasta file as dna else rna" , dest="ftDNA" , action='store_true', default=False)

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
