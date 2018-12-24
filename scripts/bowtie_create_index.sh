#!/bin/bash

# Unload conflicting modules
module unload BOWTIE1/1.0.0

# Load necessary modules
module load BOWTIE1/1.1.2

# Set user defined environment variables
jobdir=/usr/users/gjain/bin/projects/map2ncrna

#Usage: bowtie-build [options]* <reference_in> <ebwt_outfile_base>
#reference_in            comma-separated list of files with ref sequences
#ebwt_outfile_base       write Ebwt data to files with this dir/basename

# Get the input parameters
genomeFastaFile=$1
outindexdir=$2

# Create index dir if it does not exists
mkdir -p ${outindexdir}

# Generate bowtie index for the genome
echo "Final command ran:"
echo " "
echo "bowtie-build ${genomeFastaFile} ${outindexdir}/$(basename ${outindexdir})"
bowtie-build ${genomeFastaFile} ${outindexdir}/$(basename ${outindexdir})


