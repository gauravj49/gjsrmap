#!/bin/bash

# Set user defined environment variables
jobdir=/usr/users/gjain/bin/projects/map2ncrna

# Get the input parameters
genomeFastaFile=$1
outindexdir=$2

# Create index dir if it does not exists
mkdir -p ${outindexdir}

# Get the length of genome
genomeLength=$(infile=${genomeFastaFile} python - <<END
import os
print sum(len(l.rstrip()) for l in open(os.environ['infile']) if l[0] != ">")
END
)

# For smaller genome, calculate genomeSAindexNbases parameter
# stargenomeSAindexNbases = min(14, log2(GenomeLength)/2 - 1)
log2=$(echo "(l($genomeLength)/l(2))/2 - 1" | bc -l)
# convert to integer
log2=$( printf "%.0f" $log2 )
comp=$(echo "$log2 < 14" | bc)
stargenomeSAindexNbases=$(($comp==0?14:$log2))

# Generate STAR index for the genome
# Note: Very important parameter for small genome: --genomeSAindexNbases 2
echo "Final command ran:"
echo " "
echo "/usr/users/gjain/bin/libs/STAR/bin/STAR --runMode genomeGenerate --genomeSAindexNbases ${stargenomeSAindexNbases} --genomeDir ${outindexdir} --genomeFastaFiles ${genomeFastaFile} --runThreadN 7 --genomeChrBinNbits 12"
/usr/users/gjain/bin/libs/STAR/bin/STAR --runMode genomeGenerate --genomeSAindexNbases ${stargenomeSAindexNbases} --genomeDir ${outindexdir} --genomeFastaFiles ${genomeFastaFile} --runThreadN 7 --genomeChrBinNbits 12

# Move the log file to the output folder
mv Log.out ${outindexdir}

