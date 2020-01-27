# pwd
cd /home/rad/users/gaurav/projects/gjsrmap

# # Template: 
# sh 06_run_ncrna_mapping_usage.sh hsa input/fastq/test/v1hsa output/test/v1hsa
# sh 06_run_ncrna_mapping_usage.sh mmu input/fastq/test/v2mmu output/test/v2mmu

# Git tag
git tag
# v1.2.0
git show v1.2.0

git tag -a v1.9 -m "Standalone version for a single server" 

###################################################################
# 2) RUPERT
# 02.1.1) Run mapping for Rupert's AML samples (Judith Hacker) to allncrnas on hsa/hg19 
# The samples are prepared using NEB kit.
# Adapter from rupert: AGATCGGAAGAGCACACGTCTGAACTCCAGT
# Copy the samples
projdir="/media/rad/HDD1/smallrna/rupert/jhackerAML"
fastqdir="${projdir}/fastq"
# cd "/media/nas/raw/TUM_NextSeq/191218_NB501802_0234_AH2WNHBGXF/Data/Intensities/BaseCalls"
# ls *.gz | parallel --progress --eta -j 32 "rsync -arvRP {} ${fastqdir}"
# cd -
sh 06_run_ncrna_mapping_usage.sh hsa ${fastqdir} ${projdir} AGATCGGAAGAGCACACGTCTGAACTCCAGT

###################################################################


# AGATCGGAAGAGCACACGTCTGAACTCCAGT 
# \

# GGCTGGTCCGATGGTAGTGGGTTATCAGGACAGATCGGAAGAGCACACGTCTGAACTCCAGTCACTT

# GGCTGGTCCGATGGTAGTGGGTTATCAGAACTTAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC

# TTGCCTGAAGCTGATGATGAGTTAGATCGGAAGAGCCCACGTCTGAACGCCAGTCACTTAGGCAGGT