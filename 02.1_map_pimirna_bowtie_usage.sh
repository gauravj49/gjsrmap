#!/bin/bash

# # Load necessary modules
# module load CUTADAPT/1.16
# module load BEDTOOLS/2.24
# module load BOWTIE1/1.1.2
# module load BOWTIE2/2.2.5
# module load FASTQC/0.11.4
# module load SAMTOOLS/1.9

# Get input arguments
jobdir=$1         # "path_to_project_directory"
bowtieindexdir=$2 # "path_to_bowtie_index"
fastqdir=$3       # "path_to_fastq_files"
outputdir=$4      # "path_to_output_directory"
reffastafile=$5   # "path_to_reference_genome/dna/hsa_dna_mirna_unique.fa"
ifTRIM=${6:-1}    # If you want to trim the reads (default is yes)
threePadapter=${7:-"TGGAATTCTCGGGTGCCAAGG"}
queue=${8:-"fat"}  # mpi or fat

# Set up relevant filenames and directories
bamdir=${outputdir}/"bam/iteration1"
samdir=${outputdir}/"sam/iteration1"
countsdir=${outputdir}/"counts"
trimmedFastqDir=${outputdir}/"trimmed_fastq/iteration1"
unmappedFastqDir=${outputdir}/"unmapped_fastq/iteration1"
mappingLogsDir="${outputdir}/mappingLogs/iteration1"
fastqcDirI0="${outputdir}/fastqc/iteration0"
fastqcDirI1="${outputdir}/fastqc/iteration1"
scriptsdir=${jobdir}/scripts/02.1_map_pimirna/$(basename ${outputdir})

# Create necessary dirs
mkdir -p ${samdir} ${bamdir} ${countsdir} ${scriptsdir} ${trimmedFastqDir} ${unmappedFastqDir} ${fastqcDirI0} ${fastqcDirI1}

for f in ${fastqdir}/*.fastq.gz
do 
 # Get the fastq file names
 ffname=$(basename ${f} .fastq.gz)
 fastqFile=${trimmedFastqDir}/${ffname}.fastq.gz

 # Get the jobname to submit for each job
 jobname="1_$(basename ${outputdir})_${ffname}"

 # Set up output folders
 bowtieoutputdir=${mappingLogsDir}/${ffname}
 errorsFile="${bowtieoutputdir}/${ffname}_iteration1_mapping.err"
 stdoutFile="${bowtieoutputdir}/${ffname}_iteration1_mapping.out"
 
 # Create necessary dirs
 mkdir -p ${bowtieoutputdir} ${scriptsdir}
   
 # Initialize a script file
 script="${scriptsdir}/${ffname}.sh"
 touch ${script}
 echo "#!/bin/bash" > "${script}"
 echo "" >> "${script}"

 # Get fastQC analysis before trimming
 echo "# Get fastQC analysis before any trimming" >> "${script}"
 echo "fastqc ${f} -o ${fastqcDirI0} -q" >> "${script}"
 echo "" >> "${script}"

 # Trim the adapters and save the output file
   #--minimum-length N or -m N: Use this to throw away processed reads shorter than N bases.
   #--maximum-length N or -M N: Use this to throw away processed reads longer than N bases.
   #--too-long-output FILE    : Instead of throwing away the reads that are too long (according to -M), write them to FILE (in FASTA/FASTQ format).
   #--trim-qualities: Quality trimming (-q) parameter can be used to trim low-quality ends from reads before adapter removal. For this to work correctly, the quality values must be encoded as ascii(phred quality + 33). If they are encoded as ascii(phred quality + 64), you need to add --quality-base=64 to the command line.

 if [ ${ifTRIM} -ne "0" ]; then
  echo "# Trimming the files" >> "${script}"
  echo "cutadapt -q 28 -a ${threePadapter} -m 16 -M 33 --too-long-output ${unmappedFastqDir}/${ffname}_greater_than_33bp_long.fastq -o ${fastqFile} ${f}" >> "${script}"
  echo "" >> "${script}"
 else
  echo "Trimming is turned off"
  echo "# Trimming is turned off" >> "${script}"
  echo "" >> "${script}"
  fastqFile=${f}
 fi

 # Get fastQC analysis after trimming
 echo "# Get fastQC analysis after trimming" >> "${script}"
 echo "fastqc ${fastqFile} -o ${fastqcDirI1} -q" >> "${script}"
 echo "" >> "${script}"

 # Unzip fastq.gz file before mapping with bowtie
 echo "# Unzip fastq.gz file before mapping with bowtie" >> "${script}"
 echo "gunzip ${fastqFile}" >> "${script}"
 echo "" >> "${script}"

 # Get the new fastq variable
 fastqFile=${trimmedFastqDir}/${ffname}.fastq

 # Get the output file names
 bamfile="${bamdir}/${ffname}.bam"
 samfile="${samdir}/${ffname}.sam"
 unmapped_fastq="${unmappedFastqDir}/${ffname}.unmapped.fastq"

 # Mapping short reads to the reference genome using bowtie aligner
 # Bowtie command
 # -q                 query input files are FASTQ .fq/.fastq (default)
 # -v <int>           report end-to-end hits w/ <=v mismatches; ignore qualities; No more than v mismatches in the entire length of the read
 # -p/--threads <int> number of alignment threads to launch (default: 1)
 # -S/--sam           write hits in SAM format
 # -k <int>           report up to <int> good alignments per read (default: 1)
 # -a/--all           report all alignments per read (much slower than low -k)
 # -m <int>           suppress all alignments if > <int> exist (def: no limit)
 #                    -m 3 instructs bowtie to refrain from reporting any alignments for reads having more than 3 reportable alignments. The -m option
 #                    is useful when the user would like to guarantee that reported alignments are "unique", for some definition of unique.
 echo "# Mapping short reads to the pimiRNA as reference genome using bowtie aligner" >> "${script}"
 echo "bowtie -q -v 0 -p 4 -S -a -m 1 --un ${unmapped_fastq} ${bowtieindexdir}/$(basename ${bowtieindexdir}) ${fastqFile} ${samfile} " >> "${script}"
 echo "" >> "${script}"

 # Zip unmapped_fastq files for next iteration
 echo "# Zip unmapped_fastq files for next iteration" >> "${script}"
 echo "gzip ${unmapped_fastq}" >> "${script}"
 echo "gzip ${unmappedFastqDir}/${ffname}_greater_than_33bp_long.fastq" >> "${script}"
 echo "" >> "${script}"

 # Get the BAM files from the SAM files
 echo "# Get the BAM files from the SAM files" >> "${script}"
 echo "samtools view -Sb ${samfile} > ${bamfile}" >> "${script}"
 echo "" >> "${script}"

 # Sort the bam file before creating the index
 echo "# Sort the bam file before creating the index" >> "${script}"
 echo "samtools sort -o ${bamfile}.sorted ${bamfile}" >> "${script}"
 echo "mv ${bamfile}.sorted ${bamfile}" >> "${script}"
 echo "" >> "${script}"

 # Create the index of the bam file
 echo "# Create the index of the bam file" >> "${script}"
 echo "samtools index ${bamfile}" >> "${script}"
 echo "" >> "${script}"

 # Remove the sam files
 echo "# Remove the sam files" >> "${script}"
 echo "rm -rf ${samfile} " >> "${script}"
 echo "" >> "${script}"

 chmod 775 "${script}"
 echo "Submitting job for ${ffname} ..." 
 sh "${script}"
 echo " "
done

