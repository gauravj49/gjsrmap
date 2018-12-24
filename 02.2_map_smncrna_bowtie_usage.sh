#!/bin/bash

# Unload conflicting modules
module unload BOWTIE1/1.0.0
module unload SAMTOOLS/1.9
module unload SAMTOOLS/1.2

# Load necessary modules
module load CUTADAPT/1.16
module load BEDTOOLS/2.24
module load BOWTIE1/1.1.2
module load BOWTIE2/2.2.5
module load FASTQC/0.11.4
module load SAMTOOLS/1.9

# Get input arguments
jobdir=$1         # "path_to_project_directory"
bowtieindexdir=$2 # "path_to_bowtie_index"
fastqdir=$3       # "path_to_fastq_files"
outputdir=$4      # "path_to_output_directory"
reffastafile=$5   # "path_to_reference_genome/dna/hsa_dna_mirna_unique.fa"
ifTRIM=${6:-1}    # If you want to trim the reads (default is yes)
threePadapter=${7:-"TGGAATTCTCGGGTGCCAAGG"}
queue=${8:-"fat"} # mpi or fat

# Set up relevant filenames and directories
bamdir=${outputdir}/"bam/iteration2"
samdir=${outputdir}/"sam/iteration2"
countsdir=${outputdir}/"counts"
trimmedFastqDir=${outputdir}/"trimmed_fastq/iteration2"
unmappedFastqDir=${outputdir}/"unmapped_fastq/iteration2"
unmappedFastqDirI1=${outputdir}/"unmapped_fastq/iteration1"
mappingLogsDir="${outputdir}/mappingLogs/iteration2"
fastqcDirI2="${outputdir}/fastqc/iteration2"
scriptsdir=${jobdir}/scripts/02.2_map_smncrna/$(basename ${outputdir})

# Create necessary dirs
mkdir -p ${samdir} ${bamdir} ${countsdir} ${scriptsdir} ${trimmedFastqDir} ${unmappedFastqDir} ${fastqcDirI2}

# Initialize the counter to submit jobs as an array
i=0

for f in ${fastqdir}/*.fastq.gz
do 
 # Get the fastq file names
 ffname=$(basename ${f} .fastq.gz)
 fastqFile=${unmappedFastqDirI1}/${ffname}_greater_than_33bp_long.fastq

 # Get the jobname to submit each job as a job array
 # jobname="02.2_$(basename ${outputdir})_${ffname}"
 i=$((i + 1))
 jobname="2_$(basename ${outputdir})_${ffname}"
 prevjobname="1_$(basename ${outputdir})_${ffname}"

 # Set up output folders
 bowtieoutputdir=${mappingLogsDir}/${ffname}
 errorsFile="${bowtieoutputdir}/${ffname}_iteration2_mapping.err"
 stdoutFile="${bowtieoutputdir}/${ffname}_iteration2_mapping.out"
 
 # Create necessary dirs
 mkdir -p ${bowtieoutputdir} ${scriptsdir}
   
 # Initialize a script file
 script="${scriptsdir}/${ffname}.sh"
 touch ${script}
 echo "#!/bin/bash" > "${script}"
 echo "" >> "${script}"

 # Unzip fastq.gz file before mapping with bowtie
 echo "# Unzip fastq.gz file before mapping with bowtie" >> "${script}"
 echo "gunzip ${fastqFile}.gz" >> "${script}"
 echo "" >> "${script}"

  # Get fastQC analysis for iteration 2
 echo "# Get fastQC analysis for iteration 2" >> "${script}"
 echo "fastqc ${fastqFile} -o ${fastqcDirI2} -q" >> "${script}"
 echo "" >> "${script}"

 # Get the output file names
 bamfile="${bamdir}/${ffname}.bam"
 samfile="${samdir}/${ffname}.sam"
 unmapped_fastq="${unmappedFastqDir}/${ffname}.unmapped.fastq"

 # Mapping short reads to the reference genome using bowtie aligner
 # -q                 query input files are FASTQ .fq/.fastq (default)
 # -v <int>           report end-to-end hits w/ <=v mismatches; ignore qualities; No more than v mismatches in the entire length of the read
 # -p/--threads <int> number of alignment threads to launch (default: 1)
 # -S/--sam           write hits in SAM format
 # -k <int>           report up to <int> good alignments per read (default: 1)
 # -a/--all           report all alignments per read (much slower than low -k)
 # -m <int>           suppress all alignments if > <int> exist (def: no limit)
 #                    -m 3 instructs bowtie to refrain from reporting any alignments for reads having more than 3 reportable alignments. The -m option
 #                    is useful when the user would like to guarantee that reported alignments are "unique", for some definition of unique.
 echo "# Mapping short reads to the smncRNA reference genome using bowtie aligner" >> "${script}"
 echo "bowtie -q -v 1 -p 4 -S -a -m 1 --un ${unmapped_fastq} ${bowtieindexdir}/$(basename ${bowtieindexdir}) ${fastqFile} ${samfile}" >> "${script}"
 echo "" >> "${script}"

 # Move unmapped fastq file to unmapped fastq dir
 echo "# Zip unmapped_fastq files for next iteration" >> "${script}"
 echo "gzip ${unmapped_fastq}" >> "${script}"
 echo "" >> "${script}"

 # After mapping add the previously unmapped reads from iteration 1 to the unmapped reads from iteration2 that will be used to map on genome in iteration3
 echo "# After mapping add the previously unmapped reads from iteration 1 to the unmapped reads from iteration2 that will be used to map on genome in iteration3" >> "${script}"
 echo "cat ${unmapped_fastq}.gz ${unmappedFastqDirI1}/${ffname}.unmapped.fastq.gz > ${unmapped_fastq}_tmp.fastq.gz" >> "${script}" 
 echo "mv ${unmapped_fastq}_tmp.fastq.gz ${unmapped_fastq}.gz" >> "${script}" 
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
 echo "Submitting job for $f ..." 
 if [[ $queue =~ "short" ]]
 then 
  qtime="01:59"
 else
  qtime="48:00"
 fi
  bsub -w "done(${prevjobname})" -q ${queue} -n 4 -R "big" -W ${qtime} -e "${errorsFile}" -o "${stdoutFile}" -J ${jobname} sh "${script}"
 echo " "
done

