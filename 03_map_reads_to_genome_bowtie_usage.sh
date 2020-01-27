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
outputdir=$3      # "path_to_output_directory"
reffastafile=$4   # "path_to_reference_genome/dna/hsa_dna_mirna_unique.fa"
species=$5        # "hsa" or "mmu"
fastqdir=$6       # "path_to_fastq_files"
rnaClassDir=${7:-"input/annotation/rna_classes"}   # input/annotation/rna_classes
queue=${8:-"fat"} # mpi or fat

# Set up relevant filenames and directories
bamdir="${outputdir}/bam/iteration3"
samdir="${outputdir}/sam/iteration3"
finalbamdir="${outputdir}/bamFiles"
unmappedFastqDir="${outputdir}/unmapped_fastq/iteration2"
countsdir="${outputdir}/counts"
mappingLogsDir="${outputdir}/mappingLogs/iteration3"
fastqcDirI3="${outputdir}/fastqc/iteration3"
fastqcDirIm="${outputdir}/fastqc/iteration{1,2,3}"
mergedSamdir="${outputdir}/mergedSam"
scriptsdir=${jobdir}/scripts/03_map_genome_bowtie/$(basename ${outputdir})

mkdir -p ${finalbamdir} ${samdir} ${bamdir} ${countsdir} ${scriptsdir} ${trimmedFastqDir} ${unmappedFastqDir} ${fastqcDirI3} ${mergedSamdir}

# Process and submit jobs
for f in ${fastqdir}/*.fastq.gz
do
 # Get the fastq file names
 ffname=$(basename ${f} .fastq.gz)

 # Get unmapped fastq file
 fastqFile="${unmappedFastqDir}/${ffname}.unmapped.fastq"

 # Get the jobname to submit each job as a job array
 jobname="3_$(basename ${outputdir})_${ffname}"
 prevjobname="2_$(basename ${outputdir})_${ffname}"

 # Set up output folders
 bowtieoutputdir=${mappingLogsDir}/${ffname}
 errorsFile="${bowtieoutputdir}/${ffname}_gwdg_mapping.err"
 stdoutFile="${bowtieoutputdir}/${ffname}_gwdg_mapping.out"
 
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

 # Get fastQC analysis for iteration 3
 echo "# Get fastQC analysis for iteration3" >> "${script}"
 echo "fastqc ${fastqFile} -o ${fastqcDirI3} -q" >> "${script}"
 echo "" >> "${script}"

 #  Mapping all unmapped short reads to the reference genome using bowtie aligner
 samfile="${samdir}/${ffname}.sam"
 echo "# Mapping all unmapped short reads to the reference genome using bowtie aligner" >> "${script}"
 echo "bowtie -q -v 1 -p 4 -S -a -m 1 ${bowtieindexdir}/$(basename ${bowtieindexdir}) ${fastqFile} ${samfile}" >> "${script}"
 echo "" >> "${script}"

 # Get the BAM files from the SAM files for iteration3
 bamfile3="${bamdir}/${ffname}.bam"
 echo "# Get the BAM files from the SAM files for iteration3" >> "${script}"
 echo "samtools view -Sb ${samfile} > ${bamfile3}" >> "${script}"
 echo "" >> "${script}"

 # Sort the bam file before creating the index
 echo "# Sort the bam file before creating the index" >> "${script}"
 echo "samtools sort -o ${bamfile3}.sorted ${bamfile3} " >> "${script}"
 echo "mv ${bamfile3}.sorted ${bamfile3}" >> "${script}"
 echo "" >> "${script}"

 # Create the index of the bam file
 echo "# Create the index of the bam file" >> "${script}"
 echo "samtools index ${bamfile3}" >> "${script}"
 echo "" >> "${script}"

 # Merge the bam files
 # -n The input alignments are sorted by read names rather than by chromosomal coordinates
 # -r Attach an RG tag to each alignment. The tag value is inferred from file names
 # -O BAM or SAM
 echo "# Get the merged bam file from all three iterations" >> "${script}"
 bamfile1="${outputdir}/bam/iteration1/${ffname}.bam"
 bamfile2="${outputdir}/bam/iteration2/${ffname}.bam"
 mergedBam="${finalbamdir}/${ffname}.bam"
 echo "samtools merge -n -r -O BAM ${mergedBam} ${bamfile1} ${bamfile2} ${bamfile3}" >> "${script}"
 echo "" >> "${script}"

 # Sort the final merged (mapped + unmapped) bam
 echo "# Sort the final merged (mapped + unmapped) bam file before getting counts" >> "${script}"
 echo "samtools sort -o ${mergedBam}.sorted ${mergedBam} " >> "${script}"
 echo "mv ${mergedBam}.sorted ${mergedBam}" >> "${script}"
 echo "" >> "${script}"

 # Create the index of the final merged (mapped + unmapped) bam file
 echo "# Create the index of the final merged (mapped + unmapped) bam file" >> "${script}"
 echo "samtools index ${mergedBam}" >> "${script}"
 echo "" >> "${script}"

 # Get high quality (mapq >= 30) mapped sam files for counting
 mergedSam=${mergedSamdir}/${ffname}.sam
 echo "# Get high quality (mapq >= 30) mapped sam files for counting" >> "${script}"
 echo "samtools view -F 4 -q 30 ${mergedBam} > ${mergedSam}" >> "${script}"
 echo "" >> "${script}"

 # Get the counts from SAM file
 countsfile="${countsdir}/${ffname}_counts.txt" 
 echo "# Get the counts from SAM file" >> "${script}"
 echo "python scripts/get_counts_from_sam.py -if=${reffastafile} -is=${mergedSam} -of=${countsfile}" >> "${script}"
 echo "" >> "${script}"

 chmod 775 "${script}"
 echo "Submitting job for $jobname ..." 
 sh "${script}"
 echo " "
done

