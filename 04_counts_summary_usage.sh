#!/bin/bash

# # Load necessary modules
# python="/usr/users/gjain/bin/libs/python/bin/python"

# Get input arguments
jobdir=$1         # /usr/users/gjain/bin/projects/map2ncrna
outputdir=$2      # output/csf/
fastqdir=$3       # input/fastq/merged_untrimmed/test
reffastafile=$4   # input/reference_genome/dna/hsa_dna_allsmncrna_unique.fa
species=$5        # hsa or mmu
rnaClassDir=${6}  # input/annotation/rna_classes
queue=${7:-"fat"} # mpi or fat

# Set up relevant filenames and directories
projName=$(basename ${outputdir})
bamdir=${outputdir}/"bamFiles"
countsDir=${outputdir}/"counts"
summaryPlots=${outputdir}/"summary_plots/all_samples"
scriptsdir=${jobdir}/scripts/04_counts_summary_balance/${projName}
trimmedFastqDir=${outputdir}/"trimmed_fastq/iteration1"
fastqcDirI0="${outputdir}/fastqc/iteration0"
fastqcDirI1="${outputdir}/fastqc/iteration1"
fastqcDirI2="${outputdir}/fastqc/iteration2"
fastqcDirI3="${outputdir}/fastqc/iteration3"
fastqcDirIm="${outputdir}/fastqc/iteration{1,2,3}"
multiqcDirI0="${outputdir}/multiqc/iteration0"
multiqcDirI1="${outputdir}/multiqc/iteration1"
multiqcDirI2="${outputdir}/multiqc/iteration2"
multiqcDirI3="${outputdir}/multiqc/iteration3"
multiqcDirIm="${outputdir}/multiqc/allMapped"

# Set up output folders
errorsFile="${outputdir}/mappingLogs/04_counts_summary.err"
stdoutFile="${outputdir}/mappingLogs/04_counts_summary.out"

# Create necessary dirs
mkdir -p ${scriptsdir} ${summaryPlots}

# Initialize a script file
script="${scriptsdir}/${projName}.sh"
touch ${script}
echo "#!/bin/bash" > "${script}"
echo "" >> "${script}"

 # Divide counts into classes
 echo "# Divide counts into classes" >> "${script}"
 echo "sh scripts/distribute_counts_usage.sh ${countsDir} ${species} ${rnaClassDir}" >> "${script}"
 echo "" >> "${script}"

# Get MultiQC analysis before trimming
echo "# Get MultiQC analysis before any trimming" >> "${script}"
echo "multiqc -o ${multiqcDirI0} -n iteration0 ${fastqcDirI0}" >> "${script}"
echo "" >> "${script}"

# Get MultiQC analysis after iteration 1
echo "# Get MultiQC analysis after trimming" >> "${script}"
echo "multiqc -o ${multiqcDirI1} -n iteration1 ${fastqcDirI1}" >> "${script}"
echo "" >> "${script}"

# Get MultiQC analysis for iteration 2
echo "# Get MultiQC analysis for iteration 2" >> "${script}"
echo "multiqc -o ${multiqcDirI2} -n iteration2 ${fastqcDirI2}" >> "${script}"
echo "" >> "${script}"

# Get MultiQC analysis for iteration 3
echo "# Get MultiQC analysis for iteration 3" >> "${script}"
echo "multiqc -o ${multiqcDirI3} -n iteration3 ${fastqcDirI3}" >> "${script}"
echo "" >> "${script}"

# Get MultiQC analysis after mapping (iteration 1,2,3)
echo "# Get MultiQC analysis after mapping (iteration 1,2,3)" >> "${script}"
echo "multiqc -o ${multiqcDirIm} ${fastqcDirIm}" >> "${script}"
echo "" >> "${script}"

# Get summary plot for all samples
echo "" >> "${script}"
echo "# Get summary plot for all samples" >> "${script}"
echo "${python} scripts/ncRNA_mapping_summary.py -id=${countsDir} -of=${summaryPlots}/${projName}_all_samples.txt " >> "${script}"
echo "" >> "${script}"

# Get the jobname to submit each job as a job array
jobname="4_summaryPlots_${projName}"
prevjobname="3_${projName}_*"

# Submit the job
chmod 775 "${script}"
echo "Submitting job for ${projName} ..." 
bash "${script}"

echo " "

