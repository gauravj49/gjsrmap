#!/bin/bash

# Get input arguments
jobdir=$1         # "path_to_project_directory"
outputdir=$2      # "path_to_output_directory"
fastqdir=$3       # "path_to_fastq_files"
queue=${4:-"fat"} # mpi or fat

# Set up relevant filenames and directories
bamdir=${outputdir}/"bam"
samdir=${outputdir}/"sam"
trmfastqdir=${outputdir}/"unmapped_fastq"
unmfastqdir=${outputdir}/"trimmed_fastq"
mergedSamdir="${outputdir}/mergedSam"
scriptsdir=${jobdir}/scripts/05_cleanup/$(basename ${outputdir})

# Set up output folders
errorsFile="${outputdir}/mappingLogs/final_cleanup.err"
stdoutFile="${outputdir}/mappingLogs/final_cleanup.out"

# Create necessary dirs
mkdir -p ${scriptsdir}

# Initialize a script file
script="${scriptsdir}/$(basename ${outputdir}).sh"
touch ${script}
echo "#!/bin/bash" > "${script}"
echo "" >> "${script}"

# Zip the mappingLogs and fastqc folders
echo "# Zip the logs folder" >> "${script}"
#echo "tar -cJf ${outputdir}/mappingLogs.tar.xz ${outputdir}/mappingLogs" >> "${script}"
#echo "tar -cJf ${outputdir}/counts.tar.xz ${outputdir}/counts" >> "${script}"
echo "(cd ${outputdir} && tar -cJf mappingLogs.tar.xz mappingLogs)" >> "${script}"
echo "(cd ${outputdir} && tar -cJf counts.tar.xz counts)" >> "${script}"
echo "" >> "${script}"

# Cleanup some mappingLogs file
echo "# Remove intermediate folders and files" >> "${script}"
echo "rm -rf ${mergedSamdir} ${bamdir} ${samdir} ${trmfastqdir} ${unmfastqdir} ${outputdir}/mappingLogs" >> "${script}"
echo "" >> "${script}"

# Get the jobname to submit each job as a job array
projName=$(basename ${outputdir})
jobname="5_cleanup_${projName}"
prevjobname="4_summaryPlots_${projName}*"

# Submit the job
chmod 775 "${script}"
echo "Submitting job for ${projName} ..." 
bash "${script}"

echo " "

