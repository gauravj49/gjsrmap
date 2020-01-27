#!/bin/bash

# Get function arguments
SPC=${1}                                            # Species: hsa or mmu or some other species
IFD=${2}                                            # Input Fastq Dir: input/fastq/test
ORD=${3}                                            # Output Results Dir: output/test
threePadapter=${4:-"TGGAATTCTCGGGTGCCAAGG"}         # trueseq adapter
SPK=${5:-""}                                        # exiseq_spikein_dna_unique.fa or spike_rna1_unique.fa
JID=${6:-"/home/rad/users/gaurav/projects/gjsrmap"} # Job dir
NCL=${7:-"input/annotation/rna_classes"}            # ncrna folder containing ncrna class fasta
BWD=${8:-"/home/rad/packages/bowtie/indexes"}                                        # path/to/bowtie/indexes
ifTRIM=1                                            # If you want to trim the reads (default is yes)

# Get mapping related input files
if [ ${SPC} = "hsa" ]; then
 gSPC="hg19"
fi

if [ ${SPC} = "mmu" ]; then
 gSPC="mm10"
fi

if [ ${SPC} = "macFas5" ]; then
 gSPC="macFas5"
fi


bowtiepimirnaIDX="${BWD}/${SPC}_pimirna"
bowtiesmncrnaIDX="${BWD}/${SPC}_smncrna"
bowtiegenomeIDX="${BWD}/${gSPC}"
pimiRNAfasta="input/reference_genome/dna/${SPC}_dna_pimirna_unique.fa"
smncRNAfasta="input/reference_genome/dna/${SPC}_dna_smncrna_unique.fa"
NCRNAfasta="input/reference_genome/dna/${SPC}_dna_allsmncrna_unique.fa"

# If spike in sequence
if [ "${SPK}" != "" ]
then
 spk_name=$(basename $(basename ${SPK} _unique.fa) _dna)
 pimiRNAfasta="input/reference_genome/dna/${SPC}_dna_pimirna_${SPK}"
 bowtiepimirnaIDX="${BWD}/${SPC}_pimirna_${spk_name}"
 NCRNAfasta="input/reference_genome/dna/${SPC}_dna_allsmncrna_${SPK}"

 echo "############## Using Spike in RNA ###########"
 echo ""
 echo ${NCRNAfasta}
 echo ${bowtiepimirnaIDX}
 echo ${pimiRNAfasta}
 echo ""
 echo "##############  Submitting jobs  ############"
 echo ""
fi

echo "# Map to pimirna"
echo "sh 02.1_map_pimirna_bowtie_usage.sh ${JID} ${bowtiepimirnaIDX} ${IFD} ${ORD} ${pimiRNAfasta} ${ifTRIM} ${threePadapter}"
echo ""

echo "# Map smncrna (small ncrna: osncrna, rrna, snrna, snorna, premirna)"
echo "sh 02.2_map_smncrna_bowtie_usage.sh ${JID} ${bowtiesmncrnaIDX} ${IFD} ${ORD} ${smncRNAfasta} ${ifTRIM} ${threePadapter}"
echo ""

echo "# Map remaining unmapped reads to genome"
echo "sh 03_map_reads_to_genome_bowtie_usage.sh ${JID} ${bowtiegenomeIDX} ${ORD} ${NCRNAfasta} ${SPC} ${IFD} ${NCL}"
echo ""

echo "# Get balanced counts and summary for counts and balanced counts"
echo "sh 04_counts_summary_usage.sh ${JID} ${ORD} ${IFD} ${NCRNAfasta} ${SPC} ${NCL}"
echo ""

echo "# Remove unnecessary folders"
echo "sh 05_cleanup_usage.sh ${JID} ${ORD} ${IFD}"
echo ""

# # Map to pimirna
# sh 02.1_map_pimirna_bowtie_usage.sh ${JID} ${bowtiepimirnaIDX} ${IFD} ${ORD} ${pimiRNAfasta} ${ifTRIM} ${threePadapter}

# # Map smncrna (small ncrna: osncrna, rrna, snrna, snorna, premirna)
# sh 02.2_map_smncrna_bowtie_usage.sh ${JID} ${bowtiesmncrnaIDX} ${IFD} ${ORD} ${smncRNAfasta} ${ifTRIM} ${threePadapter}

# # Map remaining unmapped reads to genome
# sh 03_map_reads_to_genome_bowtie_usage.sh ${JID} ${bowtiegenomeIDX} ${ORD} ${NCRNAfasta} ${SPC} ${IFD} ${NCL}

# # Get balanced counts and summary for counts and balanced counts
# sh 04_counts_summary_usage.sh ${JID} ${ORD} ${IFD} ${NCRNAfasta} ${SPC} ${NCL}

# # Remove unnecessary folders
# sh 05_cleanup_usage.sh ${JID} ${ORD} ${IFD}
