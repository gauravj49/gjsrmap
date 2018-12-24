#!/bin/bash

countsDir=$1
species=$2
rnaClassDir=$3


rawCountsDir=${countsDir}/raw_counts 
normCountsDir=${countsDir}/normalized_counts

mkdir -p ${rawCountsDir}/allncrna ${normCountsDir}/allncrna
for r in rRNA snRNA snoRNA osncRNA mirna pirna premirna
do
 # Raw counts
 mkdir -p ${rawCountsDir}/${r}
 for f in ${countsDir}/*_counts.txt
 do
  of=${rawCountsDir}/${r}/$(basename ${f} _counts.txt)_${r}Counts.txt
  awk -F '\t' 'NR==FNR {id[$1]; next} $1 in id' ${rnaClassDir}/${species}_${r}.txt ${f} > ${of}
 done;

 # Normalized counts
 mkdir -p ${normCountsDir}/${r}
 for f in ${countsDir}/*_counts_normalized.txt
 do
  of=${normCountsDir}/${r}/$(basename ${f} _counts_normalized.txt)_${r}Counts_normalized.txt
  awk -F '\t' 'NR==FNR {id[$1]; next} $1 in id' ${rnaClassDir}/${species}_${r}.txt ${f} > ${of}
 done;
done;

for f in ${countsDir}/*_counts.txt
do
 aof=${rawCountsDir}/allncrna/$(basename ${f} _counts.txt)_allncrnaCounts.txt
 mv ${f} ${aof}
done
for f in ${countsDir}/*_counts_normalized.txt
do
 aof=${normCountsDir}/allncrna/$(basename ${f} _counts_normalized.txt)_allncrnaCounts_normalized.txt
 mv ${f} ${aof}
done;



