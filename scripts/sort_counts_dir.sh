#!/bin/bash

# Get function arguments
countsDIR=$1        # counts_dir: 

# sort the count files by first column
mkdir -p ${countsDIR}/unsorted
mv ${countsDIR}/*.txt ${countsDIR}/unsorted
for f in ${countsDIR}/unsorted/*.txt
do
 of="${countsDIR}/$(basename ${f})"
 # skip the footer which starts with "__"
 grep -vE "__no_feature|__ambiguous|__too_low_Qual|__not_aligned|__alignment_not_unique" ${f} | sort -k1,1 -u > ${of}
 grep -E  "__no_feature|__ambiguous|__too_low_Qual|__not_aligned|__alignment_not_unique" ${f} >> ${of}
done
rm -rf ${countsDIR}/unsorted

# check if all the files has same number of mirs
find ${countsDIR} -name *_*.txt -type f  -exec wc -l {} \;
mkdir -p ${countsDIR}/test_counts

for f in ${countsDIR}/*.txt
do
 b=$(basename ${f} .txt)
 cut -f1 ${f} > ${countsDIR}/test_counts/${b}.txt
done
a1=`ls ${countsDIR}/test_counts | head -1`

for f in ${countsDIR}/test_counts/*.txt
do
 diff ${countsDIR}/test_counts/${a1} $f
done
rm -rf ${countsDIR}/test_counts

