


GJSrMap: The smallRNA mapping pipeline
====================================================

## Overview: 
An open source and fully customized pipeline that:
* Maps all small non-coding RNAs to customized and comprehensive reference sequences
* Is modularized and iterative
* Run on an HPCC cluster (default) or on a computer/server
* Performs quality control on the data
* Provides a detailed summary of mapping:
	* FastQC plots for every iteration
	* MultiQC plots for every iteration
	* Summary plots and statistics of smallRNA distribution and abundance
* Provides raw and normalized (RPKM) counts 
* Detailed logs of every iteration and steps

## Description
1. gjsrmap: SmallRNA mapping and analysis pipeline schematic:  

   ![fig_gjsrmap_overview](https://user-images.githubusercontent.com/10153240/50408167-63380900-07e6-11e9-8e93-f1716a5233f0.jpg)  

	- Pre-processing of the sequences:
		- This is the iteration 0 in the above schematic diagram
		- Preprocessing of the sequences to avoid multi mapping of the reads 
		- Build custom reference sequence indexes
		- Removal of low quality reads
		- Removal of 3' adapter sequences and size reduction of the reads

	- Iterative mapping of processed and filtered reads:
		- Iteration 1: Map reads between 16 to 33 bp to custom reference sequences of mature microRNAs and piRNAs  

		- Iteration 2: Map reads greater than 32 bp to custom reference sequences of other small non-coding RNAs. These are:  
		  ```
		  | Other small non-coding RNAs | Description                                          |
		  |-----------------------------|------------------------------------------------------|
		  | rRNA                        | Ribosomal RNA                                        |
		  | scRNA                       | Small cytoplasmic RNA                                |
		  | snRNA                       | Small nuclear RNA                                    |
		  | snoRNA                      | Small nucleolar RNA                                  |
		  | premiRNA                    | microRNA precursors                                  |
		  | osncRNA                     | Other small noncoding RNA                            |
		  | - tRNA                      | - Transfer RNA                                       |
		  | - Mt-tRNA                   | - Transfer RNA located   in the mitochondrial genome |
		  | - misc_RNA                  | - Miscellaneous other RNA                            |

		  ```
		- Iteration 3: Map the unmapped reads from iteration 1 and 2 to the species reference genome

	- Count the reads and distribute them to individual smallRNA classes

	- Generate QC, mapping and summary report

	  ![fig_qc](https://user-images.githubusercontent.com/10153240/50408146-20763100-07e6-11e9-9002-9f70682b8cff.jpg)

	  	- Sequence quality information
		- Bar plot of library sizes
		- Small non-coding RNA reads distrubution
		- Profile of expressed small non-coding RNAs (miRNAs in the above figure). Plots are also generated for other classes  as well

## Dependencies:
* ``bedtools``
* ``bowtie``
* ``bowtie2``
* ``cutadapt``
* ``fastqc``
* ``matplotlib``
* ``multiqc``
* ``numpy``
* ``samtools``
* ``scipy``

## Usage:

- Run wrapper with following options:
```
SPC=${1}                                            # Species: hsa or mmu or some other species
IFD=${2}                                            # Input Fastq Dir: input/fastq/test
ORD=${3}                                            # Output Results Dir: output/test
BWD=${4}                                            # path/to/bowtie/indexes
QUE=${5:-"fat"}                                     # mpi, fat, mpi-short, fat-short, mpi-long, fat-long
SPK=${6:-""}                                        # exiseq_spikein_dna_unique.fa or spike_rna1_unique.fa
threePadapter=${7:-"TGGAATTCTCGGGTGCCAAGG"}         # trueseq adapter
JID=${8:-"$(echo $HOME)/gjsrmap"}                   # Job dir
NCL=${9:-"input/annotation/rna_classes"}            # ncrna folder containing ncrna class fasta
```

- Example command:  
`` bash 06_run_ncRNA_mapping_usage.sh <above mentioned arguments>``

## Citation

- ```Gaurav Jain, GJSrMap: The smallRNA mapping pipeline, (2018), GitHub repository, https://github.com/gauravj49/gjsrmap```

- BibTeX entry:
```
@misc{gauravj49,
  author = {Jain, Gaurav},
  title = {GJSrMap: The smallRNA mapping pipeline},
  year = {2018},
  publisher = {GitHub},
  journal = {GitHub repository},
  howpublished = {\url{https://github.com/gauravj49/gjsrmap}},
  commit = {60c75be4ed1a60eebd656148f567ef23f251f3ac}
}
```
