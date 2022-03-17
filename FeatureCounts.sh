#!/bin/bash

#SBATCH -c 64
#SBATCH -o updated-borrelia-features.slurm

#<Summary> Script to run featureCounts

featureCounts -T 64 -a ./updated-borrelia-annotation.saf -o ./updated-borrelia-features.txt -F 'SAF' \
	-O -p -d 20 ./*filtered.sam

#####----- Provided arguments to featureCounts -----#####
# T: number of cores (not threads apparently)
# a: annotation file
# o: Output path featureCounts matrix
# F: Format ('SAF'/'GTF')
# O: Option to count overlaps (if a read matches to more than one gene)
# p: Option to specify paired-end reads
# d: Min fragment (paired reads in featureCounts terminology) length (PE only option)
# final argument is the list of SAM/BAM files (using wildcard globbing
