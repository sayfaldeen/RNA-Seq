#!/bin/bash

#SBATCH -c 12

#<># Script to run FastQC and combine the data with MQC

# Run FastQC
fqc ./*_paired.fastq.gz -t 24 --outdir ./

# Run MQC to combine all FQC's into one file
/home/sayf/.local/bin/multiqc ./ -i ./ -b "After quality control and adapter trimming" -n PostQC
