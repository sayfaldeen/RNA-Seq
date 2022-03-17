#!/bin/bash

#SBATCH -c 64

#<># Script to trim the reads based on q-score

# Pull the forward reads
for inp1 in $(ls ../SampleFiles/*_*1.fastq.gz)
do
	# Get the sample names for the reverse reads
	inp2=$(echo ${inp1} | sed 's/1.fastq/2.fastq/g')
	
	# Create the output file names
	out1_pair=$(echo ${inp1} | awk -F "/" '{print $3}' | sed 's/1.fastq/1_paired.fastq/')
	out1_unpair=$(echo ${inp1} | awk -F "/" '{print $3}' | sed 's/1.fastq/1_unpaired.fastq/')

	out2_pair=$(echo ${inp2} | awk -F "/" '{print $3}' | sed 's/2.fastq/2_paired.fastq/')
	out2_unpair=$(echo ${inp2} | awk -F "/" '{print $3}' | sed 's/2.fastq/2_unpaired.fastq/')
	

	# Run Trimmomatic
	java -jar /home/sayf/bin/trimmomatic-0.36.jar PE -phred33 \
		${inp1} ${inp2} \
		${out1_pair} ${out1_unpair} \
		${out2_pair} ${out2_unpair} \
		-threads 128 MINLEN:30 SLIDINGWINDOW:4:15 \
		ILLUMINACLIP:/home/syooseph/utils/Trimmomatic-0.36/adapters/TruSeq3-PE-2.fa:2:30:10
done
