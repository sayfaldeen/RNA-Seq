#!/bin/bash

# <Summary> Script to pull out reads that aligned concordantly only once from the SAM files

#SBATCH -c 12

# $1: SAM file
# $2: Output name
# $3: Number of threads

function FilterSams {
	# Decompress sam files using multi-threading -> sort sam files by name -> \
	# pull only samples that have 'YT:Z:CP' flag (concordantly paired) -> \
	# pull samples without 'XS' flag ('XS' indicates multiple valid mappings
	pigz -d -p ${3} -c ${1} | samtools sort -n -O sam -@ ${3} |  grep "YT:Z:CP" | grep -v "XS" > ${2}
}


# Iterate through the files and filter out the reads that align concordantly more than once
for f in ./*paired.sam.gz
do
	o=$(echo $f | sed 's/paired/paired_filtered/g')
	FilterSams ${f} ${o} 24

	echo "${f} filtered into ${o}"
done
