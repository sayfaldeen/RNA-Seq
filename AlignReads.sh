#!/bin/bash

#SBATCH --mem 200G
#SBATCH -c 50
#SBATCH -o global-run.1.out
# F/R orientation -- Seed length: 10

# Load the BT2 module
module purge
module load bowtie2


#Start time, in seconds since epoch
st=$(date +%s)


# Loop through all of the fastq files of interest
for f in $(ls /home/sayf/Projects/JewettCollab/CleanedSamples/3186_log_*_paired*fastq.gz)	
do

	# Generate a name for the resulting SAM files from the original fastq file names
	output_name=$(echo $f | awk -F ".1_paired" '{print $1"_paired.sam"}' | \
		awk -F "/" '{print "./"$NF}')
	
	# Get Read 2 filename
	r2=$(echo $f | sed 's/1_paired/2_paired/')

	# Create a filename for the reads that did not align concordantly
	un_out=$(echo $f | awk -F ".1_paired" '{print $1"_unaligned.fq.gz"}')

	# Run the BowTie2 aligner
	bowtie2 -p 100 -x /home/sayf/Projects/JewettCollab/BT2/Borrelia -1 ${f} -2 ${r2} -S ${output_name} \
		--fr -L 10 -i S,0,0.5 --end-to-end \
		--un-conc-gz ${un_out}

	#####----- Provided arguments to BowTie2 -----#####
	# -1: read 1 (forward) file name
	# -2: read 2 (reverse) file name
	# -S: resulting SAM file name (name to store output)
	# --fr: indicates that reads are given in Forward/Reverse orientation
	# -L: minumum seed length (10)
	# -i: interval between seeds; given as a formula. See BowTie2 manual
	# --end-to-end: specifies that alignment be conducted in end-to-end fashion
	# --un-conc-gz: specifies a filename for reads that did not align concordantly to be printed to

	echo ${output_name}


	# Compress SAM files in parallel
	echo """#!/bin/bash
#SBATCH --output=.bowzip.out
#SBATCH --qos short
pigz -p 2 -f -v ${output_name}""" > ./.bowzip
	sbatch ./.bowzip
done

#final time, in seconds since epoch
ft=$(date +%s)

#Elapsed time
et=$[ft-st]

#Elapsed minutes
ems=$[et/60]

echo "${j} took ${ems} minutes to complete"
