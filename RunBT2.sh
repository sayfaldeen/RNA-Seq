#!/bin/bash

#SBATCH --mem 200G
#SBATCH -c 50
#SBATCH -o Borrelia-BT2-alignment.slurm
# F/R orientation -- Seed length: 10

# Load the BT2 module
module purge
module load bowtie2


#Start time, in seconds since epoch
st=$(date +%s)


for f in $(ls /home/sayf/PhD/Projects/JewettCollab/CleanedSamples/*.1_paired.fastq.gz)
do

	output_name=$(echo $f | awk -F ".1_paired.fastq.gz" '{print $1"paired.sam"}' | \
		awk -F "/" '{print "./"$NF}')
	
	# Get Read 2 filename
	r2=$(echo $f | sed 's/1_paired/2_paired/')

	# Create a filename for the reads that did not align concordantly
	un_out=$(echo $f | awk -F "1_paired.fastq.gz" '{print $1"un-concordant.fq.gz"}' | \
		awk -F "/" '{print "./"$NF}')

	# Remember to give paired reads
	bowtie2 -p 100 -x /home/sayf/AzarianLab/JewettCollab/BT2/Borrelia -1 ${f} -2 ${r2} -S ${output_name} \
		--fr -L 10 --end-to-end --very-fast \
		--un-conc-gz ${un_out}
		# --fr: Forward/Reverse orientation of R1/R2, respectively
		# -L 10: Seed length of 10
		# --end-to-end: Request end to end alignment
		# --very-fast: Less sensitive alignment but faster

	echo ${output_name}


<<troubleshoot
	# Make sure the names are all right
	if [[ -z $r2 ]]
	then
		echo "R2 does not exist"
	else
		echo "R2 exists"
	fi


	#echo "${f} -- ${r2}"
	#echo "${output_name} -- ${un_out}"
troubleshoot

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

echo "${j} took ~ ${ems} minutes to complete"


# Clean up the files

if [[ ! -d UnConcs ]]
then
	mkdir UnConcs
fi

mv *un-concordant* UnConcs
