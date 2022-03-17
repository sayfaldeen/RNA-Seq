#!/bin/bash

#<Script to create a BowTie2 index with given fa file and name> 
	# - $1: File
	# - $2: Name

# Load BowTie2
module load bowtie2


bowtie2-build ${1} ${2}
