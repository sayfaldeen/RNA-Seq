# README
- Folder containing scripts used for RNA-seq analysis

  - Scripts used for each section of the analysis are referenced and briefly explained

    - For details, read scripts; they are commented so as to be readily understandable

  - Scripts are numbered in parentheses to indicate order in work-flow 

    

## RNA-seq read-mapping pipeline

### Read trimming
- `./TrimReads.sh` (1)
	- Reads were trimmed using Trimmomatic
		- Min read length: 30
		- Sliding window of size 4 with minimum Phred score of 15
		- TruSeq PE adapters were clipped
	- No need for host-depletion since this was only from the organism of interest


### FastQ file quality control
- `./RunFqc.sh` (optional)
	- fqc and mqc programs were used to examine the quality of the sequencing files
	

### Read alignment
- `./AlignReads.sh` (2)
	- Cleaned (trimmed and clipped) paired-end reads were aligned using `BowTie2`
	- Only reads that aligned concordantly were stored for down-stream use

- `./CountExactConcs.py` (optional)
	- Counts and plots number of alignments that aligned exactly one time
		- Chosen to be conservative and not take multi-mapping reads; especially due to highly repetetive nature of organism of interest

- `./FilterSams.sh` (3)
	- Only reads that aligned concordantly **<u>EXACTLY ONCE</u>** were pulled out of the SAM files resulting from `BowTie2`
	  - Chosen to be conservative and not take multi-mapping reads, especially due to the large number of repeats within the genome
	
	    <br>

## Differential expression analysis

### Functional annotation
- `./GTF_to_SAF.py` (optional)
	- Script to turn a GTF file into an acceptable and simplified format for `featureCounts`
- `./FeatureCounts.sh` (4)
	- Using an annotation file (SAF) and alignment files (SAM) to assign reads to features (in our case, genes) using the `featureCounts` program
	
	  

### DESeq2 for differential expression
- `./RunDES2.R` (5)
	- `DESeq2` program was used to conduct differential expression analysis using the output of `featureCounts`
		- `ashr` method was used to conduct log fold-change shrinkage
			- Due to ease of comparing multiple different conditions within `DESeq2`
			- Shrinkage results are also very similar to 'apeglm' and 'normal' options for `lfcShrink()` function in `DESeq2`
			- Both `apeglm` and `ashr` exhibit low bias and high convergence with actual results when used for LFC shrinkage
				- `apeglm` was slightly more accurate with high LFC's (LFC > 4)

