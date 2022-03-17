#<Summary> Script to use DESeq2 for differential expression analysis of given count data from featureCounts

# $1: featureCounts output

#####----- Preliminary setup -----#####
# Set the correct working directory
setwd("/home/sayf/Projects/JewettCollab/DE/")

# Pull the current date and time
datetime = format(Sys.time(), "%A, %b %d, %Y @ %X")

# Import the necessary libraries
library("ggplot2") # For making plots
library("DESeq2") # For differential analysis

#####----- Import the data -----#####
ImportData <- function() {

	# Import the featureCounts output
	fc_out = read.csv("updated-borrelia-features.txt", sep="\t", skip=1, row.names="Geneid")
	# Line 1 is skipped because it only contains commands given to fC, not any data
	
	# Extract only the columns for the gene counts
	cts = as.matrix(fc_out[,grepl("sam", colnames(fc_out))])
	
	# Change the colnames to only reflect sample names rather than filenames; more concise
	colnames(cts) <- gsub(pattern="^..", replacement="", x=gsub(pattern="_paired_filtered.sam", replacement="", x=colnames(cts)))
	
	# Create a varibale to house the sample 'treatment' (condition) data
	coldata <- read.csv("sample-condition-data.csv", row.names=1)["condition"] # Samplenames (index) and condition column
	
	# re-arrange the colnames to match the rowname of coldata
	cts <- cts[, rownames(coldata)]
	
	return(list("cts" = cts, "coldata" = coldata))
	
}

RunDES <- function(cts, coldata) {
	#####----- Run the Differential expression analysis -----##### 
	
	MakeDDS <- function(cts, coldata) {
	# Create a DESeq2 data object from the count-matrix (cts)
	dds <- DESeqDataSetFromMatrix(countData = cts, colData = coldata, design = ~ condition) # 'treatment' column in coldata
	
	# Create the metadata variable for the data
	featureData <- data.frame(gene=rownames(cts)) # Must be in data.frame format
	mcols(dds) <- DataFrame(mcols(dds), featureData)
	mcols(dds)
	
	# Pre-filter; DESeq2 recommends removing genes with less than 10 total reads
	# - Mainly to reduce memory and computational cost with large datasets
	# - I will only remove genes with no reads; dataset not large and we are running independent filtering anyways
	# - Also makes the summary print-out more informative (gives total number of genes, not just the ones after the filtering)
	keep <- rowSums(counts(dds)) > 0
	dds <- dds[keep,]
	
	# Runs estimation of size factors, estimation of dispersion: neg. binomial GLM
	dds <- DESeq(dds)
	return(dds)
	}
	
	
	#####----- Run the log-fold change shrinkage (LFC) -----#####
	RunDEA <- function(dds) {
	library("ashr") # Adaptive shrinkage; utilizes the FSR (False Sign Rate); should be more conservative
	library("xlsx")
		# Also easier to use when attempting to contrast multiple conditions rather than setting one as a reference
		# Very comparable in performance to apeglm, does potentially introduce slightly more bias while shrinking large LFC's for reads with 
	ref_mut_stat_resLFC = lfcShrink(dds, contrast=c("condition", "refStat", "mutantStat"), type="ashr") # Conduct shrinkage WT (reference) to mutant for stationary phase
	ref_mut_stat_sig <- subset(na.omit(ref_mut_stat_resLFC[order(ref_mut_stat_resLFC$padj),]), padj < 0.05) # Extract only differences with an adjusted p-val < 0.05 (BH FDR)
	write.csv(ref_mut_stat_resLFC, "wt_v_mutant_stat_tot.csv")	
	write.xlsx(ref_mut_stat_sig, file="DESeq2-results.xlsx", sheetName="wt_v_mutant_stat_sigs") # Then re-order by p-value
	
	ref_mut_log_resLFC = lfcShrink(dds, contrast=c("condition", "refLog", "mutantLog"), type="ashr") # WT v genetic complement
	ref_mut_log_sig <- subset(na.omit(ref_mut_log_resLFC[order(ref_mut_log_resLFC$padj),]), padj < 0.05) # Extract only differences with an adjusted p-val < 0.05 (BH FDR)
	write.csv(ref_mut_log_resLFC, "wt_v_mutant_log_tot.csv")
	write.xlsx(ref_mut_log_sig, file="DESeq2-results.xlsx", sheetName="wt_v_mutant_log_sigs", append=T)  # Then re-order by p-value
	
	ref_comp_stat_resLFC <- lfcShrink(dds, contrast=c("condition", "refStat", "compStat"), type="ashr") # Compare WT (reference) to genetic complement for stationary phase
	ref_comp_stat_sig <- subset(na.omit(ref_comp_stat_resLFC[order(ref_comp_stat_resLFC$padj),]), padj < 0.05) # Extract only differences with an adjusted p-val < 0.05 (BH FDR)
	write.csv(ref_comp_stat_resLFC, "wt_v_complement_stat_tot.csv")	
	write.xlsx(ref_comp_stat_sig, file="DESeq2-results.xlsx", sheetName="wt_v_complement_stat_sigs", append=T)
	
	ref_comp_log_resLFC <- lfcShrink(dds, contrast=c("condition", "refLog", "compLog"), type="ashr") # Compare WT (reference) to genetic complement for log phase
	ref_comp_log_sig <- subset(na.omit(ref_comp_log_resLFC[order(ref_comp_log_resLFC$padj),]), padj < 0.05) # Extract only differences with an adjusted p-val < 0.05 (BH FDR)
	write.csv(ref_comp_log_resLFC, "wt_v_complement_log_tot.csv")
	write.xlsx(ref_comp_log_sig, file="DESeq2-results.xlsx", sheetName="wt_v_complement_log_sigs", append=T)
	
	mut_comp_stat_resLFC <- lfcShrink(dds, contrast=c("condition", "mutantStat", "compStat"), type="ashr") # Compare mutant to genetic complement for stationary phase
	mut_comp_stat_sig <- subset(na.omit(mut_comp_stat_resLFC[order(mut_comp_stat_resLFC$padj),]), padj < 0.05) # Extract only differences with an adjusted p-val < 0.05 (BH FDR)
	write.csv(mut_comp_stat_resLFC, "mutant_v_complement_stat_tot.csv")
	write.xlsx(mut_comp_stat_sig, file="DESeq2-results.xlsx", sheetName="mutant_v_complement_stat_sigs", append=T)
	
	mut_comp_log_resLFC <- lfcShrink(dds, contrast=c("condition", "mutantLog", "compLog"), type="ashr") # Compare mutant to genetic complement for log phase
	mut_comp_log_sig <- subset(na.omit(mut_comp_log_resLFC[order(mut_comp_log_resLFC$padj),]), padj < 0.05) # Extract only differences with an adjusted p-val < 0.05 (BH FDR)
	write.csv(mut_comp_log_resLFC, "mutant_v_complement_log_tot.csv")
	write.xlsx(mut_comp_log_sig, file="DESeq2-results.xlsx", sheetName="mutant_v_complement_log_sigs", append=T)
	
	# Print a summary of the results after shirnkage
		# Outliers: calculated using Cook's distance
			# - A high distance for a specific point indicates that th point is very influential and may be an outlier
		# Low counts: Number of samples removed during independent filtering
			# - Uses dynamic CPM threshold (counts/million reads)
			# - Finds the CPM that provides the highest number of genes with FDR > 0.05

	cat("Analysis conducted on: ", datetime, "\n")

	cat("\nDiffential expression analysis summary: WT vs mutant (stationary phase)\n")
	cat("=======================================================================")
	cat(summary(ref_mut_stat_resLFC, alpha=0.05))
	
	cat("\n\nDiffential expression analysis summary: WT vs mutant (log phase)\n")
	cat("================================================================")
	cat(summary(ref_mut_log_resLFC, alpha=0.05))
	
	
	cat("\n\nDiffential expression analysis summary: WT vs genetic complement (stationary phase)\n")
	cat("===================================================================================")
	cat(summary(ref_comp_stat_resLFC, alpha=0.05))

	cat("\n\nDiffential expression analysis summary: WT vs genetic complement (log phase)\n")
	cat("============================================================================")
	cat(summary(ref_comp_log_resLFC, alpha=0.05))
	
	
	cat("\n\nDiffential expression analysis summary: mutant vs genetic complement (stationary phase)\n")
	cat("=======================================================================================")
	cat(summary(mut_comp_stat_resLFC, alpha=0.05))
	
	cat("\n\nDiffential expression analysis summary: mutant vs genetic complement (log phase)\n")
	cat("================================================================================")
	cat(summary(mut_comp_log_resLFC, alpha=0.05))
	
	# Return the Shrunken LFC data
	Results <- list("dds" = dds,
					"ref_mut_stat_resLFC" = ref_mut_stat_resLFC,
					"ref_mut_log_resLFC" = ref_mut_stat_resLFC,
					"ref_comp_stat_resLFC" = ref_comp_stat_resLFC,
					"ref_comp_log_resLFC" = ref_comp_stat_resLFC,
					"mut_comp_stat_resLFC" = mut_comp_stat_resLFC,
					"mut_comp_log_resLFC" = mut_comp_stat_resLFC)

	return(Results)
	}
	
	dds = MakeDDS(cts, coldata)
	return(RunDEA(dds))
	
}

PlotDiffs <- function(resLFC, out="null-log-fold-diffs.tiff") {
	#####----- Create the 'MA'-plot; really a log fold-change scatter plot between two conditions -----#####
	# find the max and min log-fold changes and set those as the axes
	LFC_min <- min(resLFC$log2FoldChange)
	LFC_max <- max(resLFC$log2FoldChange) 
	
	# Make the MA-plot
	plotMA(resLFC, ylim=c(LFC_min, LFC_max))
	
	# Save the image
	ggsave(out, device="tiff", width=10, height=6, units="cm", dpi=300)
	
}

PlotPCA <- function(dds, out="null-pca-plot.tiff") {
	#####----- Make the PCA plot -----#####
	
	# Transform the data
	vsd <- vst(dds, blind=T) # Use variance-stabilizing transformation
	
	# Draw the PCA plot
	plotPCA(vsd, intgroup="condition")
	
	ggsave(out, device="tiff", width=10, height=6, units="cm", dpi=300)
}

# Run the functions for analysis
Data = ImportData()
Results = RunDES(cts=Data$cts, coldata=Data$coldata)


#####----- Draw and save the plots -----#####

# Plot and save the PCA results
PlotPCA(Results$dds, out="./Borrelia-riboflavin-transport-experiment-pca.tiff")

# Plot the differential abundance results
PlotDiffs(Results$ref_mut_stat_resLFC, out="wt_v_mutant_stat_diffs.tiff")
PlotDiffs(Results$ref_mut_log_resLFC, out="wt_v_mutant_log_diffs.tiff")

PlotDiffs(Results$ref_comp_stat_resLFC, out="wt_v_complement_stat_diffs.tiff")
PlotDiffs(Results$ref_comp_log_resLFC, out="wt_v_complement_log_diffs.tiff")

PlotDiffs(Results$mut_comp_stat_resLFC, out="mut_v_complement_stat_diffs.tiff")
PlotDiffs(Results$mut_comp_log_resLFC, out="mut_v_complement_log_diffs.tiff")
