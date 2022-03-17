#!/usr/bin/python3

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os
import subprocess as sp

# Create the argument parser
import argparse as ap
parser = ap.ArgumentParser(description="Script to convert GTF to SAF.")

# Add the arguments
parser.add_argument("-g", "--gtf", dest="gtf_file_name",
        help="path/to/gtf")

parser.add_argument("-o", "--output", dest="out",
        help="Specify name of output file")

args = parser.parse_args()

# Ensure that an input file was given
if not args.gtf_file_name or args.gtf_file_name.split(".")[-1].lower() != "gtf":
    parser.error("Please specify a GTF file. See --help for options.")

# If not output name was given, automatically create one
if not args.out:
    args.out = args.gtf_file_name.replace("gtf", "saf")
    print(f"No output file name was given. '{args.out}' automatically generate as output file name")

# Create the function to convert GTF to SAF
def GTF_to_SAF(gtf_file_name):
    """
    Function to convert GTF to SAF

    Parameters
    ----------
    gtf: GTF file as DF object

    Returns
    -------
    SAF file as DF object

    """
    # Clean the GTF file
    gtf_file_name = "../AssemblyFiles/GCF_000008685.2_ASM868v2_genomic.gtf"
    gtf_c = gtf_file_name.replace(".gtf", "_cleaned.gtf")
    sp.run(f"grep -v '#' {gtf_name} > {gtf_c}", shell=True)
    gtf = pd.read_csv(gtf_c, sep="\t", header=None)
    gene_gtf = gtf[gtf[2] == "gene"] # Only focus on reads labaled as genes

    # Create the SAF
    saf = gene_gtf[[0, 3, 4, 6, 8]].copy()
    saf.columns = ["Chr", "Start", "End", "Strand", "Attribute"]

    # Pull the GeneID's
    GeneIDs = saf.Attribute.fillna("none").apply(lambda x:str(re.findall('gene_id "\w{10}"', x)).split('"')[1] \
                           if x != "none" else np.nan)

    #Add the GeneID's to the SAF
    saf["GeneID"] = GeneIDs

    # Make the GeneID's the index column to conform to SAF format
    saf.set_index("GeneID", inplace=True)

    # Remove 'Attribute' column to conform to SAF format
    saf.drop("Attribute", axis=1, inplace=True)

    # Remove the cleaned GTF
    os.remove(gtf_c)

    return saf

# Run the function and save the output
GTF_to_SAF(gtf_name).to_csv("args.out", sep="\t")
