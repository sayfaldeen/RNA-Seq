#!/usr/bin/python3

import numpy as np
import pandas as pd

# Create the argument parser
import argparse as ap
parser = ap.ArgumentParser(description="Script to convert GTF to SAF.")

# Add the arguments
parser.add_argument("-t", "--tsv", dest="tsv_file_name",
        help="path/to/gtf")

parser.add_argument("-o", "--output", dest="out",
        help="Specify name of output file")

args = parser.parse_args()

# Ensure that an input file was given
if not args.tsv_file_name or args.tsv_file_name.split(".")[-1].lower() != "tsv":
    parser.error("Please specify a TSV file. See --help for options.")

# If not output name was given, automatically create one
if not args.out:
    args.out = args.tsv_file_name.replace("tsv", "saf")
    print(f"No output file name was given. '{args.out}' automatically generate as output file name")

# Create the function to convert GTF to SAF
def TSV_to_SAF(f):
    """
    Function to convert GTF to SAF

    Parameters
    ----------
    tsv: TSV file as DF object

    Returns
    -------
    SAF file as DF object

    """
    
    tsv = pd.read_csv(f, sep="\t")


    # For the SAF file, we need:
        # 1: GeneID [?] [details_old_names]
        # 2: Chr [x] [Chrom]
        # 3: Start [x] [5']
        # 4: End [x] [3']
        # 5: Strand [x]
    # Can be unordered

    saf = tsv[["details_old_names",
     "chrom",
     "5'",
     "3'",
     "strand"]]

    saf.columns = [["GeneID", "Chr", "Start", "End", "Strand"]]
    
    return saf

# Run the function and save the output
TSV_to_SAF(args.tsv_file_name).to_csv(args.out, sep="\t", 
        index=False)
