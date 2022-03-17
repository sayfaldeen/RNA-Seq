#!/share/apps/anaconda3-97/bin/python3

# <Summary> Script to count concordant alignemnts that aligned uniquely (EXACTLY 1 time)
import numpy as np
import pandas as pd
import os
import subprocess as sp
import argparse as ap
import matplotlib.pyplot as plt


# Create the argument parser
parser = ap.ArgumentParser(description="Script to count and visualize the number of reads that aligned concordantly EXACTLY ONCE, according to BowTie2.")

parser.add_argument("-i", "--input-dir", dest="idir",
        help="Input directory to scan for sam.gz files.")
parser.add_argument("-o", "--out", dest="out", default=None,
        help="Path to save image out to (.tiff, .png, .jpeg). Default: None")

args = parser.parse_args()


# Create the function to pull exact concordant alignments
def CountExactConcs(sam_name):
    """
    Function to extract the EXACT concordant alignments from the SAM.gz files

    Parameters
    ----------
    sam_name: Name of sam.gz file

    Returns
    -------
    Tuple containing (filename without extensions, exact concordant alignment percentage)
    """

    # Pull the total read numbers
    tot = int(sp.getoutput(f"pigz -d -p 4 -c {sam_name} | grep -v '@' -c"))

    # Pull the exact concordant alignments
    exact_concs = int(sp.getoutput(f"pigz -d -p 4 -c {sam_name} | grep 'YT:Z:CP' | grep -v 'XS' -c"))

    # Calculate the percentage of exact concordant alignments
    perc = round((exact_concs/tot)*100, 2)

    print(f"Percent exact concordant alignemnts for {sam_name}: {perc}%")

    return(sam_name.split("/")[-1].replace("_paired.sam.gz", ""), perc)

# Create function to graph the exact concordant alignments
def GraphExactConcs(results_dict, out=None, dpi=200):
    """
    Function to graph the results for multiple files from the CountExactConcs() function

    Parameters
    ----------
    results_dict: Dictionary containing the results from CountExactConcs()
    out: Path to save image as (.tiff, .png, .jpeg). Must point to existing directory

    Returns
    -------
    Graph that can optionally be saved

    """
    from stats import mean, median, mode

    plt.style.use("ggplot")
    plt.figure(figsize=(10,8))

    plt.bar(x=rdict.keys(), height=rdict.values(),
            color='gold', linewidth=1, edgecolor='k')

    plt.title("Reads that aligned concordantly EXACTLY once per sample\n")
    plt.ylabel("Percent of total reads\n")
    plt.xlabel("\nSample")

    plt.ylim(0,70)

    plt.xticks(rotation=45, ha='right')

    # Pull descriptive stats
    mu = round(mean(rdict.values()),2)
    med = round(median(rdict.values()),2)
    minv = round(min(rdict.values()),2)
    maxv = round(max(rdict.values()),2)

    # Pull plot mins and maxes
    xmin, xmax = plt.gca().get_xlim()
    # ymin, ymax = plt.gca().get_ylim(); manually set now

    # Write descriptive stats out to image
    plt.text(x=xmin*0.9, y=60,
            s=f"""Mean: {mu} -- Median: {med}
Min: {minv} -- Max: {maxv}""",
            va='top', ha="left",
            bbox=dict(facecolor="wheat", edgecolor='k'))

    if out:
        plt.savefig(out, dpi=dpi, bbox_inches='tight')
    else:
        plt.show()


# Pull the sam.gz files from the input dir
sams = [args.idir+x for x in os.listdir(args.idir) if x.endswith("sam.gz")]

res = set() # Prep the set to store the tuple of results

# Iterate through the sam files and calculate the alignment percs
for sam in sams:
    res.add(CountExactConcs(sam))

rdict = dict(res)
del res

# Graph the results
GraphExactConcs(rdict,args.out)
