#!/usr/bin/env python3

import pandas as pd
import argparse

parser = argparse.ArgumentParser(description="This script takes the sample.thresholds.bed.gz & sample.regions.bed.gz as input files, converts the number of bases at XX depth to percentage of bases at XX depth, adds mean depth from sample.regions.bed.gz, and also adds a column for region lengths")

parser.add_argument("-t","--thresholds_file", required=True, help="Path to sample.thresholds.bed.gz")
parser.add_argument("-r","--regions_cov_file", required=True, help="Path to sample.regions.bed.gz")
parser.add_argument("-o","--outfile", required=True, help="Path to outfile")

args = parser.parse_args()

# Read thresholds file 
thr_df = pd.read_csv(args.thresholds_file, sep="\t", comment="#", header=None)
thr_df.columns = ["chrom","start","end","region","1X","10X","30X","100X","200X"]

# List of columns to convert
columns_to_convert_numeric = ["start","end","1X","10X","30X","100X","200X"]
for col in columns_to_convert_numeric:
    thr_df[col] = pd.to_numeric(thr_df[col], errors='coerce')

# Compute region length
thr_df["region_length"] = thr_df["end"] - thr_df["start"]

# Convert counts -> fractions
for col in ["1X","10X","30X","100X","200X"]:
    thr_df["cov>="+col+"(%)"] = (thr_df[col] / thr_df["region_length"]) * 100

thr_df[['cov>=1X(%)', 'cov>=10X(%)', 'cov>=30X(%)', 'cov>=100X(%)', 'cov>=200X(%)']] = thr_df[['cov>=1X(%)', 'cov>=10X(%)', 'cov>=30X(%)', 'cov>=100X(%)', 'cov>=200X(%)']].round(3)

# Read mean coverage and add column names
mean_df = pd.read_csv(args.regions_cov_file, sep="\t", header=None)
mean_df.columns = ["chrom","start","end","region","mean_coverage"]

# Merge the dfs
summary_df = pd.merge(thr_df, mean_df[["chrom","start","end","mean_coverage"]],
                   on=["chrom","start","end"], how="left")

# Keep required columns
summary_df = summary_df[[
    "chrom","start","end","region","region_length",
    'cov>=1X(%)', 'cov>=10X(%)', 'cov>=30X(%)',
       'cov>=100X(%)', 'cov>=200X(%)',"mean_coverage"]]

# Save the output
summary_df.to_csv(args.outfile, sep="\t", index=False)



