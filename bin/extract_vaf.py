#!/usr/bin/env python3

import sys
import csv
import pandas as pd

vcf_file = sys.argv[1]
outputFile = sys.argv[2]

# Read VCF file into a DataFrame
vcf_df = pd.read_csv(vcf_file, sep='\t', comment='#', header=None)
#using comment"#" it will read anything starting with'#' should be treated as comments && header=None indicates that the file doesn't contain column headers
#basically removing #info and later specifiying the column name

#for extracting the specified first 5 columns
extracted_df = vcf_df.iloc[:, :5]

#To extract ref_count and alt_count from DP
format = vcf_df[8].iloc[0]
ad_index = format.split(':').index('AD')
sample_values = vcf_df[9]
ref_count = sample_values.str.split(':').str[ad_index].str.split(',').str[0]
alt_count = sample_values.str.split(':').str[ad_index].str.split(',').str[1]

extracted_df['ref_count'] = pd.Series(ref_count).astype(int)
extracted_df['alt_count'] = pd.Series(alt_count).astype(int)

#Add the VAF% values to the extracted DataFrame
extracted_df['VAF%'] = (extracted_df['alt_count'] / ( extracted_df['alt_count'] + extracted_df['ref_count'])) * 100

# Set cols
extracted_df.columns = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'ref_count', 'alt_count', 'VAF%']
extracted_df.to_csv(outputFile, sep =',', header=True, index=False)