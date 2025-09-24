#!/usr/bin/env python3

import pandas as pd
import re
import sys

args = sys.argv
sample = args[1]
outfile = args[2]
longgf = args[3]
geneion = args[4]
ctat = args[5]
fusionseeker = args[6]

#read in the reformatted fusioncaller outputs
df1=pd.read_csv(longgf, sep="\t")
df2=pd.read_csv(geneion, sep="\t")
df3=pd.read_csv(ctat, sep="\t")
df4=pd.read_csv(fusionseeker, sep="\t")

df_to_merge = [df1, df2, df3, df4]

#concatenate
#valid_dfs = [df for df in df_to_merge if not df.empty]
df = pd.concat(df_to_merge, ignore_index=True)


# Function to extract numeric part from chromosome (e.g., 'chr1' -> 1, 'chrX' -> 23, 'chrY' -> 24)
def chr_to_num(chr_str):
    chr_str = chr_str.replace('chr', '')
    if chr_str == 'X':
        return 23
    elif chr_str == 'Y':
        return 24
    elif chr_str == 'M':
        return 25
    else:
        return int(chr_str)

# Apply the function to get numeric chromosome values
df['chr1_num'] = df['chrom1'].apply(chr_to_num)
df['chr2_num'] = df['chrom2'].apply(chr_to_num)

# Swap rows where chrom1 is greater than chrom2
swapped = df['chr1_num'] > df['chr2_num']

# Perform the swap
df.loc[swapped, ['chrom1', 'chrom2']] = df.loc[swapped, ['chrom2', 'chrom1']].values
df.loc[swapped, ['pos1', 'pos2']] = df.loc[swapped, ['pos2', 'pos1']].values
df.loc[swapped, ['gene1', 'gene2']] = df.loc[swapped, ['gene2', 'gene1']].values
df.loc[swapped, ['chr1_num', 'chr2_num']] = df.loc[swapped, ['chr2_num', 'chr1_num']].values

# Drop helper columns
df = df.drop(columns=['chr1_num', 'chr2_num'])

#group rows which have same fusion and keep only the row with max reads
idx = df.groupby(['chrom1','gene1','chrom2','gene2'])['split'].idxmax()
df=df.loc[idx]

#sort df by supporting reads
df = df.sort_values(by='split', ascending=False)
#print df to outfile
df.to_csv(outfile,sep="\t",index=False,header=True)
#print(df)

