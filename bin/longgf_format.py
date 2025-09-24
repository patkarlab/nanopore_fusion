#!/usr/bin/env python3

import pandas as pd
import os, sys
import re

#run this script on SQK-NBD114-24_barcode03_longgf_final.tsv

args = sys.argv
sample_name=args[1]
longgf_out = args[2]
outfile=args[3]
#sample_name = re.sub ("_longgf_final.tsv", "", longgf_out)

if (os.path.exists(longgf_out)) and (os.path.getsize(longgf_out) != 0):
    df=pd.read_csv(longgf_out, sep="\t")
    
    if not df.empty:
        df[['gene1', 'gene2']] = df['gene'].str.split(':', expand=True)
        df[['chrom1', 'pos1']] = df['breakpoint1'].str.split(':', expand=True)
        df[['chrom2', 'pos2']] = df['breakpoint2'].str.split(':', expand=True)
    else:
    # Add empty columns so later concat won’t fail
        for col in ['gene1','gene2','chrom1','pos1','chrom2','pos2']:
            df[col] = []
        print("longgf output is empty (only headers). Added required empty columns.")


    df=df.rename(columns={'num_reads':'split'})
    df=df[df['split'] > 3]

    df['name']=sample_name
    df=df[['chrom1','pos1','gene1','chrom2','pos2','gene2','name','split']]
    df['span']='0'
    df['strand1']='.'
    df['strand2']='.'
    df['untemplated_insert']=''
    df['comment']='longgf'

    #outfile=os.path.split(longgf_out)[1]
    #outfile=os.path.splitext(longgf_out)[0]
    #outfile=outfile.replace("_final.tsv", "_reformat.tsv")
    df.to_csv(outfile,sep="\t",index=False,header=True)
