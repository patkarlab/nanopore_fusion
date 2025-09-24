#!/usr/bin/env python3

import pandas as pd
import os, sys
import re


#run this script on SQK-NBD114-24_barcode03_genion_final.tsv

args = sys.argv
genion_out = args[1]
sample_name = re.sub ("_genion_final.tsv", "", genion_out)
outfile = args[2]

if (os.path.exists(genion_out)) and (os.path.getsize(genion_out) != 0):
    df=pd.read_csv(genion_out, sep="\t", usecols=[0,1,5] )
    print(df)
    if not df.empty:
        df[['gene1', 'gene2']] = df['gene'].str.split('::', expand=True)
        df[['break1', 'break2', 'null']] = df['ranges'].str.split(';', expand=True)

        df[['chrom1', 'position1']] = df['break1'].str.split(':', expand=True)
        df[['pos1', 'pos1_end']] = df['position1'].str.split('-', expand=True)

        df[['chrom2', 'position2']] = df['break2'].str.split(':', expand=True)
        df[['pos2', 'pos2_end']] = df['position2'].str.split('-', expand=True)

    else:
    # Add empty columns so later concat won’t fail
        for col in ['gene1','gene2','chrom1','pos1','chrom2','pos2']:
            df[col] = []
        print("genion output is empty (only headers). Added required empty columns.")
        

    #df=df.drop(columns=['gene','ranges','null','break1', 'position1','pos1_end'])

    df=df.rename(columns={'num_reads':'split'})
    df=df[df['split'] > 3]

    df['name']=sample_name

    df=df[['chrom1','pos1','gene1','chrom2','pos2','gene2','name','split']]

    df['span']='0'
    df['strand1']='.'
    df['strand2']='.'
    df['untemplated_insert']=''
    df['comment']='genion'

    df['chrom1'] = 'chr' + df['chrom1'].astype(str)
    df['chrom2'] = 'chr' + df['chrom2'].astype(str)
    print(df)

    #outfile=os.path.split(genion_out)[1]
    #outfile=os.path.splitext(genion_out)[0]
    #outfile=outfile+"_reformat.tsv"
    df.to_csv(outfile,sep="\t",index=False,header=True)

