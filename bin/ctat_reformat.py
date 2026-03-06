#!/usr/bin/env python3

import pandas as pd
import os, sys

args = sys.argv
sample_name=args[1]
ctat_out = args[2]
outfile = args[3]

if (os.path.exists(ctat_out)) and (os.path.getsize(ctat_out) != 0):
    df=pd.read_csv(ctat_out, sep="\t", usecols=[1,2,4,5,7] )

    if not df.empty:
        df[['chrom1','pos1','strand1']]=df['LeftBreakpoint'].str.split(':', expand=True)
        df[['chrom2','pos2','strand2']]=df['RightBreakpoint'].str.split(':', expand=True)
    else:
        # Add empty columns so later concat won’t fail
        for col in ['chrom1','pos1','strand1','chrom2','pos2', 'strand2']:
            df[col] = []
        print("ctat output is empty (only headers). Added required empty columns.")


    df=df.drop(columns=['LeftBreakpoint','RightBreakpoint'])

    df=df.rename(columns={'num_LR':'split','LeftGene':'gene1','RightGene':'gene2'})

    df=df[df['split'] > 3]

    df['name']=sample_name
    df['span']='0'
    df=df[['chrom1','pos1','gene1','chrom2','pos2','gene2','name','split','span','strand1','strand2']]

    df['untemplated_insert']=''
    df['comment']='ctat'

    #outfile=os.path.split(ctat_out)[1]
    #outfile=os.path.splitext(ctat_out)[0]
    #outfile=outfile.replace("-LR.fusion_predictions.tsv", "_reformat.tsv")
    df.to_csv(outfile,sep="\t",index=False,header=True)


