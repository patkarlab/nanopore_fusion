#!/usr/bin/env python3

import pandas as pd
import os, sys


#run this script on confident_genefusion.txt in same directory

args = sys.argv
sample_name=args[1]
fusionseeker_out = args[2]
outfile=args[3]

if (os.path.exists(fusionseeker_out)) and (os.path.getsize(fusionseeker_out) != 0):
    df=pd.read_csv(fusionseeker_out, sep="\t", usecols=[0,1,2,3,4,5,6] )


    df=df.rename(columns={'Gene1':'gene1','Gene2':'gene2','NumSupp':'split','Chrom1':'chrom1', 'Breakpoint1':'pos1','Chrom2':'chrom2','Breakpoint2':'pos2'})
    df=df[df['split'] > 5]
    df=df.sort_values(by='split',ascending=False)

    df['name']=sample_name

    df=df[['chrom1','pos1','gene1','chrom2','pos2','gene2','name','split']]

    df['span']='0'
    df['strand1']='.'
    df['strand2']='.'
    df['untemplated_insert']=''
    df['comment']='fusionseeker'

    #outfile=os.path.split(fusionseeker_out)[1]
    #outfile=os.path.splitext(fusionseeker_out)[0]
    #outfile=outfile.replace("_final.tsv", "_reformat.tsv")
    df.to_csv(outfile,sep="\t",index=False,header=True)


