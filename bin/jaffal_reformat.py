#!/usr/bin/env python3

import pandas as pd
import argparse

parser = argparse.ArgumentParser(description="This script takes the sample_jaffal.tsv as input file and reformats it")

parser.add_argument("-s", "--sample_name", required=True, help="sample_id")
parser.add_argument("-i","--jaffal_out", required=True, help="Path to jaffal_out.tsv")
parser.add_argument("-o","--outfile", required=True, help="Path to outfile")

args = parser.parse_args()

#read the input file
df=pd.read_csv(args.jaffal_out, sep="\t" )

#split the 'fusion genes' columns
if not df.empty:
        df[['gene1', 'gene2']] = df['fusion genes'].str.split(':', expand=True)

else:
    # Add empty columns so later concat won’t fail
        for col in ['fusion genes',	'chrom1','base1','strand1',	'chrom2','base2','strand2',	'spanning reads',
				'contig break','classification','known']:
            df[col] = []
        print("jaffal output is empty (only headers). Added required empty columns.")

#rename columns
df.rename(columns={'base1': 'pos1', 'base2': 'pos2','spanning reads':'split'  }, inplace=True)

#add columns
df['name']=args.sample_name
df['span']='0'
df['strand1']='.'
df['untemplated_insert']=''
df['comment']='jaffal'

#select and rearrange columns
df=df[['chrom1','pos1','gene1','chrom2','pos2','gene2','name','split','span','strand1','strand2','untemplated_insert','comment']]


df.to_csv(args.outfile,sep="\t",index=False,header=True)


