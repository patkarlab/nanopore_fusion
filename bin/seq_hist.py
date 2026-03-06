#!/usr/bin/env python3

import matplotlib.pyplot as plt
import pandas as pd
import sys
import os

args= sys.argv
sample_name=args[1]
seq_summary= args[2]
hist_out=args[3]

if (os.path.exists(seq_summary)) and (os.path.getsize(seq_summary) != 0):
    df = pd.read_csv(seq_summary,sep="\t", comment='#',names=['readID','length','GC','avgQual'])
    
    bins = [0, 100, 200, 300, 400, 500, 600,700,800, 900,1000,1100,1200,1300, 1400,1500,1600,1700,1800,1900,2000]
    plt.hist(df.length, bins=bins,color='cadetblue', edgecolor='black')

    plt.xlabel('Read Length(bases)')
    plt.ylabel('Read Count')
    plt.title('Histogram of read lengths')

    plt.savefig(hist_out)

#this script limits the bin range to 2000 basepairs, can be changed by changing the bins parameter



