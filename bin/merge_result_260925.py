#!/usr/bin/env python3

import pandas as pd
import os, sys
import re

args = sys.argv
sample = args[1]
jaffal = args[2]
longgf = args[3]
geneion = args[4]
ctat = args[5]
coverage = args[6]
mosdepth_coverage = args[7]
outfile = args[8]

csvfilenames = [jaffal, longgf, geneion, ctat, coverage, mosdepth_coverage]

writer = pd.ExcelWriter(outfile)

for csvfilename in csvfilenames:
    if (os.path.exists(csvfilename)) and (os.path.getsize(csvfilename) != 0):

        sheetname = os.path.split(csvfilename)[1]

        try:
            df = pd.read_csv(csvfilename, sep="\t")
        except pd.errors.EmptyDataError:
            print('WARNING: file has no data rows:', csvfilename, file=sys.stderr)
            df = pd.DataFrame()

        if df.empty:
            print('WARNING: dataframe empty after reading:', csvfilename, file=sys.stderr)

        print('process file:', csvfilename, 'shape:', df.shape)

        new_sheet_name = os.path.splitext(sheetname)[0]
        new_sheet_name = re.sub(sample, "", new_sheet_name, flags=re.IGNORECASE)
        new_sheet_name = re.sub(r"_", "", new_sheet_name, 1)

        df.to_excel(writer, sheet_name=new_sheet_name, index=False)

writer.close()
