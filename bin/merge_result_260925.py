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
fusionseeker = args[6]
coverage = args[7]
mosdepth_coverage = args[8]
outfile = args[9]

csvfilenames = [ jaffal, longgf, geneion, ctat, fusionseeker, coverage, mosdepth_coverage]
writer = pd.ExcelWriter(outfile)
for csvfilename in csvfilenames:
        if (os.path.exists(csvfilename)) and (os.path.getsize(csvfilename) != 0):
                sheetname=os.path.split(csvfilename)[1]
                df = pd.read_csv(csvfilename, sep = "\t")
                print('process file:', csvfilename, 'shape:', df.shape)
                new_sheet_name = os.path.splitext(sheetname)[0]
                new_sheet_name = re.sub (sample,"", new_sheet_name, flags = re.IGNORECASE)
                new_sheet_name = re.sub (r"_", "", new_sheet_name, 1)
                df.to_excel(writer,sheet_name=new_sheet_name, index=False)
writer.close()

