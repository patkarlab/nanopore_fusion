#!/usr/bin/env python3
import pandas as pd
import os, sys
import re

args = sys.argv
output = args[1]	# Name of the output file

tsvfilenames=['coverage', 'longshot', 'NanoCaller', 'ClairsTO']
filenameindex=0
writer = pd.ExcelWriter(output)
for arguments in args[2:]:	# The first 2 arguments being script name and the output file
	if os.path.getsize(arguments) != 0:
		df = pd.read_csv(arguments, sep = '\t')
	else:
		df = pd.DataFrame()
	df.to_excel(writer, sheet_name=tsvfilenames[filenameindex], index = False)
	filenameindex += 1
writer.close()
