#!/usr/bin/env python3
import sys
import pandas as pd

input1 = sys.argv[1] # file extracted from somatic.vcf
input2 = sys.argv[2] # file extracted from somatic_VEP_delheaders
outputFile = sys.argv[3]

#To read the input files
df1 = pd.read_csv(input1, sep=',')
df2 = pd.read_csv(input2, sep=',')
ENST_values_to_accept = ['ENST00000318560']
# Merge based on columns
merge = df2.merge(df1, how = 'inner', left_on = ['CHROM','POS'], right_on = ['CHROM', 'POS'])
output = merge.drop_duplicates() #to remove duplicates
output.rename(columns= {'SYMBOL':'GENE'}, inplace = True)	# change column name from SYMBOL to GENE
extractedData = output[['CHROM',
						'POS',
						'REF',
						'ALT',
						'GENE',
						'ID',
						'ref_count',
						'alt_count',
						'VAF%',
						'AF',
						'MAX_AF_POPS',
						'ALLELE_NUM',
						'VARIANT_CLASS',
						'CLIN_SIG',
						'Feature',
						'Feature_type',
						'Consequence',
						'HGVSc',
						'HGVSp',
						'cDNA_position',
						'CDS_position',
						'Protein_position',
						'Amino_acids',
						'Codons',
						'ENSP'
						]]
columns_to_replace = ['AF', 'MAX_AF_POPS', 'Amino_acids', 'Codons', 'ENSP','CLIN_SIG']
extractedData[columns_to_replace] = extractedData[columns_to_replace].replace('-', '-1')
extractedData['Protein_position'] = extractedData['Protein_position'].apply(lambda x: '-1' if x == '-' else x)
extractedData['cDNA_position'] = extractedData['cDNA_position'].apply(lambda x: '-1' if x == '-' else x)
extractedData['CDS_position'] = extractedData['CDS_position'].apply(lambda x: '-1' if x == '-' else x)

#extractedData = extractedData [extractedData['Feature'].isin(ENST_values_to_accept)]
extractedData.to_csv(outputFile, sep = '\t', header=True, index=False)
