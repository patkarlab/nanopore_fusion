#!/usr/bin/bash

sample=$1
vcf_file=$2
output_name=$3

source deactivate
no_of_variants=$(grep -v "^#" ${vcf_file} | wc -l )
if [ ${no_of_variants} -gt 0 ]; then
	# Extract vaf,af,alt and ref count
	/home/diagnostics/pipelines/Nanopore_KDM/scripts/extract_vaf.py ${vcf_file} ${sample}.extracted.csv

	# Annotate using vep
	vep -i ${vcf_file} --cache -o ${sample}_vep.txt --offline --tab --force_overwrite --symbol --protein --af --max_af  --no_check_alleles --sift b --variant_class --canonical --allele_number --hgvs --shift_hgvs 1 --af_1kg --af_gnomadg --pubmed
	grep -v "##" ${sample}_vep.txt > ${sample}_vep_delheaders.txt

	# Extract the required columns from vep data
	/home/diagnostics/pipelines/Nanopore_KDM/scripts/extract_vepdata.py ${sample}_vep_delheaders.txt ${sample}.extractedvepdelheaders.csv

	# Merge vep and and vcf data
	/home/diagnostics/pipelines/Nanopore_KDM/scripts/mergeSomaticVep.py ${sample}.extracted.csv ${sample}.extractedvepdelheaders.csv ${output_name}

else
	touch ${output_name}
fi