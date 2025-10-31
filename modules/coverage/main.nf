#!/usr/bin/env nextflow

process COVERAGE {
	tag "${samples}"
	label 'process_medium'
	//publishDir "results/${samples}", mode:'copy'
	input:
		tuple val(samples), file(sorted_bam), file(sorted_bai)
		file(bed_coverage)
	output:
		tuple val(samples), path("${samples}_cov_bedtools.txt"),  emit: bam_cov
	script:
	"""
	bedtools coverage -counts -a ${bed_coverage}  -b  ${sorted_bam}  > ${samples}_cov.txt
	awk 'BEGIN {OFS="\t"; print "chr","start","end","region","counts"} {print \$1,\$2,\$3,\$4,\$6}' ${samples}_cov.txt > ${samples}_cov_bedtools.txt

	"""
}

