#!/usr/bin/env nextflow

process COVERAGE {
	tag "${samples}"
	label 'process_medium'
	publishDir "results/${samples}", mode:'copy'
	input:
		tuple val(samples), file(sorted_bam)
		file(bed_coverage)
	output:
		tuple val(samples), path("${samples}_bam_cov.txt"),  emit: bam_cov
	script:
	"""
	bedtools coverage -counts -a ${bed_coverage}  -b  ${sorted_bam}  > ${samples}_bam_cov.txt
	"""
}

