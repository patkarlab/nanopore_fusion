#!/usr/bin/env nextflow

process LONGGF {
	tag "${samples}"
	label 'process_medium'
	publishDir "results/${samples}", mode:'copy'
	input:
		tuple val(samples), file (sorted_bam)
		file(gtf)
	output:
		tuple val(samples), file ("${samples}_longgf_final.tsv"), emit: longgf_out
	script:
	"""
	LongGF  ${sorted_bam} ${gtf} 100 100 100 > ${samples}_longgf.log
	awk 'BEGIN {OFS="\\t"; print "gene", "num_reads", "breakpoint1", "breakpoint2"} /SumGF/ {print \$2, \$3, \$4, \$5}' ${samples}_longgf.log > ${samples}_longgf_final.tsv
	"""
}
