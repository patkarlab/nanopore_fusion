#!/usr/bin/env nextflow

process SEQKIT{
	tag "${samples}"
	label 'process_medium'
	//publishDir "results/${samples}", mode:'copy'
	input:
		tuple val(samples), path(reads)
	output:
		tuple val(samples), file("${samples}_seqkit.tsv"),  emit: seqkit_tsv
	script:
	"""
	seqkit fx2tab ${reads} -l -g -q -n -i -H > ${samples}_seqkit.tsv
	"""
}
