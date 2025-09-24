#!/usr/bin/env nextflow

process HISTOGRAM {
	tag "${samples}"
	label 'process_medium'
	publishDir "results/${samples}", mode:'copy'
	input:
		tuple val(samples), file(seqkit_tsv) 
	output:
		tuple val(samples), path("${samples}_read_length_histogram.png")
	script:
	"""
	seq_hist.py ${samples} ${samples}_seqkit.tsv ${samples}_read_length_histogram.png
	 """
}

