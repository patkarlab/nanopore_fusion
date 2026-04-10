#!/usr/bin/env nextflow

process SAMTOOLS_SORT {
	tag "${samples}"
	label 'process_medium'
	//publishDir "results/${samples}", mode:'copy'
	input:
		tuple val(samples), file(dorado_bam)
	output:
		tuple val(samples), path("${samples}_sorted.bam"),  emit: sorted_bam
	script:
	"""
	samtools view -bS ${dorado_bam} | samtools sort -n -o ${samples}_sorted.bam
	#samtools sort -n -o ${samples}_sorted.bam ${dorado_bam}
	"""
}

