#!/usr/bin/env nextflow

process COORD_SORT {
	tag "${samples}"
	label 'process_medium'
	publishDir "results/${samples}", mode:'copy'    
	input:
		tuple val(samples), file(dorado_bam)
	output:
		tuple val(samples), path("${samples}_coord_sorted.bam"),  emit: coord_sorted_bam
		tuple val(samples), path("${samples}_coord_sorted.bam.bai"), emit: bam_index
	script:
	"""
	samtools sort -o ${samples}_coord_sorted.bam ${dorado_bam}
	samtools index ${samples}_coord_sorted.bam
	"""
}

