#!/usr/bin/env nextflow

process MINIMAP_ALIGN {
	tag "${samples}"
	label 'process_medium'
	publishDir "results/${samples}", mode:'copy'
	input:
		tuple val(samples), path(reads)
		path(reference_genome)
	output:
		tuple val(samples), path("${samples}_hg38_minimap.paf"),  emit: minimap_paf
	script:
	"""
	minimap2 -c -t 128 -x map-ont ${reference_genome} ${reads}  > ${samples}_hg38_minimap.paf
	"""
}