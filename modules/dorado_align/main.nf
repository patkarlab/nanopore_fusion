#!/usr/bin/env nextflow

process DORADO_ALIGN {
	tag "${samples}"
	label 'process_medium'
	publishDir "results/${samples}", mode:'copy'    
	input:
		tuple val(samples), path(reads) 
		file(reference_genome)
	output:
		tuple val(samples), file("${samples}.bam"),  emit: dorado_bam
	script:   	
	"""
	dorado aligner -t 128 ${reference_genome} ${reads}  > ${samples}.bam
	"""
}
