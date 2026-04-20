#!/usr/bin/env nextflow

process MINIMAP_ALIGN {
	tag "${samples}"
	label 'process_high'
	//publishDir "results/${samples}", mode:'copy'
	input:
		tuple val(samples), path(reads)
		path(reference_genome)
	output:
		tuple val(samples), path(reads), path("${samples}_hg38_minimap.paf")
	script:
	"""

	#minimap2 -c -t 128 -x map-ont ${reference_genome} ${reads}  > ${samples}_hg38_minimap.paf
	minimap2 -x splice -k14 -c -t 128   -p 0.8 -N 2 --secondary=no ${reference_genome} ${reads}  > ${samples}_hg38_minimap.paf

	"""
}

process MINIMAP_SAM {
    tag "${samples}"
    label 'process_high'
    //publishDir "results/${samples}", mode:'copy'
    input:
        tuple val(samples), path(reads)
        path(reference_genome)
    output:
        tuple val(samples), path("${samples}_hg38_minimap.sam")
    script:
    """

    minimap2 -ax splice -k14 -t 128 --secondary=no  ${reference_genome} ${reads}  > ${samples}_hg38_minimap.sam
    """
}
