#!/usr/bin/env nextflow

process GENION {
	tag "${samples}"
	label 'process_medium'
	publishDir "results/${samples}", mode:'copy'     
	input:
		tuple val(samples), path(reads)
		tuple val(samples), file(minimap_paf)
		file(gtf)
		file(genion_cdna)
		file(genion_superdups)
	output:
		tuple val (samples), path ("${samples}_genion_final.tsv"), emit: genion_tsv
	script:
	"""
	gunzip -c ${reads} > ${samples}.fastq
	genion -i ${samples}.fastq --gtf ${gtf} --gpaf ${samples}_hg38_minimap.paf -s ${genion_cdna} -d ${genion_superdups} -o ${samples}_genion.tsv
	awk 'BEGIN {OFS="\t"; print "gene","num_reads","ffigf-score","FiN-score","normal_counts","ranges"} {print \$2,\$5,\$3,\$4,\$6,\$8}' ${samples}_genion.tsv > ${samples}_genion_final.tsv
	"""
}

