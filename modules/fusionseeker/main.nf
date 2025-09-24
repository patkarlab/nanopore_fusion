#!/usr/bin/env nextflow

process FUSIONSEEKER {
	tag "${samples}"
	label 'process_medium'
	publishDir "results/${samples}", mode:'copy'
	input:
		tuple val(samples), file( coord_sorted_bam)
		tuple val(samples), file (bam_index)
		file(reference_genome)
		file(gtfv104)
	output:
		tuple val(samples), file("${samples}_fusionseeker_final.tsv"),  emit: fus_out
	script:
	"""
	fusionseeker --bam ${samples}_coord_sorted.bam --datatype nanopore --minsupp 5 --thread 64 -o ${samples}_fusionseeker_out --gtf ${gtfv104} --ref ${reference_genome} 
	awk 'BEGIN{OFS="\t"} { print \$2,\$3,\$4,\$5,\$6,\$7,\$8}' ${samples}_fusionseeker_out/confident_genefusion.txt  > ${samples}_fusionseeker_out/${samples}_fusionseeker_final.tsv
	cp ${samples}_fusionseeker_out/${samples}_fusionseeker_final.tsv ${samples}_fusionseeker_final.tsv
	"""

}

