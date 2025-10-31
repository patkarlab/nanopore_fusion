#!/usr/bin/env nextflow

process REFORMATFUSVIZ {
	tag "${samples}"
	label 'process_medium'
	publishDir "results/${samples}", mode:'copy'
	input:
		val(samples)
		tuple val(samples), file(longgf_out), file(genion_out), file(ctat_out), file(jaffal_out)
	output:
		tuple val(samples), path("${samples}_fusviz_in.tsv") 
	script:
	"""
	longgf_format.py ${samples} ${samples}_longgf_final.tsv ${samples}_longgf_reformat.tsv
	genion_format.py ${samples}_genion_final.tsv ${samples}_genion_reformat.tsv
	ctat_reformat.py ${samples} ${samples}_ctat-LR.fusion_predictions.tsv ${samples}_ctat_reformat.tsv
	#fusionseeker_format.py ${samples} ${samples}_fusionseeker_final.tsv ${samples}_fusionseeker_reformat.tsv
	jaffal_reformat.py -s ${samples} -i ${samples}_jaffal.tsv -o ${samples}_jaffal_reformat.tsv
	merge_reformat.py ${samples} ${samples}_fusviz_in.tsv ${samples}_longgf_reformat.tsv ${samples}_genion_reformat.tsv ${samples}_ctat_reformat.tsv ${samples}_jaffal_reformat.tsv
	"""
}
