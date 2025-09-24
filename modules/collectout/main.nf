#!/usr/bin/env nextflow

process COLLECTOUT {
	tag "${samples}"
	label 'process_low'
	publishDir "results/${samples}", mode:'copy'
	input:
		val(samples)
		tuple val(samples), file(longgf_out), file(geneion_out), file(ctat_out), file(fusionseeker_out), file(coverage_out)
	output:
		tuple val(samples), file("${samples}_compiled.xlsx")
	script:
	"""   
	merge_result_script.py ${samples} ${samples}_compiled.xlsx  ${longgf_out} ${geneion_out} ${ctat_out} ${fusionseeker_out} ${coverage_out} 
	"""
}