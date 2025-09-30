#!/usr/bin/env nextflow

process COLLECTOUT {
	tag "${samples}"
	label 'process_low'
	publishDir "results/${samples}", mode:'copy'
	input:
		val(samples)
		tuple val(samples),file(jaffal_out), file(longgf_out), file(geneion_out), file(ctat_out), file(fusionseeker_out), file(coverage_out), file(mosdepth_coverage)
	output:
		tuple val(samples), file("${samples}.xlsx")
	script:
	"""   
	merge_result_260925.py ${samples} ${jaffal_out} ${longgf_out} ${geneion_out} ${ctat_out} ${fusionseeker_out} ${coverage_out}  ${mosdepth_coverage} ${samples}.xlsx
	"""
}
