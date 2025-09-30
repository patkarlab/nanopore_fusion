#!/usr/bin/env nextflow

process CTAT {
	tag "${samples}"
	label 'process_high'
	publishDir  "results/${samples}/", mode: 'copy'
	input:
		tuple val(samples), path(reads)
		file(ctat_genome_lib_dir)
	output:
		tuple val(samples), file ("${samples}_ctat-LR.fusion.tsv"), emit: ctat_out
		tuple val(samples), file ("${samples}_ctat-LR-fusion.fusion_inspector_web.html"), emit: ctat_html
	script:
	"""
	ctat-LR-fusion -T ${reads}  --genome_lib_dir ${ctat_genome_lib_dir} --vis
	cp ctat_LR_fusion_outdir/ctat-LR-fusion.fusion_predictions.abridged.tsv ${samples}_ctat-LR.fusion.tsv
	cp ctat_LR_fusion_outdir/ctat-LR-fusion.fusion_inspector_web.html ${samples}_ctat-LR-fusion.fusion_inspector_web.html
	"""
}
