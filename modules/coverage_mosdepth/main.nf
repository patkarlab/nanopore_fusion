#!/usr/bin/env nextflow

process COVERAGE_MOSDEPTH {
        tag "${samples}"
        label 'process_medium'
        //publishDir "results/${samples}/", mode:'copy'
        input:
                tuple val(samples), path(sorted_bam), path(bam_index)
				path(bed_coverage)
        output:
                tuple val(samples), path("${samples}.thresholds.bed.gz"), path("${samples}.regions.bed.gz")
		
        script:
        """
        mosdepth --by ${bed_coverage} --thresholds 1,10,30,100,200 ${samples} ${sorted_bam}
		"""
}


process MOSDEPTH_SUMMARY {
        tag "${samples}"
        label 'process_medium'
        //publishDir "results/${samples}", mode:'copy'
        input:
                tuple val(samples), path(thr_file), path(region_file)

        output:
                tuple val(samples), path("${samples}_cov_mosdepth.tsv")
        script:
        """
        mosdepth_coverage.py -t ${thr_file} -r ${region_file} -o ${samples}_cov_mosdepth.tsv

        """
}


