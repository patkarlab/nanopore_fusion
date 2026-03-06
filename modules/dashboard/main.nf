#!/usr/bin/env nextflow
process DASHBOARD {
        tag "${samples}"
        label 'process_low'
        publishDir "results/", mode:'copy'
        input:
                tuple val(samples), file(collectout_out)
				file(cytoband_file)
				file(gtf)
        output:
                tuple val(samples), file("${samples}_dashboard.html")
        script:
        """
        fusion_dashboard_v47_noIGV.py ${samples}.xlsx --cytoband ${cytoband_file} --gtf ${gtf}
        """
}

