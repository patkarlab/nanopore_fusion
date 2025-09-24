/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline/main'
include { SEQKIT } 		 from '../modules/seqkit/main'
include { HISTOGRAM }    from '../modules/histogram/main'
include { DORADO_ALIGN }	 from '../modules/dorado_align/main'
include { LONGGF }		 from '../modules/longgf/main'
include { SAMTOOLS_SORT } 	 from '../modules/samtools_sort/main'
include { COVERAGE }       	 from '../modules/coverage/main'
include { MINIMAP_ALIGN }        from '../modules/minimap_align/main'
include { GENION }	         from '../modules/genion/main'
include { CTAT}                 from '../modules/ctat/main'
include { COORD_SORT} 		from '../modules/coord_sort/main'
include { FUSIONSEEKER} 	from '../modules/fusionseeker/main'
include { COLLECTOUT}           from '../modules/collectout/main'
include { REFORMATFUSVIZ }	from '../modules/reformatFusviz/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

reference_genome = file("${params.genome}", checkIfExists: true)
gtf = file("${params.gtf}", checkIfExists: true)
genion_cdna = file("${params.cdna}", checkIfExists: true)
genion_superdups = file("${params.superdups}", checkIfExists: true)
gtfv104 = file("${params.gtf_2}", checkIfExists: true)
ctat_genome_lib_dir = file("${params.genome_lib}", checkIfExists: true)
bed_coverage = file("${params.bed_file}", checkIfExists: true )

workflow NANOPORE_FUSION {
    Channel
        .fromPath(params.input)
        .splitCsv(header:false)
        .map { row -> row[0] }
        .map { sample ->
            def fq = file("${params.sequences}/${sample}.fastq.gz", checkIfExists: false)
            tuple(sample, fq) 
        }
        .set { samples_ch }

    main: 
    //quality check for reads
    SEQKIT(samples_ch)

	//plot read length histogram
	HISTOGRAM(SEQKIT.out)

    // align reads to hg38 reference genome using dorado
    DORADO_ALIGN(samples_ch, reference_genome )
    
    // sort dorado bam by read names using samtools
    SAMTOOLS_SORT( DORADO_ALIGN.out.dorado_bam)

    //calculate coverage over bed file from samtools sorted bam 
    COVERAGE(SAMTOOLS_SORT.out.sorted_bam, bed_coverage )

    // calling fusion on dorado name-sorted bam
    LONGGF( SAMTOOLS_SORT.out.sorted_bam, gtf)
    
    // align fastq with minimap to get paf file for genion
    MINIMAP_ALIGN (samples_ch, reference_genome )

    // call fusion on paf using genion
    GENION (samples_ch, MINIMAP_ALIGN.out.minimap_paf, gtf, genion_cdna, genion_superdups)

    //ctat for RNA fusion
    CTAT(samples_ch, ctat_genome_lib_dir)

    //sort dorado bam on coordinates using samtools
    COORD_SORT(DORADO_ALIGN.out.dorado_bam)

    //call RNA fusion on dorado coord-sorted bam
    FUSIONSEEKER(COORD_SORT.out.coord_sorted_bam, COORD_SORT.out.bam_index, reference_genome, gtfv104)

    //script to collect output from all fusioncallers into one spreadsheet per sample
    COLLECTOUT (samples_ch, LONGGF.out.longgf_out.join(GENION.out.genion_tsv.join(CTAT.out.ctat_out.join(FUSIONSEEKER.out.fus_out.join(COVERAGE.out)))))

    //script to merge and aggregate fusioncaller output to be used as fusviz input
    REFORMATFUSVIZ (samples_ch, LONGGF.out.join(GENION.out.join(CTAT.out.ctat_out.join(FUSIONSEEKER.out))))



    ch_versions = Channel.empty()

    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name:  'nanopore_fusion_software_'  + 'versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }


    emit:
    versions       = ch_versions                 // channel: [ path(versions.yml) ]

}


workflow.onComplete {
    log.info ( workflow.success ? "\n\nDone! Output in the 'results/' directory \n" : "Oops .. something went wrong" )
    println "Total time taken: ${workflow.duration}"
}
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
