/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { softwareVersionsToYAML } 	from '../subworkflows/nf-core/utils_nfcore_pipeline/main'
include { BASECALL }              from '../modules/basecall/main'
include { SEQKIT } 		 			from '../modules/seqkit/main'
include { HISTOGRAM }    			from '../modules/histogram/main'
include { DORADO_ALIGN }	 		from '../modules/dorado_align/main'
include { MINIMAP_SAM}     			from '../modules/minimap_align/main'
include { LONGGF }		 			from '../modules/longgf/main'
include { SAMTOOLS_SORT } 	 		from '../modules/samtools_sort/main'
include { COVERAGE }       	 		from '../modules/coverage/main'
include { COVERAGE_MOSDEPTH }		from '../modules/coverage_mosdepth/main'
include { MOSDEPTH_SUMMARY }  		from '../modules/coverage_mosdepth/main'
include { MINIMAP_ALIGN }        	from '../modules/minimap_align/main'
include { GENION }	         		from '../modules/genion/main'
include { CTAT}                 	from '../modules/ctat/main'
include { JAFFAL}               	from '../modules/jaffal/main'
include { COORD_SORT} 				from '../modules/coord_sort/main'
//include { FUSIONSEEKER} 			from '../modules/fusionseeker/main'
include { COLLECTOUT}           	from '../modules/collectout/main'
include { DASHBOARD }				from '../modules/dashboard/main'
include { REFORMATFUSVIZ }			from '../modules/reformatFusviz/main'
include { Longshot;Nanocaller;ClairsTO;Clair3rna;Longcallr;VarDict;CombineCallers }                from '../modules/variant_calling/main.nf'
include { vcf_filter as vcf_filter_longshot;vcf_filter as vcf_filter_longcallr }         from '../modules/variant_calling/main.nf'
include { ANNOVAR as ANNOVAR_longshot;ANNOVAR as ANNOVAR_nanocaller;ANNOVAR as ANNOVAR_ClairsTO;ANNOVAR as ANNOVAR_clair3rna;ANNOVAR as ANNOVAR_longcallr } from '../modules/variant_calling/main.nf'
include { Combine_variants } 		from '../modules/variant_calling/main.nf'
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

reference_genome = file("${params.genome}", checkIfExists: true)
ind_files = file("${params.genome_dir}/${params.ind_files}.*")
gtf = file("${params.gtf}", checkIfExists: true)
genion_cdna = file("${params.cdna}", checkIfExists: true)
genion_superdups = file("${params.superdups}", checkIfExists: true)
//gtfv104 = file("${params.gtf_2}", checkIfExists: true)
ctat_genome_lib_dir = file("${params.genome_lib}", checkIfExists: true)
bed_coverage = file("${params.bed_file}", checkIfExists: true )
cytoband_file = file("${params.cytoband}", checkIfExists: true)
dorado_models_dir = file("${params.dorado_models_dir}", checkifExists: true)
samplesheet_file = file(params.input)

longshot = params.longshot
nanocaller = params.nanocaller
clairsto = params.clairsto
clair3rna = params.clair3rna
longcallr = params.longcallr


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
	//get read lengths
	SEQKIT(samples_ch )

	//plot read length histogram from seqkit tsv
	HISTOGRAM(SEQKIT.out )

	// NOT TO  BE USED - align reads to hg38 reference genome using dorado
	//DORADO_ALIGN(samples_ch, reference_genome )

	//NOT TO  BE USED - sort dorado bam by read names using samtools
	//SAMTOOLS_SORT( DORADO_ALIGN.out.dorado_bam )

	//align reads to reference using minimap
	MINIMAP_SAM(samples_ch, reference_genome )
	SAMTOOLS_SORT( MINIMAP_SAM.out )

	// calling fusion on dorado name-sorted bam
	LONGGF( SAMTOOLS_SORT.out.sorted_bam, gtf )

	// align fastq with minimap to get paf file for genion
	MINIMAP_ALIGN (samples_ch, reference_genome )

	// call fusions on paf using genion
	GENION (MINIMAP_ALIGN.out, gtf, genion_cdna, genion_superdups )

	//ctat-lr-fusion for RNA fusion
	CTAT(samples_ch, ctat_genome_lib_dir ) 

	//jaffal fusioncaller
	JAFFAL(samples_ch )

	//sort dorado bam on coordinates using samtools
	COORD_SORT(MINIMAP_SAM.out )

	//call RNA fusion on dorado coord-sorted bam
	//FUSIONSEEKER(COORD_SORT.out.coord_sorted_bam, COORD_SORT.out.bam_index, reference_genome, gtfv104 )

	//calculate coverage over target regions from sorted bam using bedtools
	COVERAGE(COORD_SORT.out, bed_coverage )

	// Calculate coverage over set thresholds in target regions from sorted, indexed bam using mosdepth
	COVERAGE_MOSDEPTH(COORD_SORT.out, bed_coverage ) | MOSDEPTH_SUMMARY

	//script to collect output from all fusioncallers into one spreadsheet per sample
	COLLECTOUT (samples_ch, JAFFAL.out.jaffal_tsv.join(LONGGF.out.longgf_out.join(GENION.out.genion_tsv.join(CTAT.out.ctat_out.join(COVERAGE.out.join(MOSDEPTH_SUMMARY.out))))))

	//script to create dashboard from fusion caller outputs
	DASHBOARD ( COLLECTOUT.out, cytoband_file, gtf)

	//variant calling and annotation
	Longshot(COORD_SORT.out, reference_genome, ind_files)
	vcf_filter_longshot(Longshot.out, bed_coverage, longshot)
	ANNOVAR_longshot(vcf_filter_longshot.out, longshot)
	Nanocaller(COORD_SORT.out, reference_genome, ind_files, bed_coverage )
	ANNOVAR_nanocaller(Nanocaller.out, nanocaller)
	ClairsTO(COORD_SORT.out, reference_genome, ind_files, bed_coverage)
	ANNOVAR_ClairsTO(ClairsTO.out, clairsto)
	Clair3rna(COORD_SORT.out, reference_genome, ind_files, bed_coverage)
	ANNOVAR_clair3rna(Clair3rna.out, clair3rna)
	Longcallr(COORD_SORT.out, reference_genome, ind_files)
	vcf_filter_longcallr(Longcallr.out, bed_coverage, longcallr)
	ANNOVAR_longcallr(vcf_filter_longcallr.out, longcallr)
	//VarDict(COORD_SORT.out, reference_genome, ind_files, bed_coverage)


	Combine_variants(ANNOVAR_ClairsTO.out.join(ANNOVAR_nanocaller.out.join(ANNOVAR_longshot.out.join(ANNOVAR_clair3rna.out.join(ANNOVAR_longcallr.out)))))

	//VEP steps
	//Clair_annotate(ClairsTO.out)
	//CombineCallers(COVERAGE.out.join(Clair_annotate.out.join(Nanocaller.out.join(Longshot.out))))

	//script to merge and aggregate fusioncaller output to be used as fusviz input
	//REFORMATFUSVIZ (samples_ch, LONGGF.out.join(GENION.out.join(CTAT.out.ctat_out.join(JAFFAL.out))))



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

workflow NANOPORE_FUSION_BASECALL {

	main:
	pod5_ch = Channel.fromPath(params.pod5_dir)
	samplesheet_ch = Channel.fromPath(params.input)
	.splitCsv(header:false)
	.map { row ->
		tuple(row[0], row[1])
	}

	// Basecalling + demultiplexing
	BASECALL(pod5_ch, dorado_models_dir, samplesheet_file)

	samples_ch = BASECALL.out
        .flatten()
        .map { fq ->
            def sample = fq.name.replace(".fastq.gz","")
            tuple(sample, fq)
        }

	//get read lengths
	SEQKIT(samples_ch )

	//plot read length histogram from seqkit tsv
	HISTOGRAM(SEQKIT.out )

	//NOT TO  BE USED - align reads
	//DORADO_ALIGN(samples_ch, reference_genome )

	//NOT TO  BE USED - sort dorado bam by read names using samtools
	//SAMTOOLS_SORT( DORADO_ALIGN.out.dorado_bam )

	//align reads to reference using minimap
        //MINIMAP_SAM(samples_ch, reference_genome )
        //SAMTOOLS_SORT( MINIMAP_SAM.out )

	// calling fusion on dorado name-sorted bam
	//LONGGF( SAMTOOLS_SORT.out.sorted_bam, gtf )

	// align fastq with minimap to get paf file for genion
	//MINIMAP_ALIGN (samples_ch, reference_genome )

	// call fusions on paf using genion
	//GENION (MINIMAP_ALIGN.out, gtf, genion_cdna, genion_superdups )

	//ctat-lr-fusion for RNA fusion
	//CTAT(samples_ch, ctat_genome_lib_dir )

	//jaffal fusioncaller
	//JAFFAL(samples_ch )

	//sort dorado bam on coordinates using samtools
	//COORD_SORT(DORADO_ALIGN.out.dorado_bam )

	//calculate coverage over target regions from sorted bam using bedtools
	//COVERAGE(COORD_SORT.out, bed_coverage )

	// Calculate coverage over set thresholds in target regions from sorted, indexed bam using mosdepth
	//COVERAGE_MOSDEPTH(COORD_SORT.out, bed_coverage ) | MOSDEPTH_SUMMARY

	//script to collect output from all fusioncallers into one spreadsheet per sample
	//COLLECTOUT (samples_ch, JAFFAL.out.jaffal_tsv.join(LONGGF.out.longgf_out.join(GENION.out.genion_tsv.join(CTAT.out.ctat_out.join(COVERAGE.out.join(MOSDEPTH_SUMMARY.out))))))

	//script to create dashboard from fusion caller outputs
	//DASHBOARD ( COLLECTOUT.out, cytoband_file, gtf)

	//variant calling and annotation
        //Longshot(COORD_SORT.out, reference_genome, ind_files)
        //vcf_filter_longshot(Longshot.out, bed_coverage, longshot)
        //ANNOVAR_longshot(vcf_filter_longshot.out, longshot)
        //Nanocaller(COORD_SORT.out, reference_genome, ind_files, bed_coverage )
        //ANNOVAR_nanocaller(Nanocaller.out, nanocaller)
        //ClairsTO(COORD_SORT.out, reference_genome, ind_files, bed_coverage)
        //ANNOVAR_ClairsTO(ClairsTO.out, clairsto)
        //Clair3rna(COORD_SORT.out, reference_genome, ind_files, bed_coverage)
        //ANNOVAR_clair3rna(Clair3rna.out, clair3rna)
        //Longcallr(COORD_SORT.out, reference_genome, ind_files)
        //vcf_filter_longcallr(Longcallr.out, bed_coverage, longcallr)
        //ANNOVAR_longcallr(vcf_filter_longcallr.out, longcallr)
        //VarDict(COORD_SORT.out, reference_genome, ind_files, bed_coverage)

        //Combine_variants(ANNOVAR_ClairsTO.out.join(ANNOVAR_nanocaller.out.join(ANNOVAR_longshot.out.join(ANNOVAR_clair3rna.out.join(ANNOVAR_longcallr.out)))))

	//script to merge and aggregate fusioncaller output to be used as fusviz input
	//REFORMATFUSVIZ (samples_ch, LONGGF.out.join(GENION.out.join(CTAT.out.ctat_out.join(JAFFAL.out))))
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
