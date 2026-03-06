#!/usr/bin/env nextflow
nextflow.enable.dsl=2


process Longshot {
	//publishDir "$PWD/Final_Output/${samples}/", mode: 'copy'
	tag "${samples}"
	//label 'process_medium'
	input:
		tuple val (samples), file(finalBam), file (finalBamBai)
	output:
		tuple val (samples), file ("${samples}_longshot.txt")
	script:
	"""
	${params.longshot_sc} ${finalBam} ${params.genome} ${samples}.vcf
	#${params.vep_script_path} ${samples}.vcf ${samples}
	${params.vep_wrapper} ${samples} ${samples}.vcf ${samples}_longshot.txt
	"""
}

process Longshot_filter {
	tag "${samples}"
    label 'process_medium'
    input:
        tuple val (samples), file(longshot_vcf)
    output:
        tuple val (samples), file ("${samples}_longshot_filter.vcf")
    script:
    """
	bcftools view  -R ${bed_coverage} -i 'FORMAT/DP >= 10' ${samples}.vcf -Ov -o ${samples}_longshot_filter.vcf


	"""
}

process Nanocaller {
		tag "${samples}"
		//label 'process_medium'
	input:
		tuple val (samples), file(finalBam), file (finalBamBai)
	output:
		tuple val (samples), file ("${samples}_nanocaller.txt")
	script:
	"""
	#NanoCaller --bam ${finalBam} --ref ${params.genome} --sequencing ont --bed ${params.bed_file} --prefix ${samples} --sample ${samples}
	${params.nanocaller_sc} ${finalBam} ${params.genome} ont ${params.bed_file}  ${samples}
	gunzip ${samples}.vcf.gz
	${params.vep_wrapper} ${samples} ${samples}.vcf ${samples}_nanocaller.txt
	"""
}

process DeepPEPPER {
		tag "${samples}"
		//label 'process_medium'
	input:
		tuple val (samples), file(finalBam), file (finalBamBai)
	output:
		tuple val (samples), file ("*.vcf")
	script:
	"""
	singularity exec --bind /home/ ${params.DeepPEPPER} run_pepper_margin_deepvariant call_variant -b \${PWD}/${finalBam} -f ${params.genome} -o \${PWD} -p ${samples} -t 4 -r chr9:133738150-133756051 --ont_r10_q20
	"""
}

process ClairsTO {
		tag "${samples}"
		label 'process_high'
	input:
		tuple val (samples), file(finalBam), file (finalBamBai)
		path(reference_genome)
		path(ind_files)
		path(bed_coverage)
	output:
		tuple val (samples), file ("${samples}.vcf")
	script:
	"""
	#path=`realpath ./`

	#${params.clairto} ${params.genome} ${finalBam} ${params.bed_file} \${path} ont_r10_dorado_hac_4khz
	/opt/bin/run_clairs_to --tumor_bam_fn ${finalBam} --ref_fn ${reference_genome} -t ${task.cpus}  --platform ont_r10_dorado_hac_4khz -b ${params.bed_file} --output_dir ./   
	bcftools concat -a snv.vcf.gz indel.vcf.gz -o ${samples}.vcf
	#${params.vep_wrapper} ${samples} ${samples}.vcf ${samples}_ClarisTo.txt
	"""
}

process Clair_annotate {
		tag "${samples}"
		//label 'process_medium'
	input:
		tuple val (samples), file(vcf)
	output:
        tuple val (samples), file ("${samples}_ClarisTo.txt")
    script:
    """
    ${params.vep_wrapper} ${samples} ${samples}.vcf ${samples}_ClarisTo.txt
    """
}

process annovar {
        tag "${Sample}"
        label 'process_medium'
        publishDir "${params.outdir}/${samples}/", mode: 'copy', pattern: '*.csv'
        input:
                tuple val (Sample), path(Vcf)
                val(variant_caller)
        output:
                tuple val (Sample), ${Sample}_${variant_caller}.out.hg38_multianno.csv
        script:
        """
        convert2annovar.pl -format vcf4 ${Vcf} --outfile ${Sample}_${variant_caller}.avinput --withzyg --includeinfo

        table_annovar.pl ${Sample}_${variant_caller}.avinput --out ${Sample}_${variant_caller}.out  --remove --protocol refGene,cytoBand,cosmic84,popfreq_all_20150413,avsnp150,intervar_20180118,1000g2015aug_all,clinvar_20170905 --operation g,r,f,f,f,f,f,f --buildver hg38 --nastring '-1' --otherinfo --csvout --thread ${task.cpus} /databases/humandb --xreffile /databases/gene_fullxref.txt

        """
}

process CombineCallers {
	publishDir "results/${samples}/", mode: 'copy'
	tag "${samples}"
	//label 'process_low'
	input:
		tuple val(samples), file(CoverageFile), file (ClairsTO), file(NanoCallerOutput), file(LongshotOutput)
	output:
		tuple val(samples), file("${samples}_variants.xlsx")
	script:
	"""
	
	merge-tsv.py ${samples}_variants.xlsx ${CoverageFile} ${LongshotOutput} ${NanoCallerOutput} ${ClairsTO}
	
	"""
}

//workflow Nano_KDM_nodemux {
//	Channel
//		.fromPath(params.input)
//		.splitCsv(header:false)
//		.flatten()
//		.map{ it }
//		.set { samples_ch }
//	Align_nodemux(samples_ch)
//	Coverage(Align_nodemux.out)
//	Longshot(Align_nodemux.out)
//	Nanocaller(Align_nodemux.out)
//	//DeepPEPPER(Align_nodemux.out)
//	ClairsTO(Align_nodemux.out)
//	CombineCallers(Coverage.out.join(ClairsTO.out.join(Nanocaller.out.join(Longshot.out))))
//}


//workflow.onComplete {
//	log.info ( workflow.success ? "\n\nDone! Output in the 'Final_Output' directory \n" : "Oops .. something went wrong" )
//}

