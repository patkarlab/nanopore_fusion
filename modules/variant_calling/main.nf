#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process Longshot {
	//publishDir "results/${samples}", mode:'copy'
	tag "${samples}"
	//label 'process_medium'
	input:
		tuple val (samples), file(finalBam), file (finalBamBai)
		path(reference_genome)
        path(ind_files)
	output:
		tuple val (samples), file ("${samples}_longshot_prefilter.vcf")
	script:
	"""
	
	longshot -b ${finalBam} -f ${reference_genome} -o ${samples}_longshot_prefilter.vcf
	"""
}

process Nanocaller {
	//publishDir "results/${samples}", mode:'copy'
	tag "${samples}"
	//label 'process_medium'
	input:
		tuple val (samples), file(finalBam), file (finalBamBai)
		path(reference_genome)
        path(ind_files)
        path(bed_coverage)
	output:
		tuple val (samples), file ("${samples}_nanocaller.vcf")
	script:
	"""
	NanoCaller --bam ${finalBam} --ref ${reference_genome} --sequencing ont --bed ${bed_coverage} --snp_model ONT-HG002_r10.3 --prefix ${samples} --sample ${samples}

	gunzip ${samples}.vcf.gz
	mv ${samples}.vcf ${samples}_nanocaller.vcf
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

process Longcallr {
	//publishDir "results/${samples}", mode:'copy'
	tag "${samples}"
	label 'process_high'
	input:
		tuple val (samples), file(finalBam), file (finalBamBai)
		path(reference_genome)
		path(ind_files)
	output:
		tuple val (samples), file ("${samples}_longcallr_prefilter.vcf")
	script:
	"""
	longcallR --bam-path ${finalBam} --ref-path ${reference_genome} --threads ${task.cpus} --platform ont \
	--min-allele-freq 0.05 --min-mapq 5 --min-qual-for-candidate 50 --somatic-allele-frac-cutoff 0.01 --somatic-allele-cnt-cutoff 2 \
	--output ${samples}_longcallr_prefilter 
	"""
}

process Clair3rna {
	//publishDir "results/${samples}", mode:'copy'
	tag "${samples}"
	label 'process_high'
	input:
		tuple val (samples), file(finalBam), file (finalBamBai)
		path(reference_genome)
		path(ind_files)
		path(bed_coverage)
	output:
		tuple val (samples), file ("${samples}_clair3rna.vcf")
	script:
	"""
	/opt/bin/run_clair3_rna \
		--bam_fn ${finalBam} --ref_fn ${reference_genome} -b ${bed_coverage} \
		--threads ${task.cpus} --platform ont_r10_dorado_cdna --snp_min_af 0.01 --min_mq 3 --qual 2 \
		--enable_variant_calling_at_sequence_head_and_tail --output_dir ./
	gunzip -c output.vcf.gz > ${samples}_clair3rna.vcf
	"""
}

process VarDict {
        tag "${samples}"
        //publishDir "${params.outdir}/${samples}/", mode : 'copy'
        label 'process_medium'
        input:
                tuple val (samples), file(finalBam), file(finalBamBai)
				path(reference_genome)
                path(ind_files)
                file(bed_coverage)
        output:
                tuple val(Sample), file ("${samples}_vardict.vcf")
        script:
        """
        VarDict -G ${reference_genome} -f 0.01 -N ${samples} -b ${finalBam} -q 10 -Q 0  -c 1 -S 2 -E 3 -g 4  ${bed_coverage}| sed '1d' | teststrandbias.R | var2vcf_valid.pl -N ${samples} -E -f 0.01 > ${samples}_vardict.vcf
        """
}

process ClairsTO {
	//publishDir "results/${samples}", mode:'copy'
	tag "${samples}"
	label 'process_high'
	input:
		tuple val (samples), file(finalBam), file (finalBamBai)
		path(reference_genome)
		path(ind_files)
		path(bed_coverage)
	output:
		tuple val (samples), file ("${samples}_clairsto.vcf")
	script:
	"""
	/opt/bin/run_clairs_to --tumor_bam_fn ${finalBam} --ref_fn ${reference_genome} -t ${task.cpus}  --platform ont_r10_dorado_hac_4khz -b ${params.bed_file} --output_dir ./   
	bcftools concat -a snv.vcf.gz indel.vcf.gz -o ${samples}_clairsto.vcf
	"""
}

process vcf_filter {
    //publishDir "results/${samples}", mode:'copy'
    tag "${samples}"
    label 'process_medium'
    input:
        tuple val (samples), file(Vcf)
    path(bed_coverage)
    val(variant_caller)
    output:
        tuple val (samples), path ("${samples}_${variant_caller}.vcf")
    script:
    """
    bcftools sort ${samples}_${variant_caller}_prefilter.vcf -Ov -o ${samples}_${variant_caller}_prefilter_sorted.vcf
	bgzip -c ${samples}_${variant_caller}_prefilter_sorted.vcf > ${samples}_${variant_caller}_prefilter_sorted.vcf.gz
    tabix -p vcf ${samples}_${variant_caller}_prefilter_sorted.vcf.gz
    bcftools view  -R ${bed_coverage} ${samples}_${variant_caller}_prefilter_sorted.vcf.gz -Ov -o ${samples}_${variant_caller}.vcf
    """
}

process ANNOVAR {
	tag "${samples}"
	label 'process_medium'
	//publishDir "${params.outdir}/${samples}/", mode: 'copy', pattern: '*.csv'
	input:
		tuple val (samples), path(Vcf)
		val(variant_caller)
	output:
		tuple val (samples), path("${samples}_${variant_caller}.out.hg38_multianno.csv")
	script:
	"""
	convert2annovar.pl -format vcf4 ${Vcf} --outfile ${samples}_${variant_caller}.avinput --withzyg --includeinfo

	table_annovar.pl ${samples}_${variant_caller}.avinput --out ${samples}_${variant_caller}.out \
	--remove --protocol refGene,cytoBand,avsnp151,intervar_20180118,1000g2015aug_all,cosmic70,clinvar_20250721,gnomad211_exome --operation g,r,f,f,f,f,f,f \
	--buildver hg38 --nastring '-1' --otherinfo --csvout --thread ${task.cpus} /databases/humandb -xreffile /databases/gene_fullxref.txt

	touch ${samples}_${variant_caller}.out.hg38_multianno.csv
	"""
}

process Combine_variants {
	publishDir "results/${samples}/", mode: 'copy'
	tag "${samples}"
	label 'process_medium'
	input:
		tuple val(samples), file (ClairsTO_anno), file(NanoCaller_anno), file(Longshot_anno), file(clair3_rna_anno), file(longcallr_anno)
	output:
		tuple val(samples), file("${samples}_variants_anno.xlsx")
		tuple val(samples), file( "${samples}_support_distribution.png")
	script:
	"""

	combine_anno_vcf.py -s ${samples}  -l ${Longshot_anno} -n ${NanoCaller_anno} -c ${ClairsTO_anno} -r ${clair3_rna_anno} \
  -g ${longcallr_anno}  -o ${samples}_variants_anno.xlsx
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
