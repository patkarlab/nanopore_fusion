process JAFFAL {

    label 'process_medium'
    publishDir "results/${samples}", mode:'copy'
    tag "${samples}"


    input:
        tuple val(samples), path(reads)

    output:
    tuple val(samples), path("${samples}_jaffal.tsv"),    emit: jaffal_tsv
    tuple val(samples), path("${samples}_jaffal.fasta")  


    script:
    """

    bpipe run /opt/JAFFA/JAFFAL.groovy ${reads} 

    mv jaffa_results.fasta ${samples}_jaffal.fasta
    
    awk -F',' 'BEGIN {OFS="\t"} {print \$2,\$3,\$4,\$5,\$6,\$7,\$8,\$11,\$16,\$17,\$18}' jaffa_results.csv > ${samples}_jaffal.tsv

    """
}
