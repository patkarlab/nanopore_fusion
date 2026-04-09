process BASECALL {

    publishDir "/home/patkarlab/pipelines/nfInputs/", mode: 'copy', pattern: "fastqs/*.fastq.gz"
    

    input:
        path (pod5_dir)
        path (dorado_models_dir)
        path (samplesheet_file)

    output:
        path ("fastqs/*.fastq.gz")

    script:
    """
    echo "Node:"
    hostname

    echo "GPU devices:"
    ls -l /dev/nvidia* || true

    echo "Driver status:"
    nvidia-smi || true

    mkdir -p fastqs

    dorado basecaller ${dorado_models_dir}/dna_r10.4.1_e8.2_400bps_hac@v5.0.0 ${pod5_dir}/ \\
      --kit-name SQK-NBD114-24 \\
      --barcode-both-ends \\
      --no-trim \\
      -x cuda:all \\
      > basecalls.bam

    dorado demux -o fastqs \\
      --kit-name SQK-NBD114-24 \\
      -t 64 \\
      --barcode-both-ends \\
      --emit-fastq \\
      basecalls.bam


    rename_demux_fastqs.sh fastqs ${samplesheet_file}
    """
}
