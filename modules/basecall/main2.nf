process BASECALL {
        publishDir "${baseDir}/sequences/fastqs", mode: 'copy'
	input:
                path (pod5_dir)
                path (dorado_models_dir)
        output:
                path ("fastqs")
        script:
        """
	echo "Node:"
        hostname

        echo "GPU devices:"
        ls -l /dev/nvidia*

        echo "Driver status:"
        nvidia-smi
	
	dorado basecaller ${dorado_models_dir}/dna_r10.4.1_e8.2_400bps_hac@v5.0.0  ${pod5_dir}/ \
	  --kit-name SQK-NBD114-24 \
	  --barcode-both-ends \
	  --no-trim \
	  -x cuda:all \
	  > basecalls.bam

	dorado demux -o fastqs --kit-name SQK-NBD114-24 -t 64 --barcode-both-ends --emit-fastq basecalls.bam


	"""
}

