#!/usr/bin/bash

export JAVA_HOME=/home/programs/nextflow/jdk-17.0.9
export PATH="$JAVA_HOME/bin:$PATH"
export NXF_VER=23.10.1


nextflow run . -entry ACTREC_NANOPORE_FUSION --input samplesheet.csv --outdir results  -profile docker -ansi-log false -resume -bg --bed_file /home/diagnostics/pipelines/nf-core/rnafusion/bedfiles/BALLlymphoid_fusion02062022_hg38.bed
