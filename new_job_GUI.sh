#!/bin/bash
#PBS -q medium
#PBS -N nanopore_fusion
#PBS -l select=1:ncpus=32:mem=64gb
#PBS -l walltime=48:00:00
#PBS -j oe

cd $PBS_O_WORKDIR || exit 1

# initialization of conda env
source ~/miniconda3/etc/profile.d/conda.sh
conda activate base

# Load modules
module load samtools/1.21
module load samtools/1.21
module load bedtools/2.27.1
module load python/3.11.0
module load matplotlib/3.10.8

export NXF_VER=25.10.2

# detect latest GUI upload folder
SEQ_BASE="/home/patkarlab/pipelines/nanopore_fusion/sequences"

LATEST_SEQ=$(find "$SEQ_BASE" -maxdepth 1 -type d -name "ShellScript_*" \
    -printf "%T@ %p\n" 2>/dev/null \
    | sort -nr \
    | head -1 \
    | cut -d' ' -f2-)

if [ -z "$LATEST_SEQ" ]; then
    echo "ERROR: No ShellScript_* directory found in $SEQ_BASE"
    exit 1
fi

echo "Using latest uploaded sequence folder:"
echo "$LATEST_SEQ"

# run pipeline (override params.sequences)
nextflow run /home/patkarlab/pipelines/nanopore_fusion/main.nf \
        -entry ACTREC_NANOPORE_FUSION \
        --input /home/patkarlab/pipelines/nanopore_fusion/samplesheet.csv \
	--sequences "$LATEST_SEQ"
        --outdir /home/patkarlab/pipelines/nanopore_fusion/results  \
        --bed_file /home/patkarlab/pipelines/nanopore_fusion/assets/RADICALv3_hg38_sortd.bed \
        -profile pbs,docker \
        -ansi-log false -resume


