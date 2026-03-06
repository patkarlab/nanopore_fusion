#!/bin/bash
#PBS -q medium
#PBS -N nextflow
#PBS -l select=1:ncpus=32:mem=64gb
#PBS -l walltime=48:00:00
#PBS -j oe
cd $PBS_O_WORKDIR

# ---- Conda initialization
source ~/miniconda3/etc/profile.d/conda.sh
conda activate base

module load samtools/1.21
module load bedtools/2.27.1
module load minimap/2.28
module load python/3.11.0
module load matplotlib/3.10.8

export NXF_VER=25.10.2


# ---- Run pipeline
  
nextflow run /home/patkarlab/pipelines/nanopore_fusion/main.nf \
	-entry ACTREC_NANOPORE_FUSION \
	--input /home/patkarlab/pipelines/nanopore_fusion/samplesheet.csv \
	--outdir /home/patkarlab/pipelines/nanopore_fusion/results  \
	--bed_file /home/patkarlab/pipelines/nanopore_fusion/assets/RADICALv3_hg38_sortd.bed \
	-profile pbs,docker \
	-ansi-log false -resume 
