#!/bin/bash
#PBS -N dorado_basecalling
#PBS -q a40
#PBS -l select=1:ncpus=16:mem=32gb:ngpus=1
#PBS -l walltime=24:00:00
#PBS -l container_image=ln1:5000/dorado/dorado:local
#PBS -j oe
#PBS -o dorado_basecalling.log

cd $PBS_O_WORKDIR

echo "Running on $(hostname)"
nvidia-smi

#dorado basecaller /home/patkarlab/pipelines/dorado_models/dna_r10.4.1_e8.2_400bps_hac@v5.0.0 \
#  /home/patkarlab/pipelines/nanopore_fusion/sequences/pod5 \
#  --kit-name SQK-NBD114-24 \
#  --barcode-both-ends \
#  --no-trim \
#  -x cuda:0 \
#  > /home/patkarlab/pipelines/nanopore_fusion/sequences/basecalls.bam

dorado demux -o fastqs_040326 --kit-name SQK-NBD114-24 -t 64 --barcode-both-ends --emit-fastq /home/patkarlab/pipelines/nanopore_fusion/sequences/basecalls.bam


echo "Finished at $(date)"
