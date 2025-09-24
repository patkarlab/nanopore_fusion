# actrec/nanopore_fusion

## Introduction

**actrec/nanopore_fusion** is a bioinformatics pipeline for calling fusions from long-read RNA sequencing data. The outputs from fusion callers(ctat-lr-fusion, genion, longgf and fusionseeker) are collected in single excel multisheet excel file, along with the coverage information. It also collates the fusion caller outputs into a single excel sheet, which can be used asinput for fusviz, which is a web-based RNA/DNA fusion visualization tool. Other outputs include a fusion visualization report in html format from ctat-lr-fusion, a read lenth histogram for the fastq inputs and sorted,indexed bam files.


## Usage

The follwing parameters need to be modified in the `params` section of the `nextflow.config` - 

- *genome* = The Reference genome for this pipeline can be downloaded from https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta

- *gtf* = The GTF file for Genion and LongGF can be downloaded from -
ftp://ftp.ensembl.org/pub/release-105/gtf/homo_sapiens/Homo_sapiens.GRCh38.105.gtf.gz

- *gtf_2* = The GTF file for fusionseeker can be downloaded from -
ftp://ftp.ensembl.org/pub/release-104/gtf/homo_sapiens/Homo_sapiens.GRCh38.104.gtf.gz

Please check if the chromosome names are matching in reference genome and the gtf files.

- *cdna* = Genion additionally needs 'cdna.selfalign.tsv' and 'genomicSuperDups.txt' as input. These files need to be generated as per instructions on Genion github -  https://github.com/vpc-ccg/genion?tab=readme-ov-file#input

- *genome_lib* = CTAT-LR-fusion requires CTAT genome lib, which can be downloaded from -
https://data.broadinstitute.org/Trinity/CTAT_RESOURCE_LIB/GRCh38_gencode_v22_CTAT_lib_Mar012021.source.tar.gz


## Running the pipeline

1. Transfer the `fastq.gz` files to the `sequences/` folder.

2. The samplesheet is `samplesheet.csv`. The sample_ids, without the file extension, should be mentioned in samplesheet in the following format-
sample1
sample2
sample3
Please check for empty lines in the samplesheet before running the pipeline.


3. The pipeline can be run by running the command-

```bash
nextflow run . -entry ACTREC_NANOPORE_FUSION --input samplesheet.csv --outdir results  -profile docker -ansi-log false  -resume  -bg
```

## Output
The outputs are saved in `results/` folder.
