#!/usr/bin/env bash

FASTQ_DIR=$1
SAMPLESHEET=$2

while read SAMPLEBARCODE
do
    BARCODE=$(echo $SAMPLEBARCODE | grep -o "barcode[0-9]\+")

    fq=$(find ${FASTQ_DIR} -name "*_${BARCODE}.fastq" | head -n 1)

    if [[ -f "$fq" ]]; then
        gzip -c "$fq" > ${FASTQ_DIR}/${SAMPLEBARCODE}.fastq.gz
    else
        echo "WARNING: No FASTQ found for $SAMPLEBARCODE"
    fi

done < $SAMPLESHEET
