#!/usr/bin/env bash

FASTQ_DIR=$1
SAMPLESHEET=$2

while read SAMPLEBARCODE
do
    BARCODE=$(echo $SAMPLEBARCODE | grep -o "barcode[0-9]\+")

    fq=$(find ${FASTQ_DIR} -name "*${BARCODE}.fastq" | head -n 1)

    if [[ -f "$fq" ]]; then
        gzip -c "$fq" > ${FASTQ_DIR}/${SAMPLEBARCODE}.fastq.gz
    else
        echo "WARNING: No FASTQ found for $SAMPLEBARCODE"
    fi

done < $SAMPLESHEET

# safety: ensure unclassified reads never get renamed
rm -f ${FASTQ_DIR}/unclassified*.fastq.gz 2>/dev/null
