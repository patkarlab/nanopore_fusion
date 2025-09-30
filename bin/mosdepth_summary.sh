#!/usr/bin/env bash
# mosdepth_summary.sh
# Convert mosdepth thresholds+regions into coverage summary table

set -euo pipefail

if [ $# -ne 2 ]; then
    echo "Usage: mosdepth_summary.sh <prefix.thresholds.bed.gz> <prefix.regions.bed.gz>" >&2
    exit 1
fi

THR=$1
REG=$2
OUT=$(basename "$THR" .thresholds.bed.gz).cov_mosdepth.tsv

zcat "$THR" "$REG" | \
awk 'BEGIN{
    OFS="\t"
    print "chrom","start","end","region","region_length",
          "cov>=1X(%)","cov>=10X(%)","cov>=30X(%)","cov>=100X(%)","cov>=200X(%)","mean_coverage"
}
FNR==NR && !/^#/{
    chrom=$1; start=$2; end=$3; region=$4;
    one=$5; ten=$6; thirty=$7; hundred=$8; twohund=$9;
    len=end-start;

    cov1=(one/len*100);
    cov10=(ten/len*100);
    cov30=(thirty/len*100);
    cov100=(hundred/len*100);
    cov200=(twohund/len*100);

    key=chrom":"start":"end;
    thr[key]=chrom OFS start OFS end OFS region OFS len OFS \
              sprintf("%.3f",cov1) OFS sprintf("%.3f",cov10) OFS sprintf("%.3f",cov30) OFS sprintf("%.3f",cov100) OFS sprintf("%.3f",cov200);
    next
}
!/^#/{
    chrom=$1; start=$2; end=$3; region=$4; mean=$5;
    key=chrom":"start":"end;
    if (key in thr) {
        print thr[key], mean;
    }
}' "$THR" "$REG" > "$OUT"

