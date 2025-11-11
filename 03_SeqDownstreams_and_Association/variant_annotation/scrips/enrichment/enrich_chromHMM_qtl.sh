#! /bin/bash
# FileName     : enrich_chromHMM_variant
# Author       : EastsunW eastsunw@foxmail.com
# Create at    : 2024-12-26 10:38
# Last Modified: 2025-03-05 20:50
# Modified By  : EastsunW
# -------------
# Description  :
# -------------

set -euo pipefail
source /home/wangdy/miniforge3/bin/activate gat

outdir=results/epi_enrichment/chromHMM/qtl
mkdir -p $outdir

for chromhmmfile in $(ls data/chromHMM/E*_15_coreMarks_hg38lift_converted.bed); do
    filename=$(basename "$chromhmmfile")
    annot_type=$(echo "$filename" | grep -oP '^E\d+')
    for variant_type in SNP; do
        # for QTL_type in meQTL; do
        for QTL_type in eQTL apaQTL sQTL meQTL; do
            gat-run.py \
                -t 5 \
                -n 1000 -F \
                --bucket-size=500 \
                --nbuckets=200000 \
                --ignore-segment-tracks \
                -s data/Variants/QTL/${variant_type}_${QTL_type}.bed \
                -a ${chromhmmfile} \
                -w data/Variants/All/${variant_type}.bg.bed \
                -L ${outdir}/chromHMM_${annot_type}.${variant_type}-${QTL_type}.log \
                -S ${outdir}/chromHMM_${annot_type}.${variant_type}-${QTL_type}.results.txt
            echo "${annot_type} ${variant_type} ${QTL_type} done"
        done
    done
done
echo "chromHMM QTL All done"
