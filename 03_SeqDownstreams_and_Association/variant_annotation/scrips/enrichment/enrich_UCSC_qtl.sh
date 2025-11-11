#! /bin/bash
# FileName     : enrich_UCSC_variant
# Author       : EastsunW eastsunw@foxmail.com
# Create at    : 2024-12-26 21:36
# Last Modified: 2025-03-12 21:55
# Modified By  : EastsunW
# -------------
# Description  :
# -------------

set -euo pipefail
source /home/wangdy/miniforge3/bin/activate gat

outdir=results/epi_enrichment/UCSC/qtl
mkdir -p $outdir

for variant_type in SV MNV InDel SNP; do
    for QTL_type in meQTL; do
        gat-run.py \
            -t 5 \
            -n 1000 -F \
            --bucket-size=500 \
            --nbuckets=200000 \
            --ignore-segment-tracks \
            -s data/Variants/QTL/${variant_type}_${QTL_type}.bed \
            -a data/Epigenomics/UCSC/UCSC.annotations.converted.bed \
            -w data/Variants/All/${variant_type}.bg.bed \
            -L ${outdir}/UCSC_${variant_type}-${QTL_type}.log \
            -S ${outdir}/UCSC_${variant_type}-${QTL_type}.results.txt
    echo "${variant_type}-${QTL_type} done"
    done
done

echo "UCSC QTL all done"
