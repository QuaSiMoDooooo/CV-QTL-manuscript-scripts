#! /bin/bash
# FileName     : enrich_HMdb_variant
# Author       : EastsunW eastsunw@foxmail.com
# Create at    : 2024-12-26 21:36
# Last Modified: 2025-01-05 17:38
# Modified By  : EastsunW
# -------------
# Description  :
# -------------

set -euo pipefail
source /home/wangdy/miniforge3/bin/activate gat

outdir=results/epi_enrichment/HMdb/qtl
mkdir -p $outdir

for variant_type in SV MNV InDel SNP; do
    for QTL_type in eQTL apaQTL sQTL; do
        gat-run.py \
            -t 5 \
            -n 1000 -F \
            --bucket-size=500 \
            --nbuckets=200000 \
            --ignore-segment-tracks \
            -s data/Variants/QTL/${variant_type}_${QTL_type}.bed \
            -a data/Epigenomics/HMdb/HMdb_all.converted.bed \
            -w data/Variants/All/${variant_type}.bg.bed \
            -L ${outdir}/HMdb.${variant_type}-${QTL_type}.log \
            -S ${outdir}/HMdb.${variant_type}-${QTL_type}.results.txt
    echo "${variant_type}-${QTL_type} done"
    done
done
echo "HMdb QTL All done"
