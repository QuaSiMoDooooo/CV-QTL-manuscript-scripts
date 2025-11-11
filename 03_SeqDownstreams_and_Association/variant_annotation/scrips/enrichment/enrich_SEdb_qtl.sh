#! /bin/bash
# FileName     : enrich_SEdb_variant
# Author       : EastsunW eastsunw@foxmail.com
# Create at    : 2024-12-26 21:36
# Last Modified: 2024-12-30 17:30
# Modified By  : EastsunW
# -------------
# Description  :
# -------------

set -euo pipefail
source /home/wangdy/miniforge3/bin/activate gat

outdir=results/epi_enrichment/SEdb/qtl
mkdir -p $outdir

for annotation in $(ls data/Epigenomics/SEdb/SEdb_*.bed); do
    filename=$(basename "$annotation")
    annot_type=$(echo "$filename" | grep -oP '(?<=SEdb_)([^_]+)(?=\.bed)')
    for variant_type in SV MNV InDel; do
        for QTL_type in eQTL apaQTL sQTL; do
            gat-run.py \
                -t 5 \
                -n 1000 -F \
                --bucket-size=500 \
                --nbuckets=200000 \
                --ignore-segment-tracks \
                -s data/Variants/QTL/${variant_type}_${QTL_type}.bed \
                -a ${annotation} \
                -w data/Variants/All/${variant_type}.bg.bed \
                -L ${outdir}/SEdb_${annot_type}.${variant_type}-${QTL_type}.log \
                -S ${outdir}/SEdb_${annot_type}.${variant_type}-${QTL_type}.results.txt
        echo "${annot_type}.${variant_type}-${QTL_type} done"
        done
    done
done
echo "SEdb QTL All done"
