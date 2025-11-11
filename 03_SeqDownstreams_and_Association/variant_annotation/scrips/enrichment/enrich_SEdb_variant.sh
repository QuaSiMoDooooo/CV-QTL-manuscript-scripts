#! /bin/bash
# FileName     : enrich_SEdb_variant
# Author       : EastsunW eastsunw@foxmail.com
# Create at    : 2024-12-26 22:01
# Last Modified: 2025-01-01 23:26
# Modified By  : EastsunW
# -------------
# Description  : 
# -------------

set -euo pipefail
source /home/wangdy/miniforge3/bin/activate gat

outdir=results/epi_enrichment/SEdb/variant
mkdir -p $outdir

for annotation in $(ls data/Epigenomics/SEdb/SEdb_*.bed); do
    filename=$(basename "$annotation")
    annot_type=$(echo "$filename" | grep -oP '(?<=SEdb_)([^_]+)(?=\.bed)')
    for variant_type in SV MNV InDel; do
        if [ "$variant_type" == "SV" ]; then
            freq_types=("all" "common" "rare" "INS" "DEL" "INV" "DUP")
        else
            freq_types=("all" "common" "rare")
        fi
        for freq_type in "${freq_types[@]}"; do
            gat-run.py \
                -t 5 \
                -n 1000 \
                --bucket-size=500 \
                --nbuckets=200000 \
                --ignore-segment-tracks \
                -s data/Variants/All/${variant_type}.${freq_type}.bed \
                -a ${annotation} \
                -w data/Common/hg38.chrlen.bed \
                -L ${outdir}/SEdb_${annot_type}_${variant_type}_${freq_type}.log \
                -S ${outdir}/SEdb_${annot_type}_${variant_type}_${freq_type}.results.txt
            echo "${annot_type} ${variant_type} ${freq_type} done"
        done
    done
done
echo "SEdb variant all done"
