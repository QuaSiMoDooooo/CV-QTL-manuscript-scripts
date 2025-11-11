#! /bin/bash
# FileName     : enrich_UCSC_variant
# Author       : EastsunW eastsunw@foxmail.com
# Create at    : 2024-12-26 22:01
# Last Modified: 2025-01-05 17:30
# Modified By  : EastsunW
# -------------
# Description  : 
# -------------

set -euo pipefail
source /home/wangdy/miniforge3/bin/activate gat

outdir=results/epi_enrichment/UCSC/variant
mkdir -p $outdir

for variant_type in SV MNV InDel SNP; do
    if [ "$variant_type" == "SV" ]; then
        freq_types=("all" "common" "rare" "INS" "DEL" "INV" "DUP")
    else
        freq_types=("all" "common" "rare")
    fi
    for freq_type in "${freq_types[@]}"; do
        gat-run.py \
            -t 5 \
            -n 1000 -F \
            --bucket-size=500 \
            --nbuckets=200000 \
            --ignore-segment-tracks \
            -s data/Variants/All/${variant_type}.${freq_type}.bed \
            -a data/Epigenomics/UCSC/UCSC.annotations.converted.bed \
            -w data/Common/hg38.chrlen.bed \
            -L ${outdir}/UCSC_${variant_type}_${freq_type}.log \
            -S ${outdir}/UCSC_${variant_type}_${freq_type}.results.txt
        echo "${variant_type} ${freq_type} done"
    done
done

echo "UCSC variant all done"
