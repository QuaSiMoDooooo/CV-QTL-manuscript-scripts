#! /bin/bash
# FileName     : enrich_chromHMM_variant
# Author       : EastsunW eastsunw@foxmail.com
# Create at    : 2024-12-26 10:38
# Last Modified: 2025-01-05 12:08
# Modified By  : EastsunW
# -------------
# Description  :
# -------------

set -euo pipefail
source /home/wangdy/miniforge3/bin/activate gat

outdir=results/epi_enrichment/chromHMM/variant
mkdir -p $outdir

for chromhmmfile in data/chromHMM/E*_15_coreMarks_hg38lift_converted.bed; do
    filename=$(basename "$chromhmmfile")
    annot_type=$(echo "$filename" | grep -oP "^E\\d+")
    for variant_type in SNP; do
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
                -a $chromhmmfile \
                -w data/Common/hg38.chrlen.bed \
                -L ${outdir}/${annot_type}_${variant_type}_${freq_type}.log \
                -S ${outdir}/${annot_type}_${variant_type}_${freq_type}.results.txt
            echo "${annot_type} ${variant_type} ${freq_type} done"
        done
    done
done
echo "chromHMM Variant all done"
