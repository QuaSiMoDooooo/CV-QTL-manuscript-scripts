#! /bin/bash
# FileName     : stat_qtls
# Author       : EastsunW eastsunw@foxmail.com
# Create at    : 2025-05-26 16:55
# Last Modified: 2025-05-26 16:58
# Modified By  : EastsunW
# -------------
# Description  : 
# -------------

set -euo pipefail
#!/bin/bash

# 变量定义
variants=("SNP" "InDel" "MNV" "SV")
phenotypes=("expression" "APA" "splicing" "methylation")
cis_types=("cis" "trans")

for cis in "${cis_types[@]}"; do
    for v in "${variants[@]}"; do
        for p in "${phenotypes[@]}"; do
            file="/home/wangdy/Projects/Weibin/Downstreams/phenotype_qtl/results/${v}-${p}/QTL_results/${cis}.filtered.txt.gz"
            if [[ -f "$file" ]]; then
                count=$(zcat "$file" | wc -l)
                count=$((count - 1))
                echo "${cis},${v},${p},$count"
            else
                echo "${cis},${v},${p},file_not_found"
            fi
        done
    done
done
