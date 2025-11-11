#! /bin/bash
# FileName     : check_overlap
# Author       : EastsunW eastsunw@foxmail.com
# Create at    : 2025-01-06 16:09
# Last Modified: 2025-01-06 16:33
# Modified By  : EastsunW
# -------------
# Description  : 
# -------------

set -euo pipefail

# 检查所有变异的plof：和至少一个CDS重合
bedtools intersect \
    -a data/Variants/All/SNP.all.bed \
    -b data/Epigenomics/UCSC/UCSC_coding.bed \
    -wa -wb \
    > results/epi_enrichment/plof/SNP.plof.bed
bedtools intersect \
    -a data/Variants/All/InDel.all.bed \
    -b data/Epigenomics/UCSC/UCSC_coding.bed \
    -wa -wb \
    > results/epi_enrichment/plof/InDel.plof.bed
bedtools intersect \
    -a data/Variants/All/MNV.all.bed \
    -b data/Epigenomics/UCSC/UCSC_coding.bed \
    -wa -wb \
    > results/epi_enrichment/plof/MNV.plof.bed
bedtools intersect \
    -a data/Variants/All/SV.all.bed \
    -b data/Epigenomics/UCSC/UCSC_coding.bed \
    -wa -wb \
    > results/epi_enrichment/plof/SV.plof.bed

# 检查完全覆盖基因区域的变异
# bedtools intersect \
#     -a data/Variants/All/SNP.all.bed \
#     -b data/Epigenomics/UCSC/UCSC_gene.bed \
#     -wa -wb -F 1.0 \
#     > results/epi_enrichment/plof/SNP.whole_gene.bed
# bedtools intersect \
#     -a data/Variants/All/InDel.all.bed \
#     -b data/Epigenomics/UCSC/UCSC_gene.bed \
#     -wa -wb -F 1.0 \
#     > results/epi_enrichment/plof/InDel.whole_gene.bed
# bedtools intersect \
#     -a data/Variants/All/MNV.all.bed \
#     -b data/Epigenomics/UCSC/UCSC_gene.bed \
#     -wa -wb -F 1.0 \
#     > results/epi_enrichment/plof/MNV.whole_gene.bed
bedtools intersect \
    -a data/Variants/All/SV.all.bed \
    -b data/Epigenomics/UCSC/UCSC_gene.bed \
    -wa -wb -F 1.0 \
    > results/epi_enrichment/plof/SV.whole_gene.bed
