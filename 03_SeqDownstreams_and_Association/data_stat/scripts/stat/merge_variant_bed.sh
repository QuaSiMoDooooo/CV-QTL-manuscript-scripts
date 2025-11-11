#! /bin/bash
# FileName     : merge_bed
# Author       : EastsunW eastsunw@foxmail.com
# Create at    : 2024-12-05 19:14
# Last Modified: 2024-12-05 19:14
# Modified By  : EastsunW
# -------------
# Description  : 
# -------------

set -euo pipefail
/home/wangdy/miniforge3/envs/test/bin/bedtools sort -i results/bed/InDel.all.bed | /home/wangdy/miniforge3/envs/test/bin/bedtools merge > results/bed/InDel.all.merged.bed
/home/wangdy/miniforge3/envs/test/bin/bedtools sort -i results/bed/InDel.common.bed | /home/wangdy/miniforge3/envs/test/bin/bedtools merge > results/bed/InDel.common.merged.bed
/home/wangdy/miniforge3/envs/test/bin/bedtools sort -i results/bed/InDel.rare.bed | /home/wangdy/miniforge3/envs/test/bin/bedtools merge > results/bed/InDel.rare.merged.bed
/home/wangdy/miniforge3/envs/test/bin/bedtools sort -i results/bed/MNV.all.bed | /home/wangdy/miniforge3/envs/test/bin/bedtools merge > results/bed/MNV.all.merged.bed
/home/wangdy/miniforge3/envs/test/bin/bedtools sort -i results/bed/MNV.common.bed | /home/wangdy/miniforge3/envs/test/bin/bedtools merge > results/bed/MNV.common.merged.bed
/home/wangdy/miniforge3/envs/test/bin/bedtools sort -i results/bed/MNV.rare.bed | /home/wangdy/miniforge3/envs/test/bin/bedtools merge > results/bed/MNV.rare.merged.bed
/home/wangdy/miniforge3/envs/test/bin/bedtools sort -i results/bed/SNP.all.bed | /home/wangdy/miniforge3/envs/test/bin/bedtools merge > results/bed/SNP.all.merged.bed
/home/wangdy/miniforge3/envs/test/bin/bedtools sort -i results/bed/SNP.common.bed | /home/wangdy/miniforge3/envs/test/bin/bedtools merge > results/bed/SNP.common.merged.bed
/home/wangdy/miniforge3/envs/test/bin/bedtools sort -i results/bed/SNP.rare.bed | /home/wangdy/miniforge3/envs/test/bin/bedtools merge > results/bed/SNP.rare.merged.bed
/home/wangdy/miniforge3/envs/test/bin/bedtools sort -i results/bed/SV.all.bed | /home/wangdy/miniforge3/envs/test/bin/bedtools merge > results/bed/SV.all.merged.bed
/home/wangdy/miniforge3/envs/test/bin/bedtools sort -i results/bed/SV.common.bed | /home/wangdy/miniforge3/envs/test/bin/bedtools merge > results/bed/SV.common.merged.bed
/home/wangdy/miniforge3/envs/test/bin/bedtools sort -i results/bed/SV.rare.bed | /home/wangdy/miniforge3/envs/test/bin/bedtools merge > results/bed/SV.rare.merged.bed
