#! /bin/bash
# FileName     : stat_QTL_count
# Author       : EastsunW eastsunw@foxmail.com
# Create at    : 2024-12-10 23:40
# Last Modified: 2025-03-04 16:37
# Modified By  : EastsunW
# -------------
# Description  : 
# -------------

set -euo pipefail

pushd /home/wangdy/Projects/Weibin/Downstreams/data_stat/data/QTL_results
{
zcat SNP_eQTL/cis.filtered.txt.gz | wc -l
zcat SNP_eQTL/trans.filtered.txt.gz | wc -l
zcat SNP_apaQTL/cis.filtered.txt.gz | wc -l
zcat SNP_apaQTL/trans.filtered.txt.gz | wc -l
zcat SNP_sQTL/cis.filtered.txt.gz | wc -l
zcat SNP_sQTL/trans.filtered.txt.gz | wc -l

zcat InDel_eQTL/cis.filtered.txt.gz | wc -l
zcat InDel_eQTL/trans.filtered.txt.gz | wc -l
zcat InDel_apaQTL/cis.filtered.txt.gz | wc -l
zcat InDel_apaQTL/trans.filtered.txt.gz | wc -l
zcat InDel_sQTL/cis.filtered.txt.gz | wc -l
zcat InDel_sQTL/trans.filtered.txt.gz | wc -l

zcat MNV_eQTL/cis.filtered.txt.gz | wc -l
zcat MNV_eQTL/trans.filtered.txt.gz | wc -l
zcat MNV_apaQTL/cis.filtered.txt.gz | wc -l
zcat MNV_apaQTL/trans.filtered.txt.gz | wc -l
zcat MNV_sQTL/cis.filtered.txt.gz | wc -l
zcat MNV_sQTL/trans.filtered.txt.gz | wc -l

zcat SV_eQTL/cis.filtered.txt.gz | wc -l
zcat SV_eQTL/trans.filtered.txt.gz | wc -l
zcat SV_apaQTL/cis.filtered.txt.gz | wc -l
zcat SV_apaQTL/trans.filtered.txt.gz | wc -l
zcat SV_sQTL/cis.filtered.txt.gz | wc -l
zcat SV_sQTL/trans.filtered.txt.gz | wc -l
} > /home/wangdy/Projects/Weibin/Downstreams/data_stat/results/stat/QTL_stat.txt
popd
