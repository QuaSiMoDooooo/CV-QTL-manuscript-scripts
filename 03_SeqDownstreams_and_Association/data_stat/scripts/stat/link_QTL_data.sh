#! /bin/bash
# FileName     : link_QTL_data
# Author       : EastsunW eastsunw@foxmail.com
# Create at    : 2024-12-14 16:41
# Last Modified: 2024-12-14 17:29
# Modified By  : EastsunW
# -------------
# Description  : 
# -------------

set -euo pipefail

data_dir=/home/wangdy/Projects/Weibin/Downstreams/phenotype_qtl/results

declare -A qtl_code
qtl_code[eQTL]=expression
qtl_code[apaQTL]=APA
qtl_code[sQTL]=splicing

genotype_filename=variant.genotype.matched.txt
phenotype_filename=phenotype.quantity.matched.txt

qtl_types=(eQTL apaQTL sQTL)
variant_types=(SNP InDel MNV SV)

for qtl_type in ${qtl_types[@]}; do
    for variant_type in ${variant_types[@]}; do
        ln -s ${data_dir}/${variant_type}-${qtl_code[${qtl_type}]}/matched_data/$genotype_filename data/QTL_data/${variant_type}_${qtl_type}.genotype.txt
        ln -s ${data_dir}/${variant_type}-${qtl_code[${qtl_type}]}/matched_data/$phenotype_filename data/QTL_data/${variant_type}_${qtl_type}.phenotype.txt
    done
done
