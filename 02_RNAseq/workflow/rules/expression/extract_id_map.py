# -*- coding: UTF-8 -*-
# 
# FileName     : extract_id_map
# Author       : EastsunW eastsunw@foxmail.com
# Create at    : 2024-11-14 09:29
# Last Modified: 2024-11-14 09:56
# Modified By  : EastsunW
# -------------
# Description  : 从GTF中提取基因ID和基因名的映射关系，只提取基因
# -------------

import re

output = open(snakemake.output[0], 'w')
output.write('gene_ID\tgene_name\n')
with open(snakemake.input[0]) as input:
    for line in input:
        if line.startswith('#'):
            continue
        line_splited = line.strip().split('\t')
        if line_splited[2] == 'gene':
            gene_id_match = re.match(r'.*gene_id "([^"]+)";', line_splited[8])
            gene_name_match = re.match(r'.*gene_name "([^"]+)";', line_splited[8])
            output.write(f'{gene_id_match.group(1)}\t{gene_name_match.group(1)}\n')
