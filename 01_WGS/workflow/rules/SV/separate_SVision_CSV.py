# -*- coding: UTF-8 -*-
#
# FileName     : separate_SVision_CSV
# Author       : EastsunW eastsunw@foxmail.com
# Create at    : 2024-10-09 09:56
# Last Modified: 2024-11-19 21:31
# Modified By  : EastsunW
# -------------
# Description  : 将SVision结果中的复杂变异分离出来
# -------------


simple_sv = open(snakemake.output['simpleSV'], 'wt')
complex_sv = open(snakemake.output['complexSV'], 'wt')

with open(snakemake.input[0], 'rt') as INPUT:
    for line in INPUT:
        if line.startswith('#'):
            simple_sv.write(line)
            complex_sv.write(line)
        else:
            line_splited = line.strip().split('\t')
            if line_splited[4] == "CSV":
                complex_sv.write(line)
            else:
                simple_sv.write(line)
