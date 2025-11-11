#! python
# -*- coding: UTF-8 -*-
# 
# FileName     : re_classify_chromhmm
# Author       : EastsunW eastsunw@foxmail.com
# Create at    : 2025-01-01 23:01
# Last Modified: 2025-01-01 23:07
# Modified By  : EastsunW
# -------------
# Description  : 
# -------------

import os

input_dir = 'data/chromHMM/'
output_suffix = '_15_coreMarks_hg38lift_converted.bed'

reform_dict = {
    "1_TssA"     : "Promoter",
    "2_TssAFlnk" : "Promoter",
    "3_TxFlnk"   : "Transcriped",
    "4_Tx"       : "Transcriped",
    "5_TxWk"     : "Transcriped",
    "6_EnhG"     : "Enhancer",
    "7_Enh"      : "Enhancer",
    "8_ZNF/Rpts" : "ZNF_repeat",
    "9_Het"      : "Heterchromatin",
    "10_TssBiv"  : "Promoter",
    "11_BivFlnk" : "Promoter",
    "12_EnhBiv"  : "Enhancer",
    "13_ReprPC"  : "Polycomb",
    "14_ReprPCWk": "Polycomb",
    "15_Quies"   : "Quiescent"
}

for root, dirs, files in os.walk(input_dir):
    for file_name in files:
        if file_name.endswith('_15_coreMarks_hg38lift_mnemonics.bed'):
            file_path = os.path.join(root, file_name)
            cell_line = file_name.split('_')[0]
            output_file = os.path.join(root, f'{cell_line}{output_suffix}')
            with open(file_path, 'r') as infile, open(output_file, 'w') as outfile:
                for line in infile:
                    fields = line.strip().split('\t')
                    if len(fields) > 3 and fields[3] in reform_dict:
                        fields[3] = reform_dict[fields[3]]
                    outfile.write('\t'.join(fields) + '\n')
