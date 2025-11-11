#! python
# -*- coding: UTF-8 -*-
# 
# FileName     : stat_variant_amount_byFreq
# Author       : EastsunW eastsunw@foxmail.com
# Create at    : 2024-12-12 09:49
# Last Modified: 2024-12-14 16:20
# Modified By  : EastsunW
# -------------
# Description  : 统计不同频率下的变异位点数量
# -------------

import os
os.chdir("/home/wangdy/Projects/Weibin/Downstreams/data_stat")

input_path = {
    "SNP": {
        "common": "results/bed/SNP.common.bed",
        "rare": "results/bed/SNP.rare.bed",
    },
    "InDel": {
        "common": "results/bed/InDel.common.bed",
        "rare": "results/bed/InDel.rare.bed",
    },
    "MNV": {
        "common": "results/bed/MNV.common.bed",
        "rare": "results/bed/MNV.rare.bed",
    },
    "SV": {
        "common": "results/bed/SV.common.bed",
        "rare": "results/bed/SV.rare.bed",
    },
}

output_path = "results/stat/Variant_stat/cohortVariant_amount_byFreq.txt"

with open(output_path, "w") as f:
    f.write("Type\tFreq\tAmount\n")
    for type in input_path:
        for freq in input_path[type]:
            amount = int(os.popen(f"wc -l < {input_path[type][freq]}").read().strip())
            f.write(f"{type}\t{freq}\t{amount}\n")
