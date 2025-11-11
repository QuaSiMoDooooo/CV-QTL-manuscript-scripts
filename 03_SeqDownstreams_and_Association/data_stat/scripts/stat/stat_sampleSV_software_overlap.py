#! python
# -*- coding: UTF-8 -*-
#
# FileName     : stat_SV_filter_remain
# Author       : EastsunW eastsunw@foxmail.com
# Create at    : 2024-11-22 09:56
# Last Modified: 2024-12-11 13:21
# Modified By  : EastsunW
# -------------
# Description  : 统计SV的样本中不同软件SV的重叠数量
# -------------

import os
import gzip
import pandas as pd
from itertools import combinations
os.chdir("/home/wangdy/Projects/Weibin/Downstreams/data_stat")


def generate_combo(names:list, prefix=None):
    n = len(names)
    keys = []
    for r in range(1, n + 1):
        for combo in combinations(range(n), r):
            if names:
                combo_str = '_'.join([names[i] for i in combo])
            else:
                combo_str = '_'.join([str(i) for i in combo])
            if prefix:
                keys.append(f"{prefix}_{combo_str}")
            else:
                keys.append(combo_str)
    return keys


sample_merge_sv_path = "data/sample_SV_merged"


output_dict = {}
for root, dirs, files in os.walk(sample_merge_sv_path):
    for filename in files:
        file_path = os.path.join(root, filename)
        sample_name = filename.split(".")[0]
        with gzip.open(file_path, "rt") as VCF:
            for line in VCF:
                if line.startswith("##"):
                    continue
                elif line.startswith("#"):
                    line_splited = line.strip().split("\t")
                    softwares = line_splited[-(len(line_splited)-9):]
                    softwares_combinations_dict = {
                        key: 0
                        for key in generate_combo(names=softwares)
                    }
                    continue
                line_splited = line.strip().split("\t")
                parse_result = "".join([
                    "0" if i == "./.:0:0:0:." else "1"
                    for i in line_splited[-(len(softwares)):]
                ])
                count_ones = parse_result.count("1")
                for combo in combinations(range((len(softwares))), count_ones):
                    if all(parse_result[i] == "1" for i in combo):
                        key = '_'.join([softwares[i] for i in combo])
                        softwares_combinations_dict[key] += 1
                        break
        output_dict[sample_name] = softwares_combinations_dict

df = pd.DataFrame.from_dict(output_dict, orient="index")
df.index.name = 'sample'
df.to_csv("results/stat/Variant_stat/sampleSV_software_overlap.txt", sep="\t")
