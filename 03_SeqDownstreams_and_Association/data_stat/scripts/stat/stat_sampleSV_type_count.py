#! python
# -*- coding: UTF-8 -*-
#
# FileName     : stat_SV_cohort_stat
# Author       : EastsunW eastsunw@foxmail.com
# Create at    : 2024-11-22 10:22
# Last Modified: 2024-11-22 10:37
# Modified By  : EastsunW
# -------------
# Description  : 按类型统计每个样本的非冗余SV数量
# -------------

import os
import gzip
import pandas as pd
os.chdir("/home/wangdy/Projects/Weibin/Downstreams/data_stat")

sample_filtered_sv_path = "data/sample_SV_filtered"

output_dict = {}
for root, dirs, files in os.walk(sample_filtered_sv_path):
    for filename in files:
        if not filename.endswith("by_region.vcf.gz"):
            continue
        file_path = os.path.join(root, filename)
        sample_name = filename.split(".")[0]
        with gzip.open(file_path, "rt") as VCF:
            sv_ins = 0
            sv_del = 0
            sv_inv = 0
            sv_dup = 0
            for line in VCF:
                if line.startswith("#"):
                    continue
                line_splited = line.strip().split("\t")
                sv_type = line_splited[2].split("_")[1]
                if sv_type == "INS":
                    sv_ins += 1
                elif sv_type == "DEL":
                    sv_del += 1
                elif sv_type == "INV":
                    sv_inv += 1
                elif sv_type == "DUP":
                    sv_dup += 1
                else:
                    continue
            output_dict[sample_name] = {
                "INS": sv_ins,
                "DEL": sv_del,
                "INV": sv_inv,
                "DUP": sv_dup
            }

df = pd.DataFrame.from_dict(output_dict, orient="index")
df.index.name = 'sample'
df.to_csv("results/stat/Variant_stat/sampleSV_type_count.txt", sep="\t")
