#! python
# -*- coding: UTF-8 -*-
# 
# FileName     : stat_sampleSV_filter_remain
# Author       : EastsunW eastsunw@foxmail.com
# Create at    : 2024-11-26 16:23
# Last Modified: 2024-12-14 16:07
# Modified By  : EastsunW
# -------------
# Description  : 统计每个样本经过每一步筛选后剩余的SV数量
# -------------

import pandas as pd
import gzip
import os

os.chdir("/home/wangdy/Projects/Weibin/Downstreams/data_stat")

sv_filtered_dir = "data/sample_SV_filtered"

step_remain_dict = {}
for root, dirs, files in os.walk(sv_filtered_dir):
    for filename in files:
        if not filename.endswith(".vcf.gz"):
            continue
        sample_name = filename.split(".")[0]
        filter_step = filename.split(".")[1]
        file_path = os.path.join(root, filename)
        with gzip.open(file_path, 'rt') as f:
            line_count = sum(1 for line in f if not line.startswith('#'))
            if sample_name not in step_remain_dict:
                step_remain_dict[sample_name] = {}
            if filter_step not in step_remain_dict[sample_name]:
                step_remain_dict[sample_name][filter_step] = line_count

step_remain_df = pd.DataFrame.from_dict(
    step_remain_dict,
    orient='index'
)
step_remain_df.index.name = 'sample'
step_remain_df.to_csv(
    "results/stat/Variant_stat/sampleSV_filter_remain.txt",
    sep="\t"
)
