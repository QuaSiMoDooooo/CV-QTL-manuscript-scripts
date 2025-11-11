# -*- coding: UTF-8 -*-
# 
# FileName     : generate_MNV_ID
# Author       : EastsunW eastsunw@foxmail.com
# Create at    : 2024-09-24 11:35
# Last Modified: 2024-11-02 10:56
# Modified By  : EastsunW
# -------------
# Description  : 给MNV的鉴定结果生成唯一的变异ID
# -------------

import pandas as pd


def extract_position_start(pos_str: str):
    positions = [int(pos) for pos in pos_str.split(",")]
    return min(positions)


def set_ID(row):
    # ID = “HN”+两位染色体+MNVType+从1开始的顺序，位数不够补0
    return "HNMNV" + str(row["#chr"]).zfill(2) + str(row["MNVType"]).zfill(2) + str(row["index"]).zfill(5)


df = pd.read_csv(
    snakemake.input[0], 
    sep="\t",
    low_memory=False
)
# 提出开始位置
df['position_start'] = df['pos'].apply(extract_position_start)
# 排序
df.sort_values(
    by=["#chr", "MNVType", "position_start", "adjust_AC"],
    ascending=[True, True, True, False],
    inplace=True
)
df.drop(columns=["position_start"], inplace=True)
# 分组指定编号
df['index'] = df.groupby(["#chr", "MNVType"]).cumcount() + 1
df['MNVID'] = df.apply(set_ID, axis=1)
df.drop(columns=["index", "rsID"], inplace=True)
df.to_csv(snakemake.output[0], sep="\t", index=False)
