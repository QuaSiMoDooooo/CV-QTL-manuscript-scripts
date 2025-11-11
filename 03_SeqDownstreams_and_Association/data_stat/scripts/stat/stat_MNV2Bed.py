#! python
# -*- coding: UTF-8 -*-
#
# FileName     : stat_MNV
# Author       : EastsunW eastsunw@foxmail.com
# Create at    : 2024-11-26 10:10
# Last Modified: 2024-11-26 10:12
# Modified By  : EastsunW
# -------------
# Description  : 将MNV原始数据转换为bed格式
# -------------

import os
os.chdir("/home/wangdy/Projects/Weibin/Downstreams/data_stat")

MNV_cohort_path = "data/cohort_MNV.txt"

mnv_all_bed    = open("results/bed/MNV.all.bed", "wt")
mnv_common_bed = open("results/bed/MNV.common.bed", "wt")
mnv_rare_bed   = open("results/bed/MNV.rare.bed", "wt")

with open(MNV_cohort_path, "rt") as MNV:
    for line in MNV:
        if line.startswith("#"):
            continue
        line_splited = line.strip().split("\t")
        maf = float(line_splited[9]) if float(line_splited[9]) <= 0.5 else (1 - float(line_splited[9]))
        mnv_all_bed.write("\t".join([
            line_splited[0] if line_splited[0].startswith(
                "chr") else "chr" + line_splited[0],
            str(int(line_splited[1].split(",")[0]) - 1),
            line_splited[1].split(",")[-1],
            line_splited[2],
            str(maf),
            line_splited[6]
        ]) + "\n")
        if maf > 0.05:
            mnv_common_bed.write("\t".join([
                line_splited[0] if line_splited[0].startswith(
                    "chr") else "chr" + line_splited[0],
                str(int(line_splited[1].split(",")[0]) - 1),
                line_splited[1].split(",")[-1],
                line_splited[2],
                str(maf),
                line_splited[6]
            ]) + "\n")
        else:
            mnv_rare_bed.write("\t".join([
                line_splited[0] if line_splited[0].startswith(
                    "chr") else "chr" + line_splited[0],
                str(int(line_splited[1].split(",")[0]) - 1),
                line_splited[1].split(",")[-1],
                line_splited[2],
                str(maf),
                line_splited[6]
            ]) + "\n")
