#! python
# -*- coding: UTF-8 -*-
#
# FileName     : stat_SV_cohortSV_count
# Author       : EastsunW eastsunw@foxmail.com
# Create at    : 2024-11-22 10:32
# Last Modified: 2024-11-22 10:32
# Modified By  : EastsunW
# -------------
# Description  :
# -------------

import os
import gzip
from vcf_helper import vcf_variant_info_parser
os.chdir("/home/wangdy/Projects/Weibin/Downstreams/data_stat")

SV_cohort_path = "data/cohort_SV.vcf.gz"

sv_all_bed = open("results/bed/SV.all.bed", "wt")
sv_common_bed = open("results/bed/SV.common.bed", "wt")
sv_rare_bed = open("results/bed/SV.rare.bed", "wt")

with gzip.open(SV_cohort_path, "rt") as VCF:
    for line in VCF:
        if line.startswith("#"):
            continue
        line_splited = line.strip().split("\t")
        parsed_info = vcf_variant_info_parser(line_splited[7])
        if parsed_info["SVTYPE"] not in ["INS", "DEL", "INV", "DUP"]:
            continue
        maf = float(parsed_info["AF"]) if float(
            parsed_info["AF"]) <= 0.5 else 1 - float(parsed_info["AF"])
        sv_all_bed.write("\t".join([
            line_splited[0] if line_splited[0].startswith(
                "chr") else "chr" + line_splited[0],
            str(int(line_splited[1]) - 1),
            parsed_info["END"],
            line_splited[2],
            str(maf),
            parsed_info["SVTYPE"],
            str(abs(int(parsed_info["AVGLEN"]))),
            parsed_info["SUPPORT"]
        ]) + "\n")
        if maf > 0.05:
            sv_common_bed.write("\t".join([
                line_splited[0] if line_splited[0].startswith(
                    "chr") else "chr" + line_splited[0],
                str(int(line_splited[1]) - 1),
                parsed_info["END"],
                line_splited[2],
                str(maf),
                parsed_info["SVTYPE"],
                str(abs(int(parsed_info["AVGLEN"]))),
                parsed_info["SUPPORT"]
            ]) + "\n")
        else:
            sv_rare_bed.write("\t".join([
                line_splited[0] if line_splited[0].startswith(
                    "chr") else "chr" + line_splited[0],
                str(int(line_splited[1]) - 1),
                parsed_info["END"],
                line_splited[2],
                str(maf),
                parsed_info["SVTYPE"],
                str(abs(int(parsed_info["AVGLEN"]))),
                parsed_info["SUPPORT"]
            ]) + "\n")
