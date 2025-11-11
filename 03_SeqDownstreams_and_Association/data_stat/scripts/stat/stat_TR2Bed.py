#! python
# -*- coding: UTF-8 -*-
# 
# FileName     : stat_TR2Bed
# Author       : EastsunW eastsunw@foxmail.com
# Create at    : 2024-12-09 09:37
# Last Modified: 2024-12-09 09:59
# Modified By  : EastsunW
# -------------
# Description  : 
# -------------

import os
import gzip
from vcf_helper import vcf_variant_info_parser
os.chdir("/home/wangdy/Projects/Weibin/Downstreams/data_stat")

tr_cohort_path = "data/cohort_TR.vcf.gz"
tr_bed = open("results/bed/TR.normal.bed", "wt")

with gzip.open(tr_cohort_path, "rt") as VCF:
    for line in VCF:
        if line.startswith("#"):
            continue
        line_splited = line.strip().split("\t")
        parsed_info = vcf_variant_info_parser(line_splited[7])
        tr_bed.write("\t".join([
            line_splited[0] if line_splited[0].startswith(
                "chr") else "chr" + line_splited[0],
            str(int(line_splited[1]) - 1),
            parsed_info["END"],
            parsed_info["TRID"]
        ]) + "\n")
