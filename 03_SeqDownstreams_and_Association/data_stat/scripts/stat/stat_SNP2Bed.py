#! python
# -*- coding: UTF-8 -*-
#
# FileName     : stat_SNP_indel_count
# Author       : EastsunW eastsunw@foxmail.com
# Create at    : 2024-11-22 10:46
# Last Modified: 2024-11-22 10:48
# Modified By  : EastsunW
# -------------
# Description  : 统计SNP队列结果中的SNP和InDel的数量，并且根据MAF进行common和rare的划分，SNP包括SNP和MNP，InDel包括insertion和deletion，忽略复等位基因变异和大于50bp的indel
# -------------

import gzip
from vcf_helper import vcf_variant_info_parser
import os
os.chdir("/home/wangdy/Projects/Weibin/Downstreams/data_stat")

SNP_cohort_path = "data/cohort_SNP.vcf.gz"

snp_all_bed      = open("results/bed/SNP.all.bed", "wt")
snp_common_bed   = open("results/bed/SNP.common.bed", "wt")
snp_rare_bed     = open("results/bed/SNP.rare.bed", "wt")

with gzip.open(SNP_cohort_path, "rt") as VCF:
    for line in VCF:
        if line.startswith("#"):
            continue
        line_splited = line.strip().split("\t")
        parsed_info = vcf_variant_info_parser(line_splited[7])
        # 忽略复等位基因变异
        if ";" in line_splited[2]:
            continue
        # SNP
        if len(line_splited[3]) == len(line_splited[4]):
            snp_all_bed.write("\t".join([
                line_splited[0],
                str(int(line_splited[1]) - 1),
                line_splited[1],
                line_splited[2],
                parsed_info["MAF"]
            ]) + "\n")
            if float(parsed_info["MAF"]) > 0.05:
                snp_common_bed.write("\t".join([
                    line_splited[0],
                    str(int(line_splited[1]) - 1),
                    line_splited[1],
                    line_splited[2],
                    parsed_info["MAF"]
                ]) + "\n")
            else:
                snp_rare_bed.write("\t".join([
                    line_splited[0],
                    str(int(line_splited[1]) - 1),
                    line_splited[1],
                    line_splited[2],
                    parsed_info["MAF"]
                ]) + "\n")
