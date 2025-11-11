#! python
# -*- coding: UTF-8 -*-
# 
# FileName     : split_Iceland
# Author       : EastsunW eastsunw@foxmail.com
# Create at    : 2024-12-25 10:21
# Last Modified: 2024-12-25 16:51
# Modified By  : EastsunW
# -------------
# Description  : 
# -------------

import os
import gzip
from vcf_helper import vcf_variant_info_parser
os.chdir("/home/wangdy/Projects/Weibin/Downstreams/variant_annotation")

def process_Iceland(input_file, output_dir):
    with gzip.open(input_file, 'rt') as input:
        output_ins = open(os.path.join(output_dir, "Iceland.ins.bed"), 'wt')
        output_del = open(os.path.join(output_dir, "Iceland.del.bed"), 'wt')
        for line in input:
            if line.startswith("#"):
                continue
            else:
                if not line.startswith("chr"):
                    line = "chr" + line
                line_splited = line.strip().split("\t")
                vcf_infos = vcf_variant_info_parser(line_splited[7])
                if "SVTYPE" in vcf_infos and vcf_infos["SVTYPE"] in ["DEL", "DUP", "INS", "INV"]:
                    if vcf_infos["SVTYPE"] == "INS":
                        output_ins.write("\t".join([
                            line_splited[0],
                            str(int(line_splited[1]) - 1),
                            vcf_infos["END"],
                            line_splited[2],
                        ]) + "\n")
                    elif vcf_infos["SVTYPE"] == "DEL":
                        output_del.write("\t".join([
                            line_splited[0],
                            str(int(line_splited[1]) - 1),
                            vcf_infos["END"],
                            line_splited[2],
                        ]) + "\n")
                else:
                    continue

if __name__ == "__main__":
    process_Iceland(
        "data/SV_compare/Iceland/Iceland.vcf.gz",
        "results/SV_compare/Iceland/converted"
    )
