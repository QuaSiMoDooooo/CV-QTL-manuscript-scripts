#! python
# -*- coding: UTF-8 -*-
# 
# FileName     : split_Han945
# Author       : EastsunW eastsunw@foxmail.com
# Create at    : 2025-06-18 22:35
# Last Modified: 2025-06-18 22:38
# Modified By  : EastsunW
# -------------
# Description  : 
# -------------

import os
import gzip
from vcf_helper import vcf_variant_info_parser
os.chdir("/home/wangdy/Projects/Weibin/Downstreams/variant_annotation")
if not os.path.exists("results/SV_compare/Han945/converted"):
    os.makedirs("results/SV_compare/Han945/converted")

def process_Han945(input_file, output_dir):
    with gzip.open(input_file, 'rt') as input:
        output_ins = open(os.path.join(output_dir, "Han945.ins.bed"), 'wt')
        output_del = open(os.path.join(output_dir, "Han945.del.bed"), 'wt')
        output_dup = open(os.path.join(output_dir, "Han945.dup.bed"), 'wt')
        output_inv = open(os.path.join(output_dir, "Han945.inv.bed"), 'wt')
        for line in input:
            if line.startswith("#"):
                continue
            else:
                if not line.startswith("chr"):
                    line = "chr" + line
                line_splited = line.strip().split("\t")
                vcf_infos = vcf_variant_info_parser(line_splited[7])
                if "SVTYPE" in vcf_infos and vcf_infos["SVTYPE"] in ["DEL", "DUP", "INS", "INV"]:
                    if vcf_infos["SVTYPE"] == "DEL":
                        output_del.write("\t".join([
                            line_splited[0],
                            str(int(line_splited[1]) - 1),
                            vcf_infos["END"],
                            line_splited[2],
                        ]) + "\n")
                    elif vcf_infos["SVTYPE"] == "DUP":
                        output_dup.write("\t".join([
                            line_splited[0],
                            str(int(line_splited[1]) - 1),
                            vcf_infos["END"],
                            line_splited[2],
                        ]) + "\n")
                    elif vcf_infos["SVTYPE"] == "INS":
                        output_ins.write("\t".join([
                            line_splited[0],
                            str(int(line_splited[1]) - 1),
                            vcf_infos["END"],
                            line_splited[2],
                        ]) + "\n")
                    elif vcf_infos["SVTYPE"] == "INV":
                        output_inv.write("\t".join([
                            line_splited[0],
                            str(int(line_splited[1]) - 1),
                            vcf_infos["END"],
                            line_splited[2],
                        ]) + "\n")
                else:
                    continue


if __name__ == "__main__":
    process_Han945(
        "data/SV_compare/Han945/Han945.vcf.gz",
        "results/SV_compare/Han945/converted"
    )
