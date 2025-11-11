#! python
# -*- coding: UTF-8 -*-
# 
# FileName     : process_HGSVC2
# Author       : EastsunW eastsunw@foxmail.com
# Create at    : 2024-12-24 23:31
# Last Modified: 2024-12-25 17:23
# Modified By  : EastsunW
# -------------
# Description  : 
# -------------

import os
import gzip
from vcf_helper import vcf_variant_info_parser
os.chdir("/home/wangdy/Projects/Weibin/Downstreams/variant_annotation")

def process_HGSVC(input_file, output_dir):
    if not os.path.exists(output_dir):
        os.makedirs(output_dir, exist_ok=True)
    with gzip.open(input_file, 'rt') as input:
        output_ins = open(os.path.join(output_dir, "HGSVC.ins.bed"), 'wt')
        output_del = open(os.path.join(output_dir, "HGSVC.del.bed"), 'wt')
        output_dup = open(os.path.join(output_dir, "HGSVC.dup.bed"), 'wt')
        output_inv = open(os.path.join(output_dir, "HGSVC.inv.bed"), 'wt')
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
                        if not "END" in vcf_infos:
                            if "SVLEN" in vcf_infos and vcf_infos["SVLEN"] != ".":
                                vcf_infos["END"] = str(int(line_splited[1]) - 1 + int(vcf_infos["SVLEN"]))
                            else:
                                continue
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
                    elif vcf_infos["SVTYPE"] == "INV":
                        output_inv.write("\t".join([
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
                else:
                    continue

if __name__ == "__main__":
    process_HGSVC(
        "data/SV_compare/HGSVC/HGSVC.sv.vcf.gz",
        "results/SV_compare/HGSVC/converted"
    )
