#! python
# -*- coding: UTF-8 -*-
#
# FileName     : process_1KG
# Author       : EastsunW eastsunw@foxmail.com
# Create at    : 2024-12-24 17:27
# Last Modified: 2024-12-24 17:27
# Modified By  : EastsunW
# -------------
# Description  :
# -------------

import os
import gzip
from os import path
from vcf_helper import vcf_variant_info_parser
import multiprocessing as mp
os.chdir("/home/wangdy/Projects/Weibin/Downstreams/variant_annotation")


chroms = ["chr" + str(i) for i in range(1, 23)] + ["chrX"]


def process_fun_1KG(input_file, output_dir, chr):
    if not os.path.exists(output_dir):
        os.makedirs(output_dir, exist_ok=True)
    sv_ins_path = path.join(output_dir, f"1KG_{chr}.ins.bed")
    sv_del_path = path.join(output_dir, f"1KG_{chr}.del.bed")
    sv_dup_path = path.join(output_dir, f"1KG_{chr}.dup.bed")
    sv_inv_path = path.join(output_dir, f"1KG_{chr}.inv.bed")
    ins_output = open(sv_ins_path, "wt")
    del_output = open(sv_del_path, "wt")
    dup_output = open(sv_dup_path, "wt")
    inv_output = open(sv_inv_path, "wt")
    with gzip.open(input_file, "rt") as input:
        for line in input:
            if line.startswith("#"):
                continue
            else:
                if not line.startswith("chr"):
                    line = "chr" + line
                line_splited = line.strip().split("\t")
                vcf_infos = vcf_variant_info_parser(line_splited[7])
                if "SVTYPE" in vcf_infos:
                    if vcf_infos["SVTYPE"] in ["INS", "DEL", "DUP", "INV"]:
                        if vcf_infos["SVTYPE"] == "INS":
                            ins_output.write("\t".join([
                                line_splited[0],
                                str(int(line_splited[1]) - 1),
                                vcf_infos["END"],
                                line_splited[2],
                            ]) + "\n")
                        elif vcf_infos["SVTYPE"] == "DEL":
                            del_output.write("\t".join([
                                line_splited[0],
                                str(int(line_splited[1]) - 1),
                                vcf_infos["END"],
                                line_splited[2],
                            ]) + "\n")
                        elif vcf_infos["SVTYPE"] == "DUP":
                            dup_output.write("\t".join([
                                line_splited[0],
                                str(int(line_splited[1]) - 1),
                                vcf_infos["END"],
                                line_splited[2],
                            ]) + "\n")
                        elif vcf_infos["SVTYPE"] == "INV":
                            inv_output.write("\t".join([
                                line_splited[0],
                                str(int(line_splited[1]) - 1),
                                vcf_infos["END"],
                                line_splited[2],
                            ]) + "\n")


if __name__ == "__main__":
    with mp.Pool(23) as pool:
        pool.starmap(
            process_fun_1KG,
            [(
                f"data/SV_compare/1KG/1KG_{chr}.alltype.vcf.gz",
                "results/SV_compare/1KG/converted",
                chr
            ) for chr in chroms
            ]
        )
