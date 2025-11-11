# -*- coding: UTF-8 -*-
#
# FileName     : SV_add_AF_cohort
# Author       : EastsunW eastsunw@foxmail.com
# Create at    : 2024-07-03 00:04
# Last Modified: 2024-07-03 00:04
# Modified By  : EastsunW
# -------------
# Description  : 给合并后的人群SV文件加上AF信息
# -------------

import gzip
from utils import modify_vcf_info, vcf_sample_info_parser, vcf_variant_info_parser

def find_nearst_subsv(sv_format: str, sv_start: int, sv_end: int, sv_genotype_str: str) -> str:
    if not isinstance(sv_start, int):
        sv_start = int(sv_start)
    if not isinstance(sv_end, int):
        sv_end = int(sv_end)
    subsv_genos = sv_genotype_str.split(";")
    distances = {}
    if len(subsv_genos) == 1:
        return subsv_genos[0]
    else:
        for subsv_geno in subsv_genos:
            subsv_genotype = vcf_sample_info_parser(sv_format, subsv_geno)
            if subsv_genotype["GT"] == "./.":
                continue
            sv_ranges = subsv_genotype["RG"].split("-")
            start2 = int(sv_ranges[0].split("_")[-1])
            end2 = int(sv_ranges[1].split("_")[-1])
            distances[subsv_geno] = abs(sv_start - start2) + abs(sv_end - end2)
        if distances == {}:
            return "./.:0:0:0:."
    sorted_distances = sorted(distances.items(), key=lambda x: x[1])
    return sorted_distances[0][0]

with gzip.open(snakemake.output[0], "w") as OUTPUT:
    with gzip.open(snakemake.input[0], "rt") as VCF:
        for line in VCF:
            if line.startswith("#"):
                OUTPUT.write(line.encode('utf-8'))
            else:
                line_splited = line.strip().split("\t")
                sv_infos = vcf_variant_info_parser(line_splited[7])
                n_sample = len(line_splited[9:])
                n_allele = 0
                for sample_genotype in line_splited[9:]:
                    try:
                        nearst_sample_sv_gt = find_nearst_subsv(
                            sv_format=line_splited[8],
                            sv_start=line_splited[1],
                            sv_end=sv_infos["END"],
                            sv_genotype_str=sample_genotype
                        )
                        sample_info = vcf_sample_info_parser(
                            line_splited[8], nearst_sample_sv_gt)
                        match sample_info["GT"]:
                            case "0/1":
                                n_allele += 1
                            case "1/1":
                                n_allele += 2
                            case "0/0":
                                continue
                            case "./.":
                                continue
                    except Exception as e:
                        print(line, e)
                        print(sample_genotype)
                        exit(1)

                ac = n_allele
                af = round(ac / n_sample / 2, 3)
                new_info = modify_vcf_info(
                    vcf_info=line_splited[7],
                    key=["AF", "AC"],
                    value=[str(af), str(ac)])
                OUTPUT.write(("\t".join(
                    line_splited[:7] + [new_info] + line_splited[8:]) + "\n").encode('utf-8'))
