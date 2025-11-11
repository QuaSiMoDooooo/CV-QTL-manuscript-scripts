#! python
# -*- coding: UTF-8 -*-
#
# FileName     : plink_bim_postprocess
# Author       : EastsunW eastsunw@foxmail.com
# Create at    : 2025-04-01 09:57
# Last Modified: 2025-04-01 09:57
# Modified By  : EastsunW
# -------------
# Description  : 将bim文件中的ref和alt进行修复
# -------------

var_dict = {}
with open(snakemake.input["position"], "rt") as INFO:
    for line in INFO:
        line = line.strip().split("\t")
        if line[0] == "ID":
            continue
        else:
            var_dict[line[0]] = {
                "chrom": line[1].replace("chr", ""),
                "pos": line[2],
                "ref": line[4].replace(",", "_") if line[4] != "." else "REF",
                "alt": line[5].replace(",", "_") if line[5] != "." else "ALT",
            }

with open(snakemake.output[0], "wt") as OUTPUT:
    with open(snakemake.input["bim"], "rt") as INPUT:
        for line in INPUT:
            line = line.strip().split("\t")
            if line[1] in var_dict:
                OUTPUT.write("\t".join([
                    var_dict[line[1]]["chrom"],
                    line[1],
                    "0",
                    var_dict[line[1]]["pos"],
                    var_dict[line[1]]["ref"],
                    var_dict[line[1]]["alt"],
                ]) + "\n")
