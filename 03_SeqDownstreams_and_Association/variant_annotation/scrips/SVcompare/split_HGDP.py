#! python
# -*- coding: UTF-8 -*-
#
# FileName     : split_HGDP
# Author       : EastsunW eastsunw@foxmail.com
# Create at    : 2025-02-27 20:56
# Last Modified: 2025-02-27 20:56
# Modified By  : EastsunW
# -------------
# Description  :
# -------------

import os
import gzip
from vcf_helper import vcf_variant_info_parser
os.chdir("/home/wangdy/Projects/Weibin/Downstreams/variant_annotation")


def process_HGDP(input_file, output_dir, svtype):
    with gzip.open(input_file, 'rt') as input:
        output = open(os.path.join(output_dir, f"HGDP.{svtype}.bed"), 'wt')
        for line in input:
            if line.startswith("#"):
                continue
            else:
                if not line.startswith("chr"):
                    line = "chr" + line
                line_splited = line.strip().split("\t")
                svsize = int(line_splited[4].split(":")[1].split("=")[1])
                output.write("\t".join([
                    line_splited[0],
                    str(int(line_splited[1]) - 1),
                    str(int(line_splited[1]) - 1 + svsize),
                    svtype + ":" + line_splited[2],
                ]) + "\n")


if __name__ == "__main__":
    for svtype in ["del", "dup", "ins", "inv"]:
        process_HGDP(
            f"data/SV_compare/HGDP/HGDP.{svtype}.vcf.gz",
            "results/SV_compare/HGDP/converted",
            svtype
        )
