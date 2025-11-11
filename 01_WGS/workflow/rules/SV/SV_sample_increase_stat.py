#! python
# -*- coding: UTF-8 -*-
#
# FileName     : SV_merge_cohort_step
# Author       : EastsunW eastsunw@foxmail.com
# Create at    : 2024-11-29 09:47
# Last Modified: 2024-11-29 09:48
# Modified By  : EastsunW
# -------------
# Description  : 统计每个样本增加后的共享SV的数量
# -------------

import gzip
import pandas as pd
from concurrent.futures import ThreadPoolExecutor, as_completed
from utils import vcf_variant_info_parser


def process_vcf(file_path):
    if file_path.endswith(".gz"):
        VCF_IN = gzip.open(file_path, 'rt')
    else:
        VCF_IN = open(file_path, 'r')
    SV_counts = {
        "INS": 0,
        "DEL": 0,
        "DUP": 0,
        "INV": 0,
        "BND": 0,
    }
    for line in VCF_IN:
        if line.startswith('#'):
            continue
        line_splited = line.strip().split('\t')
        vcf_infos = vcf_variant_info_parser(line_splited[7])
        if vcf_infos["SVTYPE"] in SV_counts:
            SV_counts[vcf_infos["SVTYPE"]] += 1
    return SV_counts

if __name__ == "__main__":
    result_dict = {
        str(i): {
            "INS": 0,
            "DEL": 0,
            "INV": 0,
            "DUP": 0,
            "BND": 0
        }
        for i in range(1, len(snakemake.input)+1)
    }
    with ThreadPoolExecutor(max_workers=snakemake.threads) as executor:
        future_to_sample = {
            executor.submit(process_vcf, snakemake.input[i]): i+1
            for i in range(len(snakemake.input))
        }
        for future in as_completed(future_to_sample):
            n_sample = future_to_sample[future]
            current_SVs = future.result()
            if isinstance(current_SVs, dict):
                result_dict[str(n_sample)].update(current_SVs)

    df = pd.DataFrame.from_dict(result_dict, orient="index")
    df.index.name = "n_sample"
    df.to_csv(snakemake.output[0], sep="\t")
