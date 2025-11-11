# -*- coding: UTF-8 -*-
#
# FileName     : SV_merge_cohort_step
# Author       : EastsunW eastsunw@foxmail.com
# Create at    : 2024-11-29 09:47
# Last Modified: 2024-11-29 09:48
# Modified By  : EastsunW
# -------------
# Description  : 逐个增加样本，输出共享的SV和所有的SV的数量
# -------------

import gzip
import pandas as pd
from concurrent.futures import ThreadPoolExecutor, as_completed
from typing import Iterable
from utils import vcf_variant_info_parser


def process_vcf(file_path):
    if file_path.endswith(".gz"):
        VCF_IN = gzip.open(file_path, 'rt')
    else:
        VCF_IN = open(file_path, 'r')
    current_SVs = {
        "INS": [],
        "DEL": [],
        "DUP": [],
        "INV": [],
        "BND": [],
    }
    for line in VCF_IN:
        if line.startswith('#'):
            continue
        line_splited = line.strip().split('\t')
        vcf_infos = vcf_variant_info_parser(line_splited[7])
        if vcf_infos["SVTYPE"] in current_SVs:
            current_SVs[vcf_infos["SVTYPE"]].append("_".join([
                line_splited[0],
                line_splited[1],
                vcf_infos["CHR2"],
                vcf_infos["END"]
            ]))
    return current_SVs


def check_overlap(single_sv: str, sv_set: Iterable[str]):
    """
    检查single_sv和sv_set中的SV是否有重叠，对于同一染色体上的SV，重叠定义为重叠部分长度大于等于两个SV各自长度的50%，对于不同染色体上的SV，重叠定义为两个SV的位置差小于等于100bp
    """
    chr1, pos1, chr2, pos2 = single_sv.split("_")
    pos1, pos2 = int(pos1), int(pos2)
    for sv_compare in sv_set:
        chr1_, pos1_, chr2_, pos2_ = sv_compare.split("_")
        pos1_, pos2_ = int(pos1_), int(pos2_)
        if chr1 == chr1_ and chr2 == chr2_:
            if chr1 == chr2:
                # 如果两个SV有重叠，而且重叠部分长度大于等于两个SV各自长度的50%，就返回True
                sv1_len = max(pos1, pos2) - min(pos1, pos2)
                sv2_len = max(pos1_, pos2_) - min(pos1_, pos2_)
                overlap_len = max(pos1, pos2, pos1_, pos2_) - min(pos1, pos2, pos1_, pos2_) -\
                    (max(pos1, pos1_) - min(pos1, pos1_)) -\
                    (max(pos2, pos2_) - min(pos2, pos2_))
                return overlap_len >= 0.5*sv1_len and overlap_len >= 0.5*sv2_len
            else:
                return max(pos1, pos1_) - min(pos1, pos1_) <= 100 and max(pos2, pos2_) - min(pos2, pos2_) <= 100
        else:
            return False


result_dict = {
    str(i): {
        "unique_INS": 0,
        "unique_DEL": 0,
        "unique_INV": 0,
        "unique_DUP": 0,
        "unique_BND": 0,
        "shared_INS": 0,
        "shared_DEL": 0,
        "shared_INV": 0,
        "shared_DUP": 0,
        "shared_BND": 0
    }
    for i in range(1, len(snakemake.input)+1)
}
with ThreadPoolExecutor(max_workers=snakemake.threads) as executor:
    unique_SVs = {
        "INS": set(),
        "DEL": set(),
        "DUP": set(),
        "INV": set(),
        "BND": set(),
    }
    shared_SVs = {
        "INS": set(),
        "DEL": set(),
        "DUP": set(),
        "INV": set(),
        "BND": set(),
    }
    future_to_sample = {
        executor.submit(process_vcf, snakemake.input[i]): i+1
        for i in range(len(snakemake.input))
    }
    for future in as_completed(future_to_sample):
        n_sample = future_to_sample[future]
        current_SVs = future.result()
        if n_sample == 1:
            for svtype in current_SVs:
                unique_SVs[svtype].update(current_SVs[svtype])
                shared_SVs[svtype].update(current_SVs[svtype])
                result_dict[str(n_sample)][f"unique_{svtype}"] = len(
                    unique_SVs[svtype])
                result_dict[str(n_sample)][f"shared_{svtype}"] = len(
                    shared_SVs[svtype])
        else:
            for svtype in current_SVs:
                for per_SV in current_SVs[svtype]:
                    if check_overlap(per_SV, shared_SVs[svtype]):
                        shared_SVs[svtype].add(per_SV)
                    else:
                        unique_SVs[svtype].add(per_SV)
                result_dict[str(n_sample)][f"unique_{svtype}"] = len(
                    unique_SVs[svtype])
                result_dict[str(n_sample)][f"shared_{svtype}"] = len(
                    shared_SVs[svtype])

df = pd.DataFrame.from_dict(result_dict, orient="index")
df.index.name = "n_sample"
df.to_csv(snakemake.output[0], sep="\t")
