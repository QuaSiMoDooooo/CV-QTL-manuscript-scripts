#! python
# -*- coding: UTF-8 -*-
#
# FileName     : stat_merge_GWAS
# Author       : EastsunW eastsunw@foxmail.com
# Create at    : 2024-12-09 22:59
# Last Modified: 2024-12-09 23:00
# Modified By  : EastsunW
# -------------
# Description  : 合并每个变异类型的三种GWAS结果
# -------------

import scipy.stats as stats
import concurrent.futures
import pandas as pd
import os
os.chdir("/home/wangdy/Projects/Weibin/Downstreams/data_stat")

variant_types = [
    "SNP",
    "InDel",
    "MNV",
    "SV"
]

biochem_types = [
    "biochem_common",
    "biochem_male",
    "biochem_female"
]


def merge_GWAS_results(result_dir):
    final_df = None
    for biochem_type in biochem_types:
        find_dir = os.path.join(result_dir, biochem_type)
        for root, dir, files in os.walk(find_dir):
            for file in files:
                if not file.endswith(".assoc.linear"):
                    continue
                marker_name = file.split(".")[0]
                temp_file = pd.read_table(
                    os.path.join(root, file),
                    sep=r"\s+",
                )
                temp_file = temp_file[temp_file["TEST"] == "ADD"]
                temp_file = temp_file.dropna(subset=["P"])
                temp_file["MARKER"] = marker_name
                # temp_file["FDR"] = stats.false_discovery_control(
                #     temp_file["P"].values,
                #     method="bh",
                # )
                temp_file = temp_file[[
                    "CHR", "SNP", "BP", "A1", "MARKER", "NMISS", "BETA", "STAT", "P"]]
                if final_df is None:
                    final_df = temp_file
                else:
                    # 按行合并
                    final_df = pd.concat([final_df, temp_file], axis=0)
    return final_df


def process_variant_type(variant_type):
    merged_df = merge_GWAS_results(f"data/{variant_type}_GWAS")
    merged_df.to_csv(
        f"results/stat/{variant_type}_GWAS.merged.txt",
        sep="\t",
        index=False
    )
    merged_df.query("P < 5e-8").to_csv(
        f"results/stat/{variant_type}_GWAS.merged.filtered.txt",
        sep="\t",
        index=False
    )


with concurrent.futures.ProcessPoolExecutor(max_workers=4) as executor:
    executor.map(process_variant_type, variant_types)
