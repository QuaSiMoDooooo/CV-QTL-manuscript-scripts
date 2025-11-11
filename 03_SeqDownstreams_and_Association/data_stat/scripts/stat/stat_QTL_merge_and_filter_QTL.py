#! python
# -------------
# FileName     : stat_QTL_filterQTL
# Author       : EastsunW eastsunw@foxmail.com
# Create at    : 2024-12-10 23:00
# Last Modified: 2024-12-10 23:02
# Modified By  : EastsunW
# -------------
# Description  : 合并cis和trans结果，并过滤FDR
# -------------

import pandas as pd
import os
os.chdir("/home/wangdy/Projects/Weibin/Downstreams/data_stat")

variant_types = ["SNP", "InDel", "MNV", "SV"]
# variant_types = ["SV"]
phenotypes = ["eQTL", "apaQTL", "sQTL", "meQTL"]
qtl_types = ["cis", "trans"]

fdr_threshold = 0.05

QTL_dir = "data/QTL_results"

stat_dict = {}
for variant_type in variant_types:
    for phenotype in phenotypes:
        if f"{variant_type}_{phenotype}" not in stat_dict:
            stat_dict[f"{variant_type}_{phenotype}"] = {}
        merge_df = None
        for qtl_type in qtl_types:
            df = pd.read_csv(
                f"{QTL_dir}/{variant_type}_{phenotype}/{qtl_type}.txt",
                sep="\t"
            ).query("FDR < @fdr_threshold")
            if merge_df is None:
                merge_df = df
            else:
                merge_df = pd.concat([merge_df, df], ignore_index=True)
            stat_dict[f"{variant_type}_{phenotype}"][qtl_type] = df.shape[0]
        merge_df.to_csv(
            f"results/stat/QTL/{variant_type}_{phenotype}.filtered.txt",
            sep="\t",
            index=False
        )

qtl_stat_df = pd.DataFrame(stat_dict)
qtl_stat_df.index.name = "QTL_type"
# 增加一行总计
qtl_stat_df.loc["Total"] = qtl_stat_df.sum()
qtl_stat_df.to_csv("results/stat/QTL/QTL_stat.txt", sep="\t")
