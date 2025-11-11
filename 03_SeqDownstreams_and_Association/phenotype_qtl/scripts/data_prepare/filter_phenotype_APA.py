#! python
# -*- coding: UTF-8 -*-
#
# FileName     : filter_APA
# Author       : EastsunW eastsunw@foxmail.com
# Create at    : 2024-12-24 15:50
# Last Modified: 2024-12-24 15:50
# Modified By  : EastsunW
# -------------
# Description  :
# -------------

import pandas as pd

apa_raw = pd.read_csv(
    snakemake.input[1],
    sep="\t",
    na_values=["", "NA"]
)
apa_filtered = apa_raw[
    (apa_raw["Gene_Name"] != "NA") &
    (apa_raw["Gene_Name"].notna())
]
apa_filtered = apa_filtered.rename(columns={"APA_ID": "ID"})
apa_filtered = apa_filtered.loc[
    :,
    ["ID"] + [col for col in apa_filtered.columns if col.endswith(".PAU")]
]
apa_filtered.columns = ["ID"] + [
    col.replace(".PAU", "")
    for col in apa_filtered.columns if col.endswith(".PAU")
]
apa_filtered = apa_filtered.set_index("ID")
if "missing_rate" in snakemake.params:
    apa_filtered["missing_rate"] = apa_filtered.isna().sum(axis=1) / \
        (apa_filtered.shape[1] - 1)
    apa_filtered = apa_filtered[apa_filtered["missing_rate"]
                                < snakemake.params["missing_rate"]]
    apa_filtered = apa_filtered.drop(columns=["missing_rate"])
apa_filtered = apa_filtered.fillna(0).reset_index()
apa_filtered.to_csv(snakemake.output[1], sep="\t", index=False, quoting=False)
