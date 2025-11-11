#! python
# -*- coding: UTF-8 -*-
#
# FileName     : filter_phenotype_AS
# Author       : EastsunW eastsunw@foxmail.com
# Create at    : 2024-12-24 16:03
# Last Modified: 2024-12-24 16:03
# Modified By  : EastsunW
# -------------
# Description  :
# -------------

import pandas as pd

as_raw = pd.read_csv(snakemake.input[1], sep="\t", na_values=["", "NA"])
as_filtered = as_raw.drop(columns=["#Chr", "start", "end"])
as_filtered = as_filtered.set_index("ID")
if "missing_rate" in snakemake.params:
    as_filtered["missing_rate"] = as_filtered.isna().sum(axis=1) / \
        (as_filtered.shape[1] - 1)
    as_filtered = as_filtered[as_filtered["missing_rate"]
                              < snakemake.params["missing_rate"]]
    as_filtered = as_filtered.drop(columns=["missing_rate"])
as_filtered = as_filtered.fillna(0).reset_index()
as_filtered.to_csv(snakemake.output[1], sep="\t", index=False, quoting=False)
