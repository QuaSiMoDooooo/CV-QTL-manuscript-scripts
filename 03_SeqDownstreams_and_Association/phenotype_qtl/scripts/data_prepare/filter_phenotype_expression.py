#! python
# -*- coding: UTF-8 -*-
#
# FileName     : filter_expression
# Author       : EastsunW eastsunw@foxmail.com
# Create at    : 2024-12-24 15:42
# Last Modified: 2024-12-24 15:48
# Modified By  : EastsunW
# -------------
# Description  :
# -------------

import pandas as pd

expression_raw = pd.read_csv(
    snakemake.input[1], sep="\t", na_values=["", "NA"])
expression_expressed = expression_raw.set_index("ID")

# Filter by median
if "median" in snakemake.params:
    expression_expressed["median"] = expression_expressed.median(axis=1)
    expression_expressed = expression_expressed[
        expression_expressed["median"] > snakemake.params["median"]
    ]
    expression_expressed = expression_expressed.drop(columns=["median"])

# Filter by mean
if "mean" in snakemake.params:
    expression_expressed["mean"] = expression_expressed.mean(axis=1)
    expression_expressed = expression_expressed[
        expression_expressed["mean"] > snakemake.params["mean"]
    ]
    expression_expressed = expression_expressed.drop(columns=["mean"])

# Filter by missing rate
if "missing_rate" in snakemake.params:
    expression_expressed["missing_rate"] = expression_expressed.isnull().sum(axis=1) / expression_expressed.shape[1]
    expression_expressed = expression_expressed[
        expression_expressed["missing_rate"] < snakemake.params["missing_rate"]
    ]
    expression_expressed = expression_expressed.drop(columns=["missing_rate"])

# Reset index and save to file
expression_expressed = expression_expressed.reset_index()
expression_expressed.to_csv(
    snakemake.output[1], sep="\t", index=False, quoting=False)
