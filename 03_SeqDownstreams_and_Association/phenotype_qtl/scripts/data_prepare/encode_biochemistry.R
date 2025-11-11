#! /usr/bin/Rscript
# -------------
# FileName     : encode_phenotypes
# Author       : EastsunW eastsunw@foxmail.com
# Create at    : 2024-09-06 17:16
# Last Modified: 2024-09-28 12:00
# Modified By  : EastsunW
# -------------
# Description  : 将原始的生化指标信息中的姓名替换为编码，并导出
# -------------

# 加载必要的包
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(data.table))

# 导入样本信息
sample_info <- fread(
    snakemake@input[["sample_info"]],
    sep = "\t",
    na.strings = c("NA", "")
)

# 编码生化指标
biochem_raw <- fread(
    snakemake@input[["data"]],
    sep = "\t",
    na.strings = c("NA", "")
)

# 用编码替换真实姓名
raw_names <- colnames(biochem_raw)[8:length(colnames(biochem_raw))]
name_codes <- sample_info$ID[match(raw_names, sample_info$姓名)]
colnames(biochem_raw)[8:length(colnames(biochem_raw))] <- name_codes

fwrite(
    biochem_raw,
    snakemake@output[[1]],
    sep = "\t"
)
