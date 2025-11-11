# -------------
# FileName     : normalize_biochemistry
# Author       : EastsunW eastsunw@foxmail.com
# Create at    : 2024-09-20 19:45
# Last Modified: 2024-11-30 13:10
# Modified By  : EastsunW
# -------------
# Description  : 使用z-score标准化对数值型的生化指标进行标准化，每行一个指标，每列一个样本
# -------------

# 加载必要的包
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(data.table))


z_norm <- function(x) {
    x <- as.numeric(x)
    (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
}

data_raw <- fread(snakemake@input[[1]])
data_normalized <- data_raw %>%
    column_to_rownames(colnames(.)[1]) %>%
    t() %>%
    as.data.frame() %>%
    mutate(across(where(is.numeric), z_norm)) %>%
    mutate_all(~replace_na(., 0)) %>%
    rownames_to_column("IID") %>%
    mutate(FID = IID) %>%
    relocate(FID, .before = "IID")
fwrite(
    data_normalized,
    snakemake@output[["quantity"]],
    sep = "\t",
    quote = FALSE,
    col.names = FALSE
)
fwrite(
    data.frame(name=colnames(data_normalized)[-c(1, 2)]) %>%
        rownames_to_column("index") %>%
        mutate_at(vars("name"), ~gsub(" ", "_", .)) %>%
        mutate_at(vars("name"), ~gsub("/", "_", .)),
    snakemake@output[["marker"]],
    sep = "\t",
    quote = FALSE,
    col.names = FALSE
)
