#! /usr/bin/Rscript
# -------------
# FileName     : encode_phenotype
# Author       : EastsunW eastsunw@foxmail.com
# Create at    : 2024-09-19 21:02
# Last Modified: 2024-11-02 19:30
# Modified By  : EastsunW
# -------------
# Description  : 将基因表达等表型进行编码，替换样本编号为新编码
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

rename_sample <- function(names) {
    unlist(sapply(
        names,
        function(oldname) {
            oldsample <- sub("(wholeblood_\\d{1,3}[AB])(.*)$", "\\1", oldname)
            suffix <- sub("(wholeblood_\\d{1,3}[AB])(.*)$", "\\2", oldname)
            if (endsWith(oldsample, "A")) {
                return(paste0(
                    sample_info[["ID"]][which(sample_info[["blood_1"]] == oldsample)],
                    suffix
                ))
            } else if (endsWith(oldsample, "B")) {
                return(paste0(
                    sample_info[["ID"]][which(sample_info[["blood_2"]] == oldsample)],
                    suffix
                ))
            } else {
                return(paste0("unknown", suffix))
            }
        },
        simplify = TRUE,
        USE.NAMES = FALSE
    ))
}

# 编码基因表达
if (snakemake@wildcards["phenotype"] == "expression") {
    fread(
        snakemake@input[["phenotype"]],
        sep = "\t",
        na.strings = c("NA", "")
    ) %>%
        rename_at(vars(starts_with("wholeblood")), rename_sample) %>%
        fwrite(
            file = snakemake@output[[1]],
            sep = "\t",
            na = "NA",
            quote = FALSE
        )
}

# 编码APA
if (snakemake@wildcards["phenotype"] == "APA") {
    fread(
        snakemake@input[["phenotype"]],
        sep = "\t",
        na.strings = c("NA", "")
    ) %>%
        rename_at(vars(starts_with("wholeblood")), rename_sample) %>%
        fwrite(
            file = snakemake@output[[1]],
            sep = "\t",
            na = "NA",
            quote = FALSE
        )
}

# 编码AS
if (snakemake@wildcards["phenotype"] == "splicing") {
    fread(
        snakemake@input[["phenotype"]],
        sep = "\t",
        na.strings = c("NA", "")
    ) %>%
        rename_at(vars(starts_with("wholeblood")), rename_sample) %>%
        fwrite(
            file = snakemake@output[[1]],
            sep = "\t",
            na = "NA",
            quote = FALSE
        )
}
