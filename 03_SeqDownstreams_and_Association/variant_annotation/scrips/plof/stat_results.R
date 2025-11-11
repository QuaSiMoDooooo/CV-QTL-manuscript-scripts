#! Rscript
# -------------
# FileName     : stat_results
# Author       : EastsunW eastsunw@foxmail.com
# Create at    : 2025-01-06 16:28
# Last Modified: 2025-03-03 23:49
# Modified By  : EastsunW
# -------------
# Description  : 将plof的结果 进行统计
# -------------

# 加载必要的包
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(data.table))
setwd("/home/wangdy/Projects/Weibin/Downstreams/variant_annotation")

data_dir <- "data/Variants"
out_dir <- "results/epi_enrichment/plof"
variant_types <- c("SNP", "InDel", "MNV", "SV")
qtl_types <- c("eQTL", "apaQTL", "sQTL")

gene_id <- fread(paste0(out_dir, "/gene_IDs.txt"), sep = "\t")

plof_df <- data.frame()
# 读取所有变异的plof结果
for (variant_type in variant_types) {
    result_path <- file.path(
        out_dir,
        paste0(
            variant_type, ".plof.bed"
        )
    )
    if (file.exists(result_path)) {
        result <- fread(
            result_path,
            sep = "\t",
            header = FALSE
        ) %>%
            select(1:8) %>%
            setnames(c(
                "variant_chr", "variant_start", "variant_end", "variant_ID",
                "feature_chr", "feature_start", "feature_end", "feature_ID"
            )) %>%
            mutate(
                feature_ID = str_extract(feature_ID, "ENST[0-9]{11}\\.[0-9]+")
            ) %>%
            mutate(
                variant_type = variant_type,
                effect_type = "pLoF"
            )
        plof_df <<- rbind(plof_df, result)
    }
}

plof_df_id <- plof_df %>%
    left_join(
        gene_id %>% select(`Transcript stable ID version`, `Gene stable ID version`, `Gene name`),
        by = c("feature_ID" = "Transcript stable ID version")
    ) %>%
    rename(
        "Gene_ID" = "Gene stable ID version",
        "Gene_name" = "Gene name"
    )

fwrite(
    plof_df_id,
    file.path(out_dir, "plof_results.txt"),
    sep = "\t",
    quote = FALSE
)

# 统计全基因覆盖的SV
whole_gene_df <- data.frame()
# 读取所有变异的plof结果
for (variant_type in variant_types) {
    result_path <- file.path(
        out_dir,
        paste0(
            variant_type, ".whole_gene.bed"
        )
    )
    if (file.exists(result_path)) {
        result <- fread(
            result_path,
            sep = "\t",
            header = FALSE
        ) %>%
            select(1:8) %>%
            setnames(c(
                "variant_chr", "variant_start", "variant_end", "variant_ID",
                "feature_chr", "feature_start", "feature_end", "feature_ID"
            )) %>%
            mutate(
                variant_type = variant_type,
                effect_type = "whole_gene_effect"
            )
        whole_gene_df <<- rbind(whole_gene_df, result)
    }
}
whole_gene_df_id <- whole_gene_df %>%
    left_join(
        gene_id %>% select(`Transcript stable ID version`, `Gene stable ID version`, `Gene name`),
        by = c("feature_ID" = "Transcript stable ID version")
    ) %>%
    rename(
        "Gene_ID" = "Gene stable ID version",
        "Gene_name" = "Gene name"
    ) %>%
    filter(!is.na(Gene_name)) %>%
    filter(Gene_name != "") %>%
    mutate(variant_type = substr(variant_ID, 4, 6))

fwrite(
    whole_gene_df_id,
    file.path(out_dir, "whole_gene_results.txt"),
    sep = "\t",
    quote = FALSE
)
