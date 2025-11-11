#! Rscript
# -------------
# FileName     : plot_manhattan_gwas
# Author       : EastsunW eastsunw@foxmail.com
# Create at    : 2025-03-05 16:40
# Last Modified: 2025-03-05 18:37
# Modified By  : EastsunW
# -------------
# Description  : 
# -------------

# 加载必要的包
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(qqman))
setwd("/home/wangdy/Projects/Weibin/Downstreams/data_stat/data/GWAS")

indicator <- "Anti-Ro52_antibody"
indicator_type <- "common"
variant_types <- c("SNP", "InDel", "MNV", "SV")

# 导入所有变异的该指标的GWAS原始结果
gwas_merged <- list()
for (variant_type in variant_types) {
    temp <- fread(
        file.path(
            paste0(variant_type, "_GWAS"),
            paste0("biochem_", indicator_type),
            paste0(indicator, ".assoc.linear")
        )
    ) %>%
        filter(
            TEST == "ADD",
            !is.na(P)
        ) %>%
        mutate(variant_type = variant_type)
    gwas_merged[[variant_type]] <- temp
}

# SNP
png(
    file.path(
        "/home/wangdy/Projects/Weibin/Downstreams/data_stat/results/plot/clinic",
        paste0("SNP_manhattan_", indicator, ".png")
    ),
    width = 6,
    height = 4,
    units = "in",
    res = 600
)
manhattan(
    gwas_merged[["SNP"]],
    col = c("#FA8070", "#83B2D3"),
    cex = 0.6,
    ylim = c(0, 12),
    main = paste0("SNP-GWAS")
)
dev.off()

# InDel
png(
    file.path(
        "/home/wangdy/Projects/Weibin/Downstreams/data_stat/results/plot/clinic",
        paste0("InDel_manhattan_", indicator, ".png")
    ),
    width = 6,
    height = 4,
    units = "in",
    res = 600
)
manhattan(
    gwas_merged[["InDel"]],
    col = c("#FA8070", "#83B2D3"),
    cex = 0.6,
    ylim = c(0, 12),
    main = paste0("InDel-GWAS")
)
dev.off()

# MNV
png(
    file.path(
        "/home/wangdy/Projects/Weibin/Downstreams/data_stat/results/plot/clinic",
        paste0("MNV_manhattan_", indicator, ".png")
    ),
    width = 6,
    height = 4,
    units = "in",
    res = 600
)
manhattan(
    gwas_merged[["MNV"]],
    col = c("#FA8070", "#83B2D3"),
    cex = 0.6,
    ylim = c(0, 12),
    main = paste0("MNV-GWAS")
)
dev.off()

# SV
png(
    file.path(
        "/home/wangdy/Projects/Weibin/Downstreams/data_stat/results/plot/clinic",
        paste0("SV_manhattan_", indicator, ".png")
    ),
    width = 6,
    height = 4,
    units = "in",
    res = 600
)
manhattan(
    gwas_merged[["SV"]],
    col = c("#FA8070", "#83B2D3"),
    cex = 0.6,
    ylim = c(0, 12),
    main = paste0("SV-GWAS")
)
dev.off()
