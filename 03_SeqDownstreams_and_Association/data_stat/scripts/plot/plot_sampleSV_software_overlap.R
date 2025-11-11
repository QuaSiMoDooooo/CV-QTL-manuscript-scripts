#! Rscript
#
# FileName     : plot_SV_software_overlap
# Author       : EastsunW eastsunw@foxmail.com
# Create at    : 2024-11-22 16:03
# Last Modified: 2024-11-22 16:03
# Modified By  : EastsunW
# -------------
# Description  :
# -------------

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(eulerr))
source("https://gitee.com/eastsunw/personal_code_notebook/raw/master/plot_tools/wdy_theme.r")
custom_colors <- c(
    "#FA7F6F",
    "#8ECFC9",
    "#FFBE7A",
    "#82B0D2",
    "#BEB8DC",
    "#E7DAD2",
    "#999999"
)
setwd("/home/wangdy/Projects/Weibin/Downstreams/data_stat")

overlap_stat <- fread("results/stat/Variant_stat/sampleSV_software_overlap.txt") %>%
    pivot_longer(
        -sample,
        names_to = "combination",
        values_to = "count"
    ) %>%
    group_by(combination) %>%
    summarise(mean = round(mean(count), 0)) %>%
    ungroup() %>%
    mutate_at(vars("combination"), ~gsub("_", "&", .)) %>%
    column_to_rownames("combination") %>%
    t() %>%
    as.data.frame()

softwares <- names(overlap_stat)[!grepl("&", names(overlap_stat))]

pdf("results/plot/sampleSV_software_overlap_euler.pdf", width = 4, height = 4, onefile = TRUE)
plot(
    euler(unlist(overlap_stat), shape = "ellipse"),
    quantities = TRUE,
    fills = list(fill = custom_colors[seq_len(length(softwares))], alpha = 0.8),
    lwd = rep(1.5, length(softwares)),
    labels = list(fontsize = 12)
)
dev.off()

pdf("results/plot/sampleSV_software_overlap_venn.pdf", width = 4, height = 4, onefile = TRUE)
plot(
    venn(unlist(overlap_stat)),
    quantities = TRUE,
    fills = list(fill = custom_colors[seq_len(length(softwares))], alpha = 0.8),
    lwd = rep(1.5, length(softwares)),
    labels = list(fontsize = 12)
)
dev.off()
