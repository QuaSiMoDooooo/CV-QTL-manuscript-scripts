#! Rscript
# -------------
# FileName     : plot_MNV_count
# Author       : EastsunW eastsunw@foxmail.com
# Create at    : 2024-11-25 22:59
# Last Modified: 2025-04-29 17:53
# Modified By  : EastsunW
# -------------
# Description  :
# -------------

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(data.table))
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

mnv_data <- fread(
    "results/bed/MNV.all.bed",
    header = FALSE,
    col.names = c("chr", "start", "end", "ID", "MAF", "type")
) %>%
    mutate(freq = ifelse(MAF > 0.05, "common", "rare")) %>%
    group_by(type) %>%
    summarise(n = n()) %>%
    rename(MNVType = type) %>%
    mutate(type = "MNV")
pdf("results/plot/cohort_MNV_freq_count.pdf", width = 4, height = 2, onefile = TRUE)
ggplot(
    data = mnv_data,
    mapping = aes(
        x = MNVType,
        y = n
    )
) +
    geom_bar(
        position = position_dodge(width = 0.5),
        stat = "identity",
        color = "black",
        fill = custom_colors[3]
    ) +
    geom_text(aes(label = n), hjust = 0.5, vjust = -0.75, position = position_dodge(width = 0.5), size = 10 / .pt) +
    wdy_theme(
        axis.line = element_blank(),
        legend.position = "none",
        strip.background = element_rect(fill = "#ebf1f5", linewidth = NA)
    ) +
    scale_x_continuous(breaks = 2:9) +
    scale_y_continuous(
        expand = c(0.05, 0),
        limits = c(0, max(mnv_data$n) * 1.2),
        labels = scales::label_number(scale = 1e-3, suffix = "k")
    ) +
    scale_fill_manual(values = custom_colors) +
    labs(
        x = "MNV Type",
        y = "Number of MNVs"
    )
dev.off()
