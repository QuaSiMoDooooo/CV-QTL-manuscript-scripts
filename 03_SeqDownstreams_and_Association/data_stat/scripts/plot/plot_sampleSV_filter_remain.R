#! Rscript
# -------------
# FileName     : SV_filter_remain
# Author       : EastsunW eastsunw@foxmail.com
# Create at    : 2024-11-22 09:46
# Last Modified: 2025-06-03 21:39
# Modified By  : EastsunW
# -------------
# Description  : 每一步过滤后剩余的SV数量统计
# -------------

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(ggbreak))
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

filter_data <- fread("results/stat/Variant_stat/sampleSV_filter_remain.txt") %>%
    pivot_longer(
        cols = -sample,
        names_to = "filter",
        values_to = "remained"
    ) %>%
    mutate_at(
        vars("filter"),
        ~ gsub("by_", "", .)
    ) %>%
    mutate_at(
        vars("filter"),
        ~ factor(., levels = rev(c("software", "depth", "length", "region")))
    )
pdf("results/plot/sapmleSV_filter_remain.pdf", width = 3.5, height = 2.5, onefile = FALSE)
ggplot(
    data = filter_data,
    mapping = aes(
        x = rev(remained),
        y = rev(filter),
        group = filter,
        fill = filter
    )
) +
    stat_summary(
        fun = mean,
        geom = "bar",
        width = 0.8,
        color = "black"
    ) +
    stat_summary(
        fun.data = mean_se,
        geom = "errorbar",
        width = 0.2,
        linewidth = 0.5
    ) +
    geom_text(
        stat = "summary",
        fun = mean,
        hjust = 1.5,
        aes(label = round(after_stat(x), 0)),
        size = 12 / .pt,
        color = "black"
    ) +
    scale_x_continuous(expand = c(0, 0), breaks = c(0, 100)) +
    scale_x_break(
        c(100, 20000), scales = 5,
        space = 0.6,
        ticklabels = seq(20500, 21500, by = 1000),
    ) +
    wdy_theme(
        base_size = 12,
        base_line_size = 1,
        legend.position = "none",
        axis.line.x = element_blank(),
        axis.line.y = element_blank(),
        axis.line.x.top = element_blank(),
        axis.text.x.top = element_blank(),
        axis.ticks.x.top = element_blank(),
        axis.line.y.right = element_blank(),
        axis.ticks.y.right = element_blank(),
        axis.text.y.right = element_blank()
    ) +
    scale_fill_manual(values = custom_colors) +
    labs(
        x = "Remaining SVs",
        y = "Filter steps"
    )
dev.off()
