#! Rscript
# -------------
# FileName     : plot_error_stat
# Author       : EastsunW eastsunw@foxmail.com
# Create at    : 2024-12-30 15:49
# Last Modified: 2024-12-30 16:31
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
    "#99999"
)
setwd("/home/wangdy/Projects/Weibin/Downstreams/data_stat")

stat_data <- fread("data/error_stat.txt") %>%
    select(sample, insertion_error, deletion_error, mismatch_error, total_error) %>%
    rename_at(vars(ends_with("_error")), ~ str_remove(., "_error")) %>%
    pivot_longer(
        cols = -sample,
        names_to = "type",
        values_to = "percentage"
    )

error_data <- stat_data %>%
    group_by(type) %>%
    summarise(
        mean = mean(percentage),
        sd = sd(percentage)
    )

pdf(
    "results/plot/sampleSeq_mapping_error.pdf",
    width = 3.1,
    height = 2.5,
    onefile = TRUE
)
ggplot(
    data = error_data,
    mapping = aes(
        x = type,
        y = mean,
        fill = type
    )
) +
    geom_bar(
        width = 0.6,
        stat = "identity",
        color = "black",
    ) +
    geom_errorbar(
        mapping = aes(
            ymin = mean - sd,
            ymax = mean + sd
        ),
        width = 0.2
    ) +
    scale_y_continuous(
        limits = c(0, 1.5),
        expand = c(0, 0)
    ) +
    scale_fill_manual(values = custom_colors) +
    wdy_theme(
        base_size = 10,
        axis.line = element_blank(),
        legend.position = "none"
    ) +
    labs(
        x = "Error Type",
        y = "Percentage (%)"
    )
dev.off()
