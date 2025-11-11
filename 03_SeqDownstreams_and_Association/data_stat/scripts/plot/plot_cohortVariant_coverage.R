#! Rscript
# -------------
# FileName     : plot_cohortVariant_coverage
# Author       : EastsunW eastsunw@foxmail.com
# Create at    : 2025-03-06 21:59
# Last Modified: 2025-03-07 08:43
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

variant_coverage <- fread("results/stat/Variant_stat/cohortVariant_coverage.txt") %>%
    filter(Frequancy != "all") %>%
    mutate_at(
        vars("Frequancy"),
        ~ factor(., levels = c("common", "rare"))
    ) %>%
    mutate_at(
        vars("Variant"),
        ~ factor(., levels = rev(c("SNP", "InDel", "MNV", "SV")))
    ) %>%
    mutate(
        Amount = round(coverage / 1e6, 1)
    )

variant_summary <- variant_coverage %>%
    group_by(Variant) %>%
    summarise(max = max(Amount))

pdf(
    "results/plot/cohortVariant_coverage_byFreq.pdf",
    width = 3,
    height = 3,
    onefile = TRUE
)
ggplot(
    data = variant_coverage,
    mapping = aes(
        x = Amount,
        y = Variant,
        fill = Frequancy
    )
) +
    geom_bar(
        stat = "identity",
        position = position_dodge(0.9)
    ) +
    geom_text(
        aes(label = paste0(Amount, "Mb")),
        position = position_dodge(0.9),
        hjust = 0,
        vjust = 0.5
    ) +
    scale_fill_manual(values = c(
        "common" = custom_colors[2],
        "rare" = custom_colors[3]
    )) +
    scale_x_continuous(
        expand = expansion(0, 0),
        limits = c(0, 81),
        n.breaks = 3
    ) +
    labs(
        x = "Concatenated length",
        y = "Variant type",
        fill = "Frequency"
    ) +
    wdy_theme(
        base_size = 12,
        base_line_size = 1,
        legend.position = "none",
        legend.direction = "vertical",
        legend.position.inside = c(0.99, 0.99),
        legend.justification = c(1, 1),
        legend.background = element_blank(),
        axis.line = element_blank(),
        axis.text.x.top = element_blank(),
        axis.ticks.x.top = element_blank(),
        strip.background = element_rect(
            fill = "#ebf1f5",
            linewidth = NA
        )
    )
dev.off()
