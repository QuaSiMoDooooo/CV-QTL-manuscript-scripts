#! Rscript
# -------------
# FileName     : plot_variant_amount_byFreq
# Author       : EastsunW eastsunw@foxmail.com
# Create at    : 2024-12-12 09:46
# Last Modified: 2024-12-14 16:23
# Modified By  : EastsunW
# -------------
# Description  : 画一个条形图，展示不同频率的变异数量
# -------------

# 加载必要的包
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


variant_counts <- fread("results/stat/Variant_stat/cohortVariant_amount_byFreq.txt") %>%
    mutate_at(
        vars("Freq"),
        ~ factor(., levels = c("common", "rare"))
    ) %>%
    mutate_at(
        vars("Type"),
        ~ factor(., levels = rev(c("SNP", "InDel", "MNV", "SV")))
    ) %>%
    mutate(
        Amount = Amount / 1000
    )

variant_summary <- variant_counts %>%
    group_by(Type) %>%
    summarise(
        max = max(Amount),
        min = min(Amount),
    ) %>%
    column_to_rownames("Type")

temp_plot <- ggplot(
    data = variant_counts,
    mapping = aes(
        x = Amount,
        y = Type,
        fill = Freq
    )
) +
    geom_bar(
        stat = "identity",
        position = position_dodge(0.9)
    ) +
    geom_text(
        aes(label = paste0(round(Amount, 0), "k")),
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
        limits = c(0, max(variant_summary$max * 1.2)),
        n.breaks = 3,
        breaks = c(0, 100, 200)
    ) +
    labs(
        x = "Number of variant (k)",
        y = "Variant type",
        fill = "Frequency"
    ) +
    wdy_theme(
        base_size = 12,
        base_line_size = 1,
        legend.position = "inside",
        legend.direction = "vertical",
        legend.position.inside = c(0.99, 0.01),
        legend.justification = c(1, 0),
        legend.background = element_blank(),
        axis.line = element_blank(),
        axis.text.x.top = element_blank(),
        axis.ticks.x.top = element_blank(),
        strip.background = element_rect(
            fill = "#ebf1f5",
            linewidth = NA
        )
    )
legend <- cowplot::get_plot_component(temp_plot, "guide-box-inside", return_all = TRUE)
breaked_plot <- temp_plot +
    scale_x_break(
        breaks = c(300, 480),
        scales = 1,
        expand = c(0, 0),
        ticklabels = c(600, 900, 1200)
    ) +
    scale_x_break(
        breaks = c(1550, 5800),
        scales = 0.9,
        expand = c(0, 0),
        ticklabels = c(6000, 9000, 12000)
    ) +
    theme(legend.position = "none")
final_plot <- ggplotify::as.ggplot(print(breaked_plot)) +
    ggimage::geom_subview(
        x = 0.9, y = 0.4,
        subview = legend
    )
pdf(
    "results/plot/cohortVariant_amount_byFreq.pdf",
    width = 7,
    height = 3,
    onefile = TRUE
)
final_plot
dev.off()
