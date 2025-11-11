#! Rscript
# -------------
# FileName     : plot_sequencing_stat
# Author       : EastsunW eastsunw@foxmail.com
# Create at    : 2024-11-26 17:05
# Last Modified: 2024-12-16 19:12
# Modified By  : EastsunW
# -------------
# Description  : 画出覆盖度、深度、比对率的信息
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


# 测序深度
pdf(
    "results/plot/sampleSeq_sequencing_depth.pdf",
    width = 3.5,
    height = 2,
    onefile = TRUE
)
fread(
    "data/depth.summary.txt",
    sep = "\t"
) %>%
    rename(Depth = depth) %>%
    arrange(desc(Depth)) %>%
    mutate_at(vars("sample"), ~ factor(., levels = sample)) %>%
    ggplot(
        mapping = aes(
            x = sample,
            y = Depth
        )
    ) +
    geom_bar(
        stat = "identity",
        width = 1,
        fill = custom_colors[1],
        color = NA
    ) +
    geom_hline(
        yintercept = 15,
        linetype = "dashed",
        color = "red"
    ) +
    scale_y_continuous(
        expand = expansion(0, 0)
    ) +
    wdy_theme(
        base_size = 12,
        base_line_size = 1,
        axis.line = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid = element_blank()
    ) +
    labs(
        x = "Sample",
        y = "Depth",
    )
dev.off()

# 比对率
pdf(
    "results/plot/sampleSeq_alignment_rate.pdf",
    width = 3.5,
    height = 2,
    onefile = TRUE
)
fread(
    "data/mappingrate.summary.txt",
    sep = "\t"
) %>%
    rename(Mapping_rate = mapping_rate) %>%
    arrange(desc(Mapping_rate)) %>%
    mutate_at(vars("sample"), ~ factor(., levels = sample)) %>%
    ggplot(
        mapping = aes(
            x = sample,
            y = Mapping_rate
        )
    ) +
    geom_bar(
        stat = "identity",
        width = 1,
        fill = custom_colors[2]
    ) +
    scale_y_continuous(
        expand = expansion(0, 0),
        limits = c(0, 1)
    ) +
    wdy_theme(
        base_size = 12,
        base_line_size = 1,
        axis.line = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid = element_blank()
    ) +
    labs(
        x = "Sample",
        y = "Mapping rate",
    )
dev.off()

# 覆盖度
pdf(
    "results/plot/sampleSeq_alignment_coverage.pdf",
    width = 3.5,
    height = 2,
    onefile = TRUE
)
fread(
    "data/coverage.summary.txt",
    sep = "\t"
) %>%
    rename(Coverage = coverage) %>%
    arrange(desc(Coverage)) %>%
    mutate_at(vars("sample"), ~ factor(., levels = sample)) %>%
    ggplot(
        mapping = aes(
            x = sample,
            y = Coverage
        )
    ) +
    geom_bar(
        stat = "identity",
        width = 1,
        fill = custom_colors[3]
    ) +
    scale_y_continuous(
        expand = expansion(0, 0),
        limits = c(0, 1)
    ) +
    wdy_theme(
        base_size = 12,
        base_line_size = 1,
        axis.line = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid = element_blank()
    ) +
    labs(
        x = "Sample",
        y = "Coverage",
    )
dev.off()
