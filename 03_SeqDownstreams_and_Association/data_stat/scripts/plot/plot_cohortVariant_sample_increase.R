#! Rscript
# -------------
# FileName     : plot_cohortSNP_sample_increase
# Author       : EastsunW eastsunw@foxmail.com
# Create at    : 2024-12-10 15:54
# Last Modified: 2024-12-23 16:42
# Modified By  : EastsunW
# -------------
# Description  : 画SNP数量随样本数量增加的变化
# -------------

# 加载必要的包
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


# 读取数据
snp_increase <- fread(
    "data/cohort_SNP_sample_increase.txt",
    sep = "\t",
    header = TRUE
) %>%
    pivot_longer(-n_sample, names_to = "variant", values_to = "count") %>%
    group_by(variant) %>%
    mutate(max = max(count)) %>%
    mutate(percent = count / max)
sv_increase <- fread(
    "data/cohort_SV_sample_increase.txt",
    sep = "\t",
    header = TRUE
) %>%
    pivot_longer(-n_sample, names_to = "variant", values_to = "count") %>%
    group_by(variant) %>%
    mutate(max = max(count)) %>%
    mutate(percent = count / max)
sv_all <- fread(
    "data/cohort_SV_sample_increase.txt",
    sep = "\t",
    header = TRUE
) %>%
    pivot_longer(-n_sample, names_to = "variant", values_to = "count") %>%
    group_by(n_sample) %>%
    summarise(count = sum(count)) %>%
    mutate(variant = "SV_all") %>%
    mutate(max = max(count)) %>%
    mutate(percent = count / max)
merged_df <- rbind(
    snp_increase,
    sv_all,
    sv_increase
) %>%
    mutate_at(
        vars("variant"),
        ~ factor(., levels = c("SNP", "InDel", "SV_all", "INS", "DEL", "INV", "DUP", "BND"))
    ) %>%
    group_by(variant) %>%
    mutate(max = max(count)) %>%
    mutate(percent = count / max) %>%
    mutate(half_sample = percent[n_sample == 74])

pdf(
    "results/plot/cohortVariant_sample_increase.pdf",
    width = 9,
    height = 3.5,
    onefile = TRUE
)
ggplot(
    data = merged_df,
    mapping = aes(
        x = n_sample,
        y = percent,
        color = variant
    )
) +
    geom_line(linewidth = 1) +
    geom_segment(
        data = merged_df %>% filter(n_sample == 30),
        mapping = aes(
            x = 1,
            xend = n_sample,
            y = percent,
            yend = percent
        ),
        linetype = "dashed",
        color = "red",
        linewidth = 0.5
    ) +
    geom_segment(
        data = merged_df %>% filter(n_sample == 30),
        mapping = aes(
            x = n_sample,
            xend = n_sample,
            y = 0,
            yend = percent
        ),
        linetype = "dashed",
        color = "red",
        linewidth = 0.5
    ) +
    geom_point(
        data = merged_df %>% filter(n_sample == 30),
        mapping = aes(
            x = n_sample,
            y = percent,
            color = variant
        ),
        size = 3,
        color = "black"
    ) +
    scale_x_continuous(
        expand = c(0, 0),
        limits = c(1, 148),
        breaks = c(1, 30, 74, 148)
    ) +
    scale_y_continuous(
        expand = c(0, 0),
        limits = c(0, 1),
        breaks = c(0, 0.5, 0.75, 1),
        labels = scales::label_number(scale = 100, suffix = "%")
    ) +
    scale_fill_manual(values = custom_colors) +
    facet_wrap(~variant, scales = "free_y", nrow = 2) +
    guides(col = guide_legend(nrow = 1, byrow = TRUE)) +
    wdy_theme(
        axis.line = element_blank(),
        legend.position = "none",
        legend.direction = "horizontal",
        legend.position.inside = c(0.99, 0.25),
        legend.justification = c(0.5, 0.5),
        legend.background = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size = 10),
    ) +
    labs(
        x = "Number of Samples",
        y = "Number of variant",
        fill = "Type"
    )
dev.off()
