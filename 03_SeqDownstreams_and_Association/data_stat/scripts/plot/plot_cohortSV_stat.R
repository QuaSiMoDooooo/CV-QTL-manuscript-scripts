#! Rscript
# -------------
# FileName     : plot_SV_cohort_count
# Author       : EastsunW eastsunw@foxmail.com
# Create at    : 2024-11-25 23:23
# Last Modified: 2024-12-14 16:11
# Modified By  : EastsunW
# -------------
# Description  :
# -------------

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(grid))
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

# 导入数据
cohort_sv_stat <- fread(
    "results/bed/SV.all.bed",
    header = FALSE,
    col.names = c("chr", "start", "end", "ID", "MAF", "sv_type", "sv_len", "support")
) %>%
    mutate(freq = ifelse(MAF > 0.05, "common", "rare"))

# 不同类型的SV长度密度和分布，只展示3000bp以内的INS和DEL，8000以内的INV和DUP

# INS
density_ins <- density(cohort_sv_stat$sv_len[cohort_sv_stat$sv_type == "INS" & cohort_sv_stat$sv_len < 10000])
peaks_ins <- pracma::findpeaks(density_ins$y)
alu_ins <- c(
    density_ins$x[peaks_ins[which.min(abs(density_ins$x[peaks_ins[, 2]] - 300)), 2]],
    density_ins$y[peaks_ins[which.min(abs(density_ins$x[peaks_ins[, 2]] - 300)), 2]]
)
line_ins <- c(
    density_ins$x[peaks_ins[which.min(abs(density_ins$x[peaks_ins[, 2]] - 6000)), 2]],
    density_ins$y[peaks_ins[which.min(abs(density_ins$x[peaks_ins[, 2]] - 6000)), 2]]
)
length_ins <- cohort_sv_stat %>%
    select(ID, sv_type, sv_len) %>%
    filter(sv_type == "INS")
median_ins <- c(
    median(length_ins$sv_len),
    max(density_ins$y) * 0.75
)
full_plot_ins <- ggplot(
    data = length_ins,
    mapping = aes(
        x = sv_len,
        fill = sv_type
    )
) +
    geom_density() +
    geom_vline(
        xintercept = median_ins[1],
        color = "black",
        linetype = "dashed"
    ) +
    geom_text(
        x = median_ins[1],
        y = median_ins[2],
        label = median_ins[1],
        color = "black",
        vjust = -0.5,
        hjust = -0.05,
        size = 10 / .pt
    ) +
    geom_vline(
        xintercept = alu_ins[1],
        color = "black",
        linetype = "dashed"
    ) +
    geom_text(
        x = alu_ins[1],
        y = alu_ins[2],
        label = "Alu",
        color = "black",
        vjust = -0.5,
        hjust = -0.05,
        size = 10 / .pt
    ) +
    geom_vline(
        xintercept = line_ins[1],
        color = "black",
        linetype = "dashed"
    ) +
    geom_text(
        x = line_ins[1],
        y = line_ins[2],
        label = "LINE",
        color = "black",
        vjust = -0.5,
        hjust = -0.05,
        size = 10 / .pt
    ) +
    scale_x_continuous(expand = c(0, 0), limits = c(0, 10000)) +
    scale_y_continuous(expand = c(0.1, 0), n.breaks = 3) +
    wdy_theme(
        base_size = 12,
        base_line_size = 1,
        legend.position = "none",
        axis.line = element_blank(),
        axis.title = element_blank(),
        panel.spacing.x = unit(5, "pt"),
        strip.background = element_rect(fill = "#ebf1f5", linewidth = NA)
    ) +
    scale_fill_manual(values = custom_colors[1]) +
    labs(
        x = "SV length (bp)",
        y = "Density"
    )
main_plot_ins <- full_plot_ins +
    coord_cartesian(xlim = c(50, 2000))
sub_plot_ins <- full_plot_ins +
    coord_cartesian(
        xlim = c(line_ins[1] - 300, line_ins[1] + 300),
        ylim = c(0, line_ins[2] * 1.5),
    ) +
    scale_x_continuous(
        n.breaks = 3,
        breaks = c(floor(line_ins[1] / 100 - 2) * 100, round(line_ins[1], 0), ceiling(line_ins[1] / 100 + 2) * 100),
        expand = c(0, 0),
        limits = c(0, 10000)
    )
final_plot_ins <- main_plot_ins +
    annotation_custom(
        grob = ggplotGrob(sub_plot_ins),
        xmin = 800, xmax = 1800,
        ymin = 0.05 * max(density_ins$y), ymax = 0.99 * max(density_ins$y)
    ) +
    facet_wrap(~sv_type, ncol = 1, scales = "free", strip.position = "right")

# DEL
density_del <- density(cohort_sv_stat$sv_len[cohort_sv_stat$sv_type == "DEL" & cohort_sv_stat$sv_len < 10000])
peaks_del <- pracma::findpeaks(density_del$y)
alu_del <- c(
    density_del$x[peaks_del[which.min(abs(density_del$x[peaks_del[, 2]] - 300)), 2]],
    density_del$y[peaks_del[which.min(abs(density_del$x[peaks_del[, 2]] - 300)), 2]]
)
line_del <- c(
    density_del$x[peaks_del[which.min(abs(density_del$x[peaks_del[, 2]] - 6000)), 2]],
    density_del$y[peaks_del[which.min(abs(density_del$x[peaks_del[, 2]] - 6000)), 2]]
)
length_del <- cohort_sv_stat %>%
    select(ID, sv_type, sv_len) %>%
    filter(sv_type == "DEL")
median_del <- c(
    median(length_del$sv_len),
    max(density_del$y) * 0.75
)
full_plot_del <- ggplot(
    data = length_del,
    mapping = aes(
        x = sv_len,
        fill = sv_type
    )
) +
    geom_density() +
    geom_vline(
        xintercept = median_del[1],
        color = "black",
        linetype = "dashed"
    ) +
    geom_text(
        x = median_del[1],
        y = median_del[2],
        label = median_del[1],
        color = "black",
        vjust = -0.5,
        hjust = -0.05,
        size = 10 / .pt
    ) +
    geom_vline(
        xintercept = alu_del[1],
        color = "black",
        linetype = "dashed"
    ) +
    geom_text(
        x = alu_del[1],
        y = alu_del[2],
        label = "Alu",
        color = "black",
        vjust = -0.5,
        hjust = -0.05,
        size = 10 / .pt
    ) +
    geom_vline(
        xintercept = line_del[1],
        color = "black",
        linetype = "dashed"
    ) +
    geom_text(
        x = line_del[1],
        y = line_del[2],
        label = "LINE",
        color = "black",
        vjust = -0.5,
        hjust = -0.05,
        size = 10 / .pt
    ) +
    scale_x_continuous(expand = c(0, 0), limits = c(0, 10000)) +
    scale_y_continuous(expand = c(0.1, 0), n.breaks = 3) +
    wdy_theme(
        base_size = 12,
        base_line_size = 1,
        legend.position = "none",
        axis.line = element_blank(),
        axis.title = element_blank(),
        panel.spacing.x = unit(5, "pt"),
        strip.background = element_rect(fill = "#ebf1f5", linewidth = NA)
    ) +
    scale_fill_manual(values = custom_colors[2]) +
    labs(
        x = "SV length (bp)",
        y = "Density"
    )
main_plot_del <- full_plot_del +
    coord_cartesian(xlim = c(50, 2000))
sub_plot_del <- full_plot_del +
    coord_cartesian(
        xlim = c(line_del[1] - 300, line_del[1] + 300),
        ylim = c(0, line_del[2] * 1.5),
    ) +
    scale_x_continuous(
        n.breaks = 3,
        breaks = c(floor(line_del[1] / 100 - 2) * 100, round(line_del[1], 0), ceiling(line_del[1] / 100 + 2) * 100),
        expand = c(0, 0),
        limits = c(0, 10000)
    )
final_plot_del <- main_plot_del +
    annotation_custom(
        grob = ggplotGrob(sub_plot_del),
        xmin = 800, xmax = 1800,
        ymin = 0.05 * max(density_del$y), ymax = 0.99 * max(density_del$y)
    ) +
    facet_wrap(~sv_type, ncol = 1, scales = "free", strip.position = "right")

# inv
density_inv <- density(cohort_sv_stat$sv_len[cohort_sv_stat$sv_type == "INV" & cohort_sv_stat$sv_len < 100000])
length_inv <- cohort_sv_stat %>%
    select(ID, sv_type, sv_len) %>%
    filter(sv_type == "INV")
median_inv <- c(
    median(length_inv$sv_len),
    max(density_inv$y) * 0.75
)
plot_inv <- ggplot(
    data = length_inv,
    mapping = aes(
        x = sv_len,
        fill = sv_type
    )
) +
    geom_density() +
    geom_vline(
        xintercept = median_inv[1],
        color = "black",
        linetype = "dashed"
    ) +
    geom_text(
        x = median_inv[1],
        y = median_inv[2],
        label = median_inv[1],
        color = "black",
        vjust = 1,
        hjust = -0.15,
        size = 12 / .pt
    ) +
    scale_x_continuous(expand = c(0, 0), limits = c(0, 100000)) +
    scale_y_continuous(expand = c(0.1, 0), n.breaks = 3) +
    coord_cartesian(xlim = c(50, 100000)) +
    wdy_theme(
        base_size = 12,
        base_line_size = 1,
        legend.position = "none",
        axis.line = element_blank(),
        axis.title = element_blank(),
        panel.spacing.x = unit(5, "pt"),
        strip.background = element_rect(fill = "#ebf1f5", linewidth = NA)
    ) +
    scale_fill_manual(values = custom_colors[3]) +
    labs(
        x = "SV length (bp)",
        y = "Density"
    ) +
    facet_wrap(~sv_type, ncol = 1, scales = "free", strip.position = "right")

# DUP
density_dup <- density(cohort_sv_stat$sv_len[cohort_sv_stat$sv_type == "DUP" & cohort_sv_stat$sv_len < 100000])
length_dup <- cohort_sv_stat %>%
    select(ID, sv_type, sv_len) %>%
    filter(sv_type == "DUP")
median_dup <- c(
    median(length_dup$sv_len),
    max(density_dup$y) * 0.75
)
plot_dup <- ggplot(
    data = length_dup,
    mapping = aes(
        x = sv_len,
        fill = sv_type
    )
) +
    geom_density() +
    geom_vline(
        xintercept = median_dup[1],
        color = "black",
        linetype = "dashed"
    ) +
    geom_text(
        x = median_dup[1],
        y = median_dup[2],
        label = median_dup[1],
        color = "black",
        vjust = 1,
        hjust = -0.15,
        size = 12 / .pt
    ) +
    scale_x_continuous(expand = c(0, 0), limits = c(0, 100000)) +
    scale_y_continuous(expand = c(0.1, 0), n.breaks = 3) +
    coord_cartesian(xlim = c(50, 100000)) +
    wdy_theme(
        base_size = 12,
        base_line_size = 1,
        legend.position = "none",
        axis.line = element_blank(),
        axis.title = element_blank(),
        panel.spacing.x = unit(5, "pt"),
        strip.background = element_rect(fill = "#ebf1f5", linewidth = NA)
    ) +
    scale_fill_manual(values = custom_colors[4]) +
    labs(
        x = "SV length (bp)",
        y = "Density"
    ) +
    facet_wrap(~sv_type, ncol = 1, scales = "free", strip.position = "right")

pdf("results/plot/cohortSV_length.pdf", width = 5, height = 6, onefile = TRUE)
merged_plot <- ggarrange(
    final_plot_ins, final_plot_del, plot_inv, plot_dup,
    ncol = 1,
    align = "v"
)
annotate_figure(
    merged_plot,
    left = textGrob("Density", rot = 90, gp = gpar(color = "black", fontsize = 12, fontface = "bold")),
    bottom = textGrob("SV length (bp)", gp = gpar(color = "black", fontsize = 12, fontface = "bold"))
)
dev.off()

# SV的样本支持数
sv_avgsupport <- cohort_sv_stat %>%
    group_by(sv_type) %>%
    summarise(mean_support = mean(support)) %>%
    distinct(sv_type, .keep_all = TRUE) %>%
    arrange(desc(mean_support)) %>%
    mutate_at(
        vars("sv_type"),
        ~ factor(., levels = sv_type)
    )
pdf("results/plot/cohortSV_meanSupport.pdf", width = 3.5, height = 3, onefile = TRUE)
ggplot(
    data = sv_avgsupport,
    mapping = aes(
        x = sv_type,
        y = mean_support,
        fill = sv_type
    )
) +
    geom_bar(stat = "identity", width = 0.8, color = "black") +
    geom_text(
        aes(label = round(mean_support, 0)),
        size = 12 / .pt,
        color = "black",
        vjust = -1
    ) +
    wdy_theme(
        base_size = 12,
        base_line_size = 1,
        legend.position = "none",
        axis.line = element_blank(),
        strip.background = element_rect(fill = "#ebf1f5", linewidth = NA)
    ) +
    scale_fill_manual(values = custom_colors) +
    scale_y_continuous(
        expand = c(0, 0),
        limits = c(0, 45)
    ) +
    labs(
        x = "SV type",
        y = "Agerage supprted samples"
    )
dev.off()
