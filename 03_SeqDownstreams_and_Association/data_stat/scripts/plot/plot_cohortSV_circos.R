#! Rscript
# -------------
# FileName     : plot_cohortVariant_circos
# Author       : EastsunW eastsunw@foxmail.com
# Create at    : 2024-11-26 21:30
# Last Modified: 2024-12-17 16:01
# Modified By  : EastsunW
# -------------
# Description  : 绘制变异的全基因组分布circos图
# -------------

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(circlize))
suppressPackageStartupMessages(library(ComplexHeatmap))
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

gene_density <- fread("results/stat/circos_density/gene_density.txt")
repeat_density <- fread("results/stat/circos_density/repeat_density.txt")
sv_all_density <- fread("results/stat/circos_density/SV_all_density.txt")
sv_ins_density <- fread("results/stat/circos_density/SV_INS_density.txt")
sv_del_density <- fread("results/stat/circos_density/SV_DEL_density.txt")
sv_inv_density <- fread("results/stat/circos_density/SV_INV_density.txt")
sv_dup_density <- fread("results/stat/circos_density/SV_DUP_density.txt")

pdf(
    "results/plot/cohortSV_circos.pdf",
    width = 8,
    height = 8,
    onefile = TRUE
)
circos.par(start.degree = 90)
circos.initializeWithIdeogram(
    chromosome.index = paste0("chr", c(as.character(1:22), "X")),
    plotType = c("ideogram", "labels"),
    species = "hg38"
)
circos.genomicHeatmap(
    gene_density,
    heatmap_height = 0.05,
    side = "inside",
    col = colorRamp2(
        seq(0, 1, 0.25),
        c("#ffffff", "#b3eba2", "#40ca9c", "#1b63b4", "#002153")
    ),
    connection_height = NULL
)
circos.genomicHeatmap(
    repeat_density,
    heatmap_height = 0.05,
    side = "inside",
    col = colorRamp2(
        seq(0, 1, 0.25),
        c("#ffffff", "#e8e9a7", "#f09c3d", "#c01c16", "#4b0000")
    ),
    connection_height = NULL
)
circos.genomicTrack(
    sv_all_density,
    track.height = 0.06,
    bg.border = NA,
    ylim = c(0, 1),
    panel.fun = function(region, value, ...) {
        circos.genomicLines(
            region, value,
            lwd = 1.5,
            type = "l", area = TRUE,
            border = custom_colors[5],
            col = custom_colors[5]
        )
    }
)
circos.genomicTrack(
    sv_ins_density,
    track.height = 0.06,
    bg.border = NA,
    ylim = c(0, 1),
    panel.fun = function(region, value, ...) {
        circos.genomicLines(
            region, value,
            lwd = 1.5,
            type = "l", area = TRUE,
            border = custom_colors[1],
            col = custom_colors[1]
        )
    }
)
circos.genomicTrack(
    sv_del_density,
    track.height = 0.06,
    bg.border = NA,
    ylim = c(0, 1),
    panel.fun = function(region, value, ...) {
        circos.genomicLines(
            region, value,
            lwd = 1.5,
            type = "l", area = TRUE,
            border = custom_colors[2],
            col = custom_colors[2]
        )
    }
)
circos.genomicTrack(
    sv_inv_density,
    track.height = 0.06,
    bg.border = NA,
    ylim = c(0, 1),
    panel.fun = function(region, value, ...) {
        circos.genomicLines(
            region, value,
            lwd = 1.5,
            type = "l", area = TRUE,
            border = custom_colors[3],
            col = custom_colors[3]
        )
    }
)
circos.genomicTrack(
    sv_dup_density,
    track.height = 0.06,
    bg.border = NA,
    ylim = c(0, 1),
    panel.fun = function(region, value, ...) {
        circos.genomicLines(
            region, value,
            lwd = 1.5,
            type = "l", area = TRUE,
            border = custom_colors[4],
            col = custom_colors[4]
        )
    }
)
legend_variant <- Legend(
    labels = c("SV-all", "SV-INS", "SV-DEL", "SV-INV", "SV-DUP"),
    type = "grid",
    border = NULL,
    background = NULL,
    legend_gp = gpar(col = custom_colors[c(5, 1:4)], fill = custom_colors[c(5, 1:4)]),
    title_position = "topleft",
    title = "Variant type",
    title_gp = gpar(fontsize = 12, fontface = "bold"),
    labels_gp = gpar(fontsize = 10),
    ncol = 2,
    by_row = TRUE
)
legend_gene <- Legend(
    col_fun = colorRamp2(
        c(0, 0.2, 0.5, 0.8, 1),
        c("#ffffff", "#fff27c", "#3EBB72", "#2B708D", "#45085B")
    ),
    direction = "horizontal",
    title = "Gene density"
)
legend_repeat <- Legend(
    col_fun = colorRamp2(
        c(0, 0.2, 0.5, 0.8, 1),
        c("#ffffff", "#ffcdbb", "#ca1a1a", "#570b0b", "#000000")
    ),
    direction = "horizontal",
    title = "Repeat density"
)
legend_list <- packLegend(legend_gene, legend_repeat, legend_variant)
draw(legend_list)
dev.off()
