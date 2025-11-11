#!/usr/bin/env Rscript

# -------------
# FileName     : 09_circos.R
# Author       : Tian Wu fabien012390wt@163.com
# Modified By  : Tian Wu
# -------------
# Description  : Generate circos plots for trans-QTL hotspot regions visualization
# -------------

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(ggrepel))
suppressPackageStartupMessages(library(patchwork))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(parallel))
suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(vroom))
suppressPackageStartupMessages(library(pryr))
suppressPackageStartupMessages(library(circlize))

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
options(stringsAsFactors = FALSE)
mc <- getOption('mc.cores', 36)

circos_plot <- function(bed, filename) {
    len <- dim(bed)[1]
    pdf(filename, width = 10, height = 10)
    circos.initializeWithIdeogram(species = "hg38", plotType = NULL)
    
    circos.track(
      ylim = c(0, 1), 
      panel.fun = function(x, y) {
        chr <- CELL_META$sector.index
        chr <- gsub("^chr", "", chr)
        xlim <- CELL_META$xlim
        ylim <- CELL_META$ylim
        circos.rect(xlim[1], 0, xlim[2], 1, col = c("#E7DBD2"))
        circos.text(
          mean(xlim), mean(ylim), chr, cex = 1.5,
          col = "black", facing = "inside",
          niceFacing = TRUE
        )
      }, 
      track.height = 0.06, bg.border = NA
    )
    
    circos.genomicLabels(
      bed[1, , drop = FALSE], labels.column = 4, side = "outside",
      col = "#FA7F73", line_col = "#FA7F73",
      labels_height = max(strwidth(bed$gene)) + 0.03,
      padding = 0.1, cex = 1.2
    )
    
    circos.genomicIdeogram(track.height = mm_h(4))
    
    links <- data.frame(from = rep(1, len-1), to = seq(2, len))
    circos.genomicLink(
      bed[links$from, ], bed[links$to, ], 
      col = adjustcolor("#FA7F73", alpha.f = 0.5),
      border = NA, lwd = 2
    )
    
    dev.off()
}

gene_df <- data.frame(
    chr = "19", start = 19900913, end = 19963464, gene = "ZNF93"
)

bed_df <- read.table("promoter_region.bed", header = FALSE, sep = "\t")

df <- data.frame(
    chr = c(gene_df$chr, bed_df$V1),
    start = c(gene_df$start, bed_df$V2),
    end = c(gene_df$end, bed_df$V3),
    gene = c(gene_df$gene, seq(1, nrow(bed_df)))
)

df$chr <- ifelse(grepl("^chr", df$chr), df$chr, paste0("chr", df$chr))
df$end <- df$end + 1000000

circos_plot(df, "circos.pdf")

