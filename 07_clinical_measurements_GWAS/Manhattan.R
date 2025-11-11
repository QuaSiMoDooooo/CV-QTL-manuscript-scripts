#!/usr/bin/env Rscript

# -------------
# FileName     : Manhattan.R
# Author       : Tian Wu fabien012390wt@163.com
# Modified By  : Tian Wu
# -------------
# Description  : Generate Manhattan plots for clinical measurements GWAS analysis
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

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
options(stringsAsFactors = FALSE)
mc <- getOption('mc.cores', 36)

gwas_dir = "/home/wtian/project/HZAU_cohort_meth/wdy_assist_rerun/consolidation_data_and_res_softlinks/all_gwas"
gwas_paths = list.files(gwas_dir, full.names = TRUE, recursive = FALSE, pattern = "*.tsv")

gwas = do.call(
  rbind,
  lapply(
    gwas_paths,
    function(x) {
      gwas = fread(x, header = TRUE, sep = "\t")
      gwas$variant_type = str_split(basename(x), "\\.", simplify = TRUE)[, 1]
      return(gwas)
    }
  )
)

gwas$variant_type = str_replace_all(gwas$variant_type, "SNP", "SNV")
gwas$variant_type = factor(gwas$variant_type, levels = c("SNV", "MNV", "InDel", "SV"))
gwas$Indicator = str_replace_all(gwas$Indicator, ".PHENO1", "")
colnames(gwas)[1] = "CHROM"

variant_colors <- c(
  SNV = "#F37C73",
  MNV = "#F9B180",
  InDel = "#8FCFCA",
  SV = "#82B1D3"
)

gwas[, is_shared := duplicated(ID) | duplicated(ID, fromLast = TRUE)]
gwas[, logP := -log10(P)]
setorder(gwas, CHROM, POS)
gwas[, chr := factor(CHROM, levels = 1:22)]
gwas[, pos_index := .I]
axis_dt <- gwas[, .(center = mean(pos_index)), by = chr]

biochem = fread("/home/wtian/project/HZAU_cohort_meth/wdy_assist/14_cases/consolidation_res/Biochemistry.all.txt")
Indicator_pathology = biochem$Indicator[biochem$Category == "Pathology"]
Indicator_pathology = str_replace_all(Indicator_pathology, " |/", "_")
gwas_pathology = gwas[gwas$Indicator %in% Indicator_pathology, ]
top_gwas_pathology = gwas_pathology[order(logP), .SD[1:min(.N, 10)], by = .(Indicator, variant_type)]

source("https://gitee.com/eastsunw/personal_code_notebook/raw/master/plot_tools/wdy_theme.r")
p <- ggplot(top_gwas_pathology, aes(x = pos_index, y = logP)) +
    geom_point(
    aes(color = variant_type, shape = is_shared),
    size = 3, alpha = 0.8, stroke = 0
    ) +
    scale_color_manual(values = variant_colors) +
    scale_shape_manual(values = c("FALSE" = 16, "TRUE" = 17)) +
    scale_x_continuous(
    breaks = axis_dt$center,
    labels = axis_dt$chr,
    expand = expansion(add = c(0.01, 0.01))
    ) +
    scale_y_continuous(limits = c(3.5, ceiling(max(top_gwas_pathology$logP)))) + 
    labs(
    x = "Chromosome",
    y = expression(-log[10](italic(P))),
    color = "Variant Type",
    shape = "Shared Across Phenotypes"
    ) +
    guides(
    color = guide_legend(override.aes = list(size = 4)),
    shape = guide_legend(override.aes = list(size = 4))
    ) + 
    wdy_theme() +
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1, size = 9)
    )

ggsave(p, filename = "manhattan.pdf", width = 12, height = 4, dpi = 300)