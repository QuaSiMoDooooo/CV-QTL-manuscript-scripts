#!/usr/bin/env Rscript

# -------------
# FileName     : cis_eQTL_case.R
# Author       : Tian Wu fabien012390wt@163.com
# Modified By  : Tian Wu
# -------------
# Description  : Analyze cis-eQTL characteristics by variant type (SNP, MNV, InDel, SV)
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

setwd("/home/wtian/project/HZAU_cohort_meth/wdy_assist_rerun/01_quick_stats_for_paper/QTL_reg_stats")
options(stringsAsFactors = FALSE)
mc <- getOption('mc.cores', 36)

qtl_dirs <- list.dirs("/home/wtian/project/HZAU_cohort_meth/wdy_result/phenotype_qtl_results_bim3_removaNA/", recursive = FALSE)
qtl_dirs <- qtl_dirs[!grepl("biochemistry|common", qtl_dirs) & grepl("expression", qtl_dirs)]

cis_eqtl_paths <- list.files(qtl_dirs, pattern = "cis.filtered.txt.gz", full.names = TRUE, recursive = TRUE)
var_types <- c("SNP", "MNV", "InDel", "SV")

read_func <- function(v) {
    vartype <- var_types[v]
    path <- cis_eqtl_paths[grepl(vartype, cis_eqtl_paths)]
    tmp <- as.data.frame(fread(path))
    tmp$vartype <- vartype
    tmp <- select(tmp, vartype, everything())
    return(tmp)
}

all_cis_eqtl_list <- lapply(1:length(var_types), read_func)
all_cis_eqtl <- do.call(rbind, all_cis_eqtl_list)

get_varianttype_specific_phenotypes <- function(df) {
    if (!all(c("vartype", "Phenotype") %in% colnames(df))) {
        stop("Data frame must contain 'vartype' and 'Phenotype' columns")
    }
    pheno_variant_table <- aggregate(vartype ~ Phenotype, data = df, function(x) unique(x))
    single_varianttype_phenos <- pheno_variant_table$Phenotype[sapply(pheno_variant_table$vartype, length) == 1]
    multi_varianttype_phenos <- pheno_variant_table$Phenotype[sapply(pheno_variant_table$vartype, length) > 1]
    return(list(
        single_varianttype_phenos = single_varianttype_phenos,
        multi_varianttype_phenos = multi_varianttype_phenos
    ))
}

res <- get_varianttype_specific_phenotypes(all_cis_eqtl)
single_vartype_phe <- res$single_varianttype_phenos
multi_vartype_phe <- res$multi_varianttype_phenos

single_vartype_phe_qtl <- all_cis_eqtl %>% filter(Phenotype %in% single_vartype_phe)
multi_vartype_phe_qtl <- all_cis_eqtl %>% filter(Phenotype %in% multi_vartype_phe)

single_vartype_phe_qtl_dedup <- single_vartype_phe_qtl %>% distinct(Phenotype, .keep_all = TRUE)

library(ggplot2)

df <- data.frame(
    Category = "cis-eQTL Genes",
    Type = factor(c("Only SNV", "Only MNV", "Only InDel", "Only SV", "Multiple variants"),
                levels = c("Only SNV", "Only MNV", "Only InDel", "Only SV", "Multiple variants")),
    Count = c(572, 17, 16, 10, 2829 - (572 + 17 + 16 + 10))
)

df$Percent <- df$Count / sum(df$Count)

custom_colors <- c(
    "Only SNV" = "#FB7F72",
    "Only MNV" = "#FFBD80",
    "Only InDel" = "#80C5BF",
    "Only SV" = "#82B1D1",
    "Multiple variants" = "#E7DBD2"
)

source("https://gitee.com/eastsunw/personal_code_notebook/raw/master/plot_tools/wdy_theme.r")

p <- ggplot(df, aes(x = Category, y = Percent, fill = Type)) +
    geom_bar(stat = "identity", width = 0.5) +
    scale_fill_manual(values = custom_colors) +
    scale_y_continuous(
        labels = scales::percent_format(accuracy = 1),
        limits = c(0, 1),
        breaks = seq(0, 1, by = 0.2),
        expand = expansion(add = c(0, 0.02))
    ) +
    expand_limits(y = c(0)) +
    labs(x = NULL, y = "Percentage (%)", fill = "Variant Type") +
    wdy_theme() +
    theme(legend.position = "right")

ggsave("cis-eQTL_variant_type.pdf", p, width = 4, height = 3, dpi = 300)

single_vartype_phe_qtl$abs_beta <- abs(single_vartype_phe_qtl$Beta)
mycolor <- c("#F37C72", "#F9B27C", "#82C6C0", "#78A7C4")

pdf("cis_eQTL_abs_beta_boxplot.pdf", width = 4, height = 3)
ggpubr::ggboxplot(single_vartype_phe_qtl, x = "vartype", y = "abs_beta", color = "vartype", palette = mycolor) +
    labs(x = NULL, y = "The absolute value of Beta") +
    theme(
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 10),
        legend.position = "none"
    )
dev.off()

pdf("cis_eQTL_abs_beta_horizon_boxplot.pdf", width = 3, height = 4)
ggpubr::ggboxplot(
    single_vartype_phe_qtl,
    x = "vartype",
    y = "abs_beta",
    color = "vartype",
    palette = mycolor
) +
    labs(x = NULL, y = "The absolute value of Beta") +
    theme(
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 10),
        legend.position = "none"
    ) +
    coord_flip()
dev.off()