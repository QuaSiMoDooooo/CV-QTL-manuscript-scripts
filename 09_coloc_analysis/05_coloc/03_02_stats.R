#!/usr/bin/env Rscript

# -------------
# FileName     : 03_02_stats.R
# Author       : Tian Wu fabien012390wt@163.com
# Modified By  : Tian Wu
# -------------
# Description  : Statistical results
# -------------

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(parallel))
suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(vroom))
# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
options(stringsAsFactors = FALSE)
# options(scipen = 999)
mc <- getOption('mc.cores', 36)
setwd("/home/wtian/project/HZAU_cohort_meth/wdy_assist_rerun/11_coloc/05_coloc")
getwd()

# gwas
gwas_qtl_coloc = fread("03_gwas_qtl_colocResultsig.tsv")
# head(gwas_qtl_coloc)
# dim(gwas_qtl_coloc)
# trait_phe_uniq_counts = length(unique(paste(gwas_qtl_coloc$trait, gwas_qtl_coloc$phe, sep = "_")))
# trait_phe_uniq_counts
# trait_counts = length(unique(gwas_qtl_coloc$trait))
# trait_counts
# phe_counts = length(unique(gwas_qtl_coloc$phe))
# phe_counts
# gwas_qtl_coloc = gwas_qtl_coloc %>% mutate(trait_phe = paste(trait, phe, sep = "_"))
# # SV
# gwas_qtl_coloc_SV = gwas_qtl_coloc %>% filter(var_type == "SV")
# trait_phe_uniq_counts_SV = length(unique(paste(gwas_qtl_coloc_SV$trait, gwas_qtl_coloc_SV$phe, sep = "_")))
# trait_phe_uniq_counts_SV
# trait_counts_SV = length(unique(gwas_qtl_coloc_SV$trait))
# trait_counts_SV
# trait_SV = unique(gwas_qtl_coloc_SV$trait)
# phe_counts_SV = length(unique(gwas_qtl_coloc_SV$phe))
# phe_counts_SV
# phe_SV = unique(gwas_qtl_coloc_SV$phe)
# # SNP
# gwas_qtl_coloc_SNP = gwas_qtl_coloc %>% filter(var_type == "SNP")
# trait_phe_uniq_counts_SNP = length(unique(paste(gwas_qtl_coloc_SNP$trait, gwas_qtl_coloc_SNP$phe, sep = "_")))
# trait_phe_uniq_counts_SNP
# trait_counts_SNP = length(unique(gwas_qtl_coloc_SNP$trait))
# trait_counts_SNP
# trait_SNP = unique(gwas_qtl_coloc_SNP$trait)
# phe_counts_SNP = length(unique(gwas_qtl_coloc_SNP$phe))
# phe_counts_SNP
# phe_SNP = unique(gwas_qtl_coloc_SNP$phe)
# # levels
# tb_t = table(trait_SV %in% trait_SNP)
# tb_p = table(phe_SV %in% phe_SNP)
# trait_level  = round(tb_t[2] / length(trait_SV) * 100, 2)
# phe_level = round(tb_p[2] / length(phe_SV) * 100, 2)
# trait_level
# phe_level
# df_gwas = data.frame(
#     trait_phe_uniq_counts = trait_phe_uniq_counts,
#     trait_counts = trait_counts,
#     phe_counts = phe_counts,
#     trait_phe_uniq_counts_SV = trait_phe_uniq_counts_SV,
#     trait_counts_SV = trait_counts_SV,
#     phe_counts_SV = phe_counts_SV,
#     trait_phe_uniq_counts_SNP = trait_phe_uniq_counts_SNP,
#     trait_counts_SNP = trait_counts_SNP,
#     phe_counts_SNP = phe_counts_SNP,
#     trait_level = trait_level,
#     phe_level = phe_level
# )
# df_gwas
# write.table(df_gwas, "03_gwas_qtl_colocResultsig_stats.tsv", sep = "\t", quote = F, row.names = F)

head(gwas_qtl_coloc)
trait_sv = unique(gwas_qtl_coloc$trait[gwas_qtl_coloc$var_type == "SV"])
trait_snv = unique(gwas_qtl_coloc$trait[gwas_qtl_coloc$var_type == "SNP"])
table(trait_sv %in% trait_snv)
length(unique(gwas_qtl_coloc$trait))
round(table(trait_sv %in% trait_snv)[1] / length(unique(gwas_qtl_coloc$trait)) * 100, 2)

tarit_common = unique(gwas_qtl_coloc$trait[gwas_qtl_coloc$var_type == "SV"]) %>% intersect(unique(gwas_qtl_coloc$trait[gwas_qtl_coloc$var_type == "SNP"]))
counts = 0
for ( t in 1:length(tarit_common)) {
    tarit_sv_asso = gwas_qtl_coloc$phe[gwas_qtl_coloc$var_type == "SV" & gwas_qtl_coloc$trait == tarit_common[t]]
    tarit_snv_asso = gwas_qtl_coloc$phe[gwas_qtl_coloc$var_type == "SNP" & gwas_qtl_coloc$trait == tarit_common[t]]
    if (any(! tarit_sv_asso %in% tarit_snv_asso)) {
        counts = counts + 1
    }
}
counts
round(counts / length(tarit_common) * 100, 2)

# eqtl
eqtl_coloc = fread("03_eqtl_qtl_colocResultsig.tsv")
# head(eqtl_coloc)
# dim(eqtl_coloc)
# gene_phe_uniq_counts = length(unique(paste(eqtl_coloc$gene, eqtl_coloc$phe, sep = "_")))
# gene_phe_uniq_counts
# gene_counts = length(unique(eqtl_coloc$gene))
# gene_counts
# phe_counts = length(unique(eqtl_coloc$phe))
# phe_counts
# eqtl_coloc = eqtl_coloc %>% mutate(gene_phe = paste(gene, phe, sep = "_"))
# # SV
# eqtl_coloc_SV = eqtl_coloc %>% filter(var_type == "SV")
# gene_phe_uniq_counts_SV = length(unique(paste(eqtl_coloc_SV$gene, eqtl_coloc_SV$phe, sep = "_")))
# gene_phe_uniq_counts_SV
# gene_counts_SV = length(unique(eqtl_coloc_SV$gene))
# gene_counts_SV
# gene_SV = unique(eqtl_coloc_SV$gene)
# phe_counts_SV = length(unique(eqtl_coloc_SV$phe))
# phe_counts_SV
# phe_SV = unique(eqtl_coloc_SV$phe)
# # SNP
# eqtl_coloc_SNP = eqtl_coloc %>% filter(var_type == "SNP")
# gene_phe_uniq_counts_SNP = length(unique(paste(eqtl_coloc_SNP$gene, eqtl_coloc_SNP$phe, sep = "_")))
# gene_phe_uniq_counts_SNP
# gene_counts_SNP = length(unique(eqtl_coloc_SNP$gene))
# gene_counts_SNP
# gene_SNP = unique(eqtl_coloc_SNP$gene)
# phe_counts_SNP = length(unique(eqtl_coloc_SNP$phe))
# phe_counts_SNP
# phe_SNP = unique(eqtl_coloc_SNP$phe)
# # levels
# tb_g = table(gene_SV %in% gene_SNP)
# tb_p = table(phe_SV %in% phe_SNP)
# gene_level  = round(tb_g[2] / length(gene_SV) * 100, 2)
# phe_level = round(tb_p[2] / length(phe_SV) * 100, 2)
# gene_level
# phe_level
# df_eqtl = data.frame(
#     gene_phe_uniq_counts = gene_phe_uniq_counts,
#     gene_counts = gene_counts,
#     phe_counts = phe_counts,
#     gene_phe_uniq_counts_SV = gene_phe_uniq_counts_SV,
#     gene_counts_SV = gene_counts_SV,
#     phe_counts_SV = phe_counts_SV,
#     gene_phe_uniq_counts_SNP = gene_phe_uniq_counts_SNP,
#     gene_counts_SNP = gene_counts_SNP,
#     phe_counts_SNP = phe_counts_SNP,
#     gene_level = gene_level,
#     phe_level = phe_level
# )
# df_eqtl
# write.table(df_eqtl, "03_eqtl_qtl_colocResultsig_stats.tsv", sep = "\t", quote = F, row.names = F)

head(eqtl_coloc)
gene_sv = unique(eqtl_coloc$gene[eqtl_coloc$var_type == "SV"])
gene_snv = unique(eqtl_coloc$gene[eqtl_coloc$var_type == "SNP"])
table(gene_sv %in% gene_snv)
length(unique(eqtl_coloc$gene))
round(table(gene_sv %in% gene_snv)[1] / length(unique(eqtl_coloc$gene)) * 100, 2)

gene_common = unique(eqtl_coloc$gene[eqtl_coloc$var_type == "SV"]) %>% intersect(unique(eqtl_coloc$gene[eqtl_coloc$var_type == "SNP"]))
counts = 0
for ( g in 1:length(gene_common)) {
    gene_sv_asso = eqtl_coloc$phe[eqtl_coloc$var_type == "SV" & eqtl_coloc$gene == gene_common[g]]
    gene_snv_asso = eqtl_coloc$phe[eqtl_coloc$var_type == "SNP" & eqtl_coloc$gene == gene_common[g]]
    if (any(! gene_sv_asso %in% gene_snv_asso)) {
        counts = counts + 1
    }
}
counts
round(counts / length(gene_common) * 100, 2)
