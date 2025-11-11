#!/usr/bin/env Rscript

# -------------
# FileName     : split_QTL_molphe.R
# Author       : Tian Wu fabien012390wt@163.com
# Modified By  : Tian Wu
# -------------
# Description  : Split QTL results by molecular phenotype for colocalization analysis
# -------------

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(parallel))
suppressPackageStartupMessages(library(GenomicRanges))

options(stringsAsFactors = FALSE)
options(scipen = 999)
setwd("/home/wtian/project/HZAU_cohort_meth/wdy_assist_rerun/11_coloc/04_rerun_QTL")

qtl_dir = "/home/wtian/project/HZAU_cohort_meth/wdy_assist_rerun/11_coloc/04_rerun_QTL/results"
qtl_subdir_list = list.dirs(qtl_dir, recursive = FALSE)
var_type_list = unique(str_split(basename(qtl_subdir_list), "-", simplify = TRUE)[,1])
molphe_name_list = unique(str_split(basename(qtl_subdir_list), "-", simplify = TRUE)[,2])

for (s in 1:length(qtl_subdir_list)) {
    qtl_subdir = qtl_subdir_list[s]
    qtl = fread(paste0(qtl_subdir, "/QTL_results/cis.txt"), header = TRUE, sep = "\t")
    if (!dir.exists(paste0(qtl_subdir, "/QTL_results/cis_each"))) {
        dir.create(paste0(qtl_subdir, "/QTL_results/cis_each"))
    }
    uniq_phe_list = unique(qtl$Phenotype)
    mclapply(1:length(uniq_phe_list), function(p) {
        phe = uniq_phe_list[p]
        outfile = paste0(qtl_subdir, "/QTL_results/cis_each/", phe, ".txt")
        if (!file.exists(outfile)) {
            qtl_phe = qtl[qtl$Phenotype == phe,]
            fwrite(qtl_phe, outfile, sep = "\t", quote = FALSE, row.names = FALSE)
        }
    }, mc.cores = 8)
}