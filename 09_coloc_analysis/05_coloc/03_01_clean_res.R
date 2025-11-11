#!/usr/bin/env Rscript

# -------------
# FileName     : 03_01_clean_res.R
# Author       : Tian Wu fabien012390wt@163.com
# Modified By  : Tian Wu
# -------------
# Description  : Clean and filter colocalization results for GWAS and eQTL analyses
# -------------

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(vroom))
suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(xQTLbiolinks))
suppressPackageStartupMessages(library(parallel))

setwd("/home/wtian/project/HZAU_cohort_meth/wdy_assist_rerun/11_coloc/05_coloc")

mc <- getOption("mc.cores", 128)

merge_tsv_files <- function(directory) {
  tsv_files <- list.files(directory, pattern = "\\.tsv$", full.names = TRUE, recursive = TRUE)
  
  data_list <- mclapply(tsv_files, function(file) {
    read.delim(file, header = TRUE, sep = "\t")
  }, mc.cores = mc)
  
  merged_data <- bind_rows(data_list)
  sorted_data <- merged_data %>%
    arrange(trait, var_type, molphe_name, phe)
  return(sorted_data)
}

colocResultAll <- merge_tsv_files("./01_02_gwas_coloc_res/")
colocResultsig <- filter(colocResultAll, PP.H4.abf > 0.75 & hypr_posterior > 0.5)
fwrite(colocResultsig, file = "03_gwas_qtl_colocResultsig.tsv", sep = "\t", quote = FALSE)

merge_tsv_files <- function(directory) {
  tsv_files <- list.files(directory, pattern = "\\.tsv$", full.names = TRUE, recursive = TRUE)
  
  data_list <- mclapply(tsv_files, function(file) {
    read.delim(file, header = TRUE, sep = "\t")
  }, mc.cores = mc)
  
  merged_data <- bind_rows(data_list)
  sorted_data <- merged_data %>%
    arrange(gene, var_type, molphe_name, phe)
  return(sorted_data)
}

colocResultAll <- merge_tsv_files("./02_02_eqtl_coloc_res/")
colocResultsig <- filter(colocResultAll, PP.H4.abf > 0.75 & hypr_posterior > 0.5)
fwrite(colocResultsig, file = "03_eqtl_qtl_colocResultsig.tsv", sep = "\t", quote = FALSE)
