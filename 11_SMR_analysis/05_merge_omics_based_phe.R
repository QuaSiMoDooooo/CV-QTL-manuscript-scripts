#!/usr/bin/env Rscript

# -------------
# FileName     : 05_merge_omics_based_phe.R
# Author       : Tian Wu fabien012390wt@163.com
# Modified By  : Tian Wu
# -------------
# Description  : Merge omics-based phenotype data from SMR analysis results
# -------------

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(vroom))
suppressPackageStartupMessages(library(stringr))
setwd("/home/wtian/project/HZAU_cohort_meth/wdy_assist_rerun/07_QTL_GWAS_SMR")

input_dir = "/home/wtian/project/HZAU_cohort_meth/wdy_assist_rerun/07_QTL_GWAS_SMR/04_smr_clean"
dirs <- list.dirs(input_dir, full.names = TRUE, recursive = FALSE)
cis_dirs = grep("cis", dirs, value = TRUE)
cis_name = basename(cis_dirs)
mols = str_replace(cis_name, "cis_", "")
cis_flt_path = paste0(cis_dirs, "/",cis_name, "_flt.tsv")

cis_flt = data.frame(matrix(ncol = 24, nrow = 0))
for (i in 1:length(cis_flt_path)) {
  tmp <- read_tsv(cis_flt_path[i])
  tmp$mol_phe <- mols[i]
  cis_flt <- rbind(cis_flt, tmp)
}

write.table(cis_flt, file = "05_smr_all_cis_flt.tsv", row.names = FALSE, sep = "\t", quote = FALSE)

traits_dir = "05_smr_per_trait"
if (!dir.exists(traits_dir)) {
  dir.create(traits_dir)
}
cis_flt <- cis_flt %>% select(mol_phe, everything())
for (trait in unique(cis_flt$trait)) {
  trait_data <- cis_flt[cis_flt$trait == trait, ]
  write.table(trait_data, file = paste0(traits_dir, "/", trait, "_cis_flt.tsv"), row.names = FALSE, sep = "\t", quote = FALSE)
}