#!/usr/bin/env Rscript

# -------------
# FileName     : 04_clean_smr_results.R
# Author       : Tian Wu fabien012390wt@163.com
# Modified By  : Tian Wu
# -------------
# Description  : Clean and filter SMR analysis results
# -------------

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(vroom))
suppressPackageStartupMessages(library(stringr))
setwd("/home/wtian/project/HZAU_cohort_meth/wdy_assist_rerun/07_QTL_GWAS_SMR")

input_dir = "/home/wtian/project/HZAU_cohort_meth/wdy_assist_rerun/07_QTL_GWAS_SMR/03_smr"
dirs <- list.dirs(input_dir, full.names = TRUE, recursive = FALSE)
dirs = dirs[grepl("cis", dirs)]

output_dir = "04_smr_clean"
for (dir in dirs) {
    qtl = basename(dir)
    output_subdir = paste0(output_dir, "/", qtl)
    output_all = paste0(output_subdir,"/",qtl,"_all_results.tsv")
    if (!dir.exists(output_subdir)) {
        dir.create(output_subdir, recursive = TRUE)
        files <- list.files(dir, pattern = "*.smr", full.names = TRUE)
        combined_data <- data.frame()
        for (file in files) {
            file_content <- read.table(file, sep = "\t", header = TRUE)
            trait_name <- sub("\\.smr$", "", basename(file))
            file_content <- file_content %>% mutate(trait = trait_name)
            combined_data <- bind_rows(combined_data, file_content)
        }
        fwrite(combined_data, output_all, quote = FALSE, sep = "\t", row.names = FALSE)

        combined_data$bonferroni_p <- p.adjust(combined_data$p_SMR, method = "bonferroni")
        combined_data <- combined_data %>% arrange(bonferroni_p)
        combined_data_flt <- combined_data %>% filter(bonferroni_p < 0.05, is.na(p_HEIDI) | p_HEIDI > 0.05)
        fwrite(combined_data_flt, paste0(output_subdir,"/",qtl,"_flt.tsv"), quote = FALSE, sep = "\t", row.names = FALSE)
    }
}